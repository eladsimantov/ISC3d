import socket
from datetime import datetime
import string
import json
import os

FRAME_HEADER = 0x55
FRAME_LENGTH = 11
WT55_PREFIX_LENGTH = 12
WT55_TIMESTAMP_LENGTH = 6
WT55_TRAILER_LENGTH = 2


class WITMotionPacketReader:
    def __init__(self):
        self.buffer = bytearray()

    def feed(self, data):
        self.buffer.extend(data)
        packets = []

        while len(self.buffer) > 0:
            wt55_idx = self.buffer.find(b"WT55")
            bin_idx = self.buffer.find(b"\x55")

            first_idx = -1
            is_wt55 = False

            if wt55_idx != -1 and bin_idx != -1:
                if wt55_idx < bin_idx:
                    first_idx = wt55_idx
                    is_wt55 = True
                else:
                    first_idx = bin_idx
                    is_wt55 = False
            elif wt55_idx != -1:
                first_idx = wt55_idx
                is_wt55 = True
            elif bin_idx != -1:
                first_idx = bin_idx
                is_wt55 = False

            if first_idx == -1:
                # No headers found. Clear the buffer except maybe the last 3 bytes
                # in case a partial "WT55" header is being received.
                keep = min(3, len(self.buffer))
                if keep > 0:
                    del self.buffer[:-keep]
                else:
                    self.buffer.clear()
                break

            # If there is garbage before the first header, discard it
            if first_idx > 0:
                del self.buffer[:first_idx]
                continue

            # Now, the buffer starts with the header (at index 0)
            if is_wt55:
                newline_idx = self.buffer.find(b"\r\n")
                if newline_idx == -1:
                    # Partial record. Wait for more data.
                    if len(self.buffer) > 4096:
                        # Safety fallback to prevent memory leaks if \r\n is missing
                        del self.buffer[:4]
                    break

                record = bytes(self.buffer[:newline_idx])
                del self.buffer[:newline_idx + 2]

                if record:
                    packets.extend(self._decode_record(record))
            else:
                # Binary frame: starts with 0x55. Must be at least 11 bytes.
                if len(self.buffer) < FRAME_LENGTH:
                    break

                frame = bytes(self.buffer[:FRAME_LENGTH])
                if self._is_valid_frame(frame):
                    decoded = self._decode_frame(frame)
                    if decoded is not None:
                        packets.append(decoded)
                    del self.buffer[:FRAME_LENGTH]
                else:
                    # Invalid frame. Discard the 0x55 byte to search for the next one.
                    del self.buffer[:1]

        return packets

    def _decode_record(self, record):
        packets = []
        packets.append(
            {
                "kind": "raw",
                "name": "record",
                "type": None,
                "values": {
                    "byte_count": len(record),
                    "hex": record.hex(" "),
                    "ascii": self._extract_ascii(record),
                },
            }
        )

        wt55_packet = self._decode_wt55_record(record)
        if wt55_packet is not None:
            packets.append(wt55_packet)

        packets.extend(self._extract_standard_frames(record))
        return packets

    def _decode_wt55_record(self, record):
        if not record.startswith(b"WT55"):
            return None

        minimum_length = WT55_PREFIX_LENGTH + WT55_TIMESTAMP_LENGTH + WT55_TRAILER_LENGTH
        if len(record) < minimum_length:
            return None

        prefix = record[:WT55_PREFIX_LENGTH].decode("ascii", errors="replace")
        timestamp_bytes = record[WT55_PREFIX_LENGTH:WT55_PREFIX_LENGTH + WT55_TIMESTAMP_LENGTH]
        trailer = record[-WT55_TRAILER_LENGTH:] if len(record) >= WT55_TRAILER_LENGTH else b""
        payload = record[WT55_PREFIX_LENGTH + WT55_TIMESTAMP_LENGTH:len(record) - WT55_TRAILER_LENGTH]

        if len(payload) % 2 != 0:
            return {
                "kind": "wt55",
                "name": "wt55_record",
                "type": None,
                "values": {
                    "device": prefix,
                    "timestamp": self._decode_timestamp(timestamp_bytes),
                    "payload_hex": payload.hex(" "),
                    "trailer_hex": trailer.hex(" "),
                },
            }

        words = [self._to_signed_16(payload[i], payload[i + 1]) for i in range(0, len(payload), 2)]
        values = {
            "device": prefix,
            "timestamp": self._decode_timestamp(timestamp_bytes),
            "trailer_hex": trailer.hex(" "),
            "raw_payload_hex": payload.hex(" "),
        }

        field_names = [
            "ms_counter",
            "acc_x", "acc_y", "acc_z",
            "gyro_x", "gyro_y", "gyro_z",
            "mag_x", "mag_y", "mag_z",
            "angle_roll", "angle_pitch", "angle_yaw",
            "quat_w", "quat_x", "quat_y",
        ]

        if len(words) == len(field_names):
            values.update(dict(zip(field_names, words)))
        else:
            for index, word in enumerate(words):
                values[f"value_{index}"] = word

        return {
            "kind": "wt55",
            "name": "wt55_record",
            "type": None,
            "values": values,
        }

    def _extract_standard_frames(self, data):
        packets = []
        index = 0

        while index <= len(data) - FRAME_LENGTH:
            if data[index] != FRAME_HEADER:
                index += 1
                continue

            frame = data[index:index + FRAME_LENGTH]
            if not self._is_valid_frame(frame):
                index += 1
                continue

            decoded = self._decode_frame(frame)
            if decoded is not None:
                packets.append(decoded)
            index += FRAME_LENGTH

        return packets

    @staticmethod
    def _decode_timestamp(timestamp_bytes):
        if len(timestamp_bytes) != WT55_TIMESTAMP_LENGTH:
            return None

        year, month, day, hour, minute, second = timestamp_bytes
        return f"20{year:02d}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:02d}"

    @staticmethod
    def _extract_ascii(data):
        text = ''.join(chr(byte) if chr(byte) in string.printable and byte not in (10, 13) else ' ' for byte in data)
        text = ' '.join(text.split())
        return text if text else None

    def _is_valid_frame(self, frame):
        if len(frame) != FRAME_LENGTH:
            return False
        if frame[0] != FRAME_HEADER:
            return False
        checksum = sum(frame[:10]) & 0xFF
        return checksum == frame[10]

    def _decode_frame(self, frame):
        frame_type = frame[1]
        ax = self._to_signed_16(frame[2], frame[3])
        ay = self._to_signed_16(frame[4], frame[5])
        az = self._to_signed_16(frame[6], frame[7])
        extra = self._to_signed_16(frame[8], frame[9])

        if frame_type == 0x51:
            return {
                "name": "acceleration",
                "type": frame_type,
                "values": {
                    "ax_g": ax / 32768.0 * 16.0,
                    "ay_g": ay / 32768.0 * 16.0,
                    "az_g": az / 32768.0 * 16.0,
                    "temperature_c": extra / 100.0,
                },
            }

        if frame_type == 0x52:
            return {
                "name": "angular_velocity",
                "type": frame_type,
                "values": {
                    "gx_dps": ax / 32768.0 * 2000.0,
                    "gy_dps": ay / 32768.0 * 2000.0,
                    "gz_dps": az / 32768.0 * 2000.0,
                    "temperature_c": extra / 100.0,
                },
            }

        if frame_type == 0x53:
            return {
                "name": "angle",
                "type": frame_type,
                "values": {
                    "roll_deg": ax / 32768.0 * 180.0,
                    "pitch_deg": ay / 32768.0 * 180.0,
                    "yaw_deg": az / 32768.0 * 180.0,
                    "temperature_c": extra / 100.0,
                },
            }

        if frame_type == 0x54:
            return {
                "name": "magnetometer",
                "type": frame_type,
                "values": {
                    "mx": ax,
                    "my": ay,
                    "mz": az,
                    "temperature_c": extra / 100.0,
                },
            }

        return {
            "name": f"unknown_0x{frame_type:02X}",
            "type": frame_type,
            "values": {
                "raw_1": ax,
                "raw_2": ay,
                "raw_3": az,
                "raw_4": extra,
            },
        }

    @staticmethod
    def _to_signed_16(low, high):
        value = low | (high << 8)
        if value & 0x8000:
            value -= 0x10000
        return value


def print_packet(timestamp, sensor_ip, sensor_port, packet):
    if packet["kind"] == "raw":
        print(f"[{timestamp}] IMU {sensor_ip}:{sensor_port} raw record")
        print(f"  bytes: {packet['values']['byte_count']}")
        print(f"  hex:   {packet['values']['hex']}")
        if packet["values"]["ascii"]:
            print(f"  ascii: {packet['values']['ascii']}")
        print()
        return

    if packet["kind"] == "wt55":
        print(f"[{timestamp}] IMU {sensor_ip}:{sensor_port} WT55 record")
        print(f"  device: {packet['values'].get('device')}")
        if packet["values"].get("timestamp"):
            print(f"  sensor_time: {packet['values']['timestamp']}")
        print(f"  trailer: {packet['values'].get('trailer_hex')}")
        if packet["values"].get("raw_payload_hex"):
            print(f"  payload: {packet['values']['raw_payload_hex']}")
        for key in [
            "acc_x", "acc_y", "acc_z",
            "gyro_x", "gyro_y", "gyro_z",
            "angle_roll", "angle_pitch", "angle_yaw",
            "mag_x", "mag_y", "mag_z",
            "quat_w", "quat_x", "quat_y", "quat_z",
        ]:
            if key in packet["values"]:
                print(f"  {key}: {packet['values'][key]}")
        for key, value in packet["values"].items():
            if key.startswith("value_"):
                print(f"  {key}: {value}")
        print()
        return

    print(f"[{timestamp}] IMU {sensor_ip}:{sensor_port} {packet['name']} (0x{packet['type']:02X})")
    for key, value in packet["values"].items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value}")
    print()


if __name__ == "__main__":
    import sys
    import time
    config_path = os.path.join(os.path.dirname(__file__), "sensor_config.json")
    with open(config_path, "r") as f:
        config = json.load(f)

    pc_ip = config.get("udp", {}).get("ip", "0.0.0.0")
    pc_port = config.get("udp", {}).get("port", 1399)
    segment_by_ip = config.get("witmotion", {}).get("segment_by_ip", {})
    configured_ips = set(segment_by_ip.keys())

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        sock.bind((pc_ip, pc_port))
    except OSError as e:
        print(f"\n[ERROR] Could not bind to port {pc_port} on {pc_ip}.")
        print("This usually means another program is already using this port (e.g., main.py or another instance of readWIT.py).")
        print("Please make sure all other GUI/sensor scripts are closed and try again.\n")
        sys.exit(1)

    readers = {}

    # Initialize dashboard dictionary based on configured sensors
    latest = {}
    for ip, segment in segment_by_ip.items():
        latest[segment] = {
            "acc": [0.0, 0.0, 0.0],
            "gyro": [0.0, 0.0, 0.0],
            "angle": [0.0, 0.0, 0.0],
            "ip": ip,
            "last_seen": None
        }

    # Enable ANSI escape codes on Windows console
    os.system('')

    print(f"Listening for IMUs on {pc_ip}:{pc_port}...")
    print(f"Configured IMU IPs: {configured_ips}")
    print("Starting Live IMU Dashboard. Press Ctrl+C to exit.\n")

    last_line_count = 0
    last_print_time = 0
    unknown_ips = set()

    while True:
        try:
            data, addr = sock.recvfrom(4096)
        except KeyboardInterrupt:
            print("\nExiting dashboard.")
            break

        sensor_ip, sensor_port = addr
        now = time.time()

        if sensor_ip not in configured_ips:
            if sensor_ip not in unknown_ips:
                unknown_ips.add(sensor_ip)
                # Move cursor down to not overwrite history, print warning, and reset dashboard position
                sys.stdout.write("\n" * last_line_count)
                print(f"[WARNING] UNKNOWN SENSOR DISCOVERED: {sensor_ip}:{sensor_port} (Edit sensor_config.json to add it!)")
                print()
                last_line_count = 0
            continue

        segment = segment_by_ip[sensor_ip]

        if sensor_ip not in readers:
            readers[sensor_ip] = WITMotionPacketReader()

        packets = readers[sensor_ip].feed(data)

        for packet in packets:
            vals = packet.get("values", {})
            if packet.get("kind") == "wt55":
                if "acc_x" in vals:
                    latest[segment]["acc"] = [
                        vals.get("acc_x", 0) / 32768.0 * 16.0,
                        vals.get("acc_y", 0) / 32768.0 * 16.0,
                        vals.get("acc_z", 0) / 32768.0 * 16.0
                    ]
                if "gyro_x" in vals:
                    latest[segment]["gyro"] = [
                        vals.get("gyro_x", 0) / 32768.0 * 2000.0,
                        vals.get("gyro_y", 0) / 32768.0 * 2000.0,
                        vals.get("gyro_z", 0) / 32768.0 * 2000.0
                    ]
                if "angle_roll" in vals:
                    latest[segment]["angle"] = [
                        vals.get("angle_roll", 0) / 32768.0 * 180.0,
                        vals.get("angle_pitch", 0) / 32768.0 * 180.0,
                        vals.get("angle_yaw", 0) / 32768.0 * 180.0
                    ]
                latest[segment]["last_seen"] = now
            else:
                if packet.get("name") == "acceleration":
                    latest[segment]["acc"] = [vals.get("ax_g", 0.0), vals.get("ay_g", 0.0), vals.get("az_g", 0.0)]
                elif packet.get("name") == "angular_velocity":
                    latest[segment]["gyro"] = [vals.get("gx_dps", 0.0), vals.get("gy_dps", 0.0), vals.get("gz_dps", 0.0)]
                elif packet.get("name") == "angle":
                    latest[segment]["angle"] = [vals.get("roll_deg", 0.0), vals.get("pitch_deg", 0.0), vals.get("yaw_deg", 0.0)]
                latest[segment]["last_seen"] = now

        # Throttle dashboard updates to 10Hz to prevent terminal lag
        if now - last_print_time >= 0.1:
            active_segments = []
            for segment, sdata in latest.items():
                is_active = sdata["last_seen"] and (now - sdata["last_seen"] < 2.0)
                if is_active:
                    active_segments.append((segment, sdata))

            # Move cursor up by the number of lines printed last time
            if last_line_count > 0:
                sys.stdout.write(f"\033[{last_line_count}A")

            # Print active segments
            for segment, sdata in active_segments:
                acc_str = f"Acc: [{sdata['acc'][0]:6.2f}, {sdata['acc'][1]:6.2f}, {sdata['acc'][2]:6.2f}] G"
                gyro_str = f"Gyro: [{sdata['gyro'][0]:7.1f}, {sdata['gyro'][1]:7.1f}, {sdata['gyro'][2]:7.1f}] dps"
                angle_str = f"Angle: [{sdata['angle'][0]:7.1f}, {sdata['angle'][1]:7.1f}, {sdata['angle'][2]:7.1f}] deg"
                print(f"{segment:<6} ({sdata['ip']}): \033[92mACTIVE\033[0m  | {acc_str} | {gyro_str} | {angle_str}\033[K")

            # Print blank lines for any deactivated segments to clear them
            excess = last_line_count - len(active_segments)
            for _ in range(excess):
                print("\033[K")

            # Move cursor back up over the blank lines so the next run starts at the active count
            if excess > 0:
                sys.stdout.write(f"\033[{excess}A")

            last_line_count = len(active_segments)
            sys.stdout.flush()
            last_print_time = now