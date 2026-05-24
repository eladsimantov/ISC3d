import socket
from datetime import datetime
import string


PC_IP = "192.168.1.19"
PC_PORT = 1399

SENSOR_IPS = {"192.168.1.93", "192.168.1.94", "192.168.1.95"}

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

        while True:
            newline_index = self.buffer.find(b"\r\n")
            if newline_index < 0:
                break

            record = bytes(self.buffer[:newline_index])
            del self.buffer[:newline_index + 2]

            if record:
                packets.extend(self._decode_record(record))

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
            "acc_x", "acc_y", "acc_z",
            "gyro_x", "gyro_y", "gyro_z",
            "angle_roll", "angle_pitch", "angle_yaw",
            "mag_x", "mag_y", "mag_z",
            "quat_w", "quat_x", "quat_y", "quat_z",
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


sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.bind((PC_IP, PC_PORT))

reader = WITMotionPacketReader()

print(f"Listening for IMUs on {PC_IP}:{PC_PORT}...\n")

while True:
    data, addr = sock.recvfrom(4096)
    sensor_ip, sensor_port = addr

    if sensor_ip not in SENSOR_IPS:
        continue

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
    packets = reader.feed(data)

    if not packets:
        print(f"[{timestamp}] IMU {sensor_ip}:{sensor_port} received {len(data)} bytes")
        print(f"  hex: {data.hex(' ')}")
        print()
        continue

    for packet in packets:
        print_packet(timestamp, sensor_ip, sensor_port, packet)