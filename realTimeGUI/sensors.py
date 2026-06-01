import asyncio
import math
import random
import json
import os
from readWIT import WITMotionPacketReader

# Fourier coefficients based on your fits
THIGH_P = {'a0': 0.0861, 'a1': -0.3439, 'b1': 0.1067, 'a2': -0.0555, 'b2': -0.0430, 'w': 2.004}
SHANK_P = {'a0': -0.2208, 'a1': -0.4512, 'b1': -0.2530, 'a2': 0.1945,  'b2': -0.1844, 'w': 1.659}
FOOT_P  = {'a0': 1.3709, 'a1': -0.3465, 'b1': -0.3631, 'a2': 0.2424,  'b2': -0.2449, 'a3': 0.0756, 'b3': 0.0689, 'w': 1.7673}

def eval_fourier(x_phase, p):
    """Evaluates the normalized Fourier fit for a given gait phase (0.0 to 1.0)."""
    z = (x_phase - 0.5) / 0.2916
    # add white noise to parameters to simulate variability
    p = {k: v * (1 + 0.005 * (2 * random.random() - 1)) for k, v in p.items()}
    xw = z * p['w']
    y = p['a0'] + p['a1']*math.cos(xw) + p['b1']*math.sin(xw) + p['a2']*math.cos(2*xw) + p['b2']*math.sin(2*xw)
    if 'a3' in p:
        y += p['a3']*math.cos(3*xw) + p['b3']*math.sin(3*xw)
    return y*180.0 / math.pi  # Convert to degrees

async def mock_sensor_stream():
    """Simulates real-time gait data using exact Fourier fits."""
    t = 0.0
    dt = 0.01             # 100Hz
    gait_cycle_sec = 1.0  # Duration of one full step cycle
    
    while True:
        await asyncio.sleep(dt)
        t += dt
        
        # Phase cycles from 0.0 to 1.0. Right leg is shifted by 0.5 (180 degrees)
        phase_L = (t % gait_cycle_sec) / gait_cycle_sec
        phase_R = (phase_L + 0.5) % 1.0
        
        # Left Leg
        t_L = eval_fourier(phase_L, THIGH_P)
        s_L = eval_fourier(phase_L, SHANK_P)
        f_L = eval_fourier(phase_L, FOOT_P)
        
        # Right Leg
        t_R = eval_fourier(phase_R, THIGH_P)
        s_R = eval_fourier(phase_R, SHANK_P)
        f_R = eval_fourier(phase_R, FOOT_P)
        
        yield (t_L, s_L, f_L), (t_R, s_R, f_R)

def get_elevation_angle(roll_deg):
    """
    Calculates the anatomical elevation angle of a leg segment.
    Since the IMUs are mounted such that rotation around the X-axis (Roll)
    represents segment movement, we add a 90-degree offset to calculate
    the elevation angle relative to the gravity vector.
    """
    return roll_deg + 90.0

class SensorProtocol(asyncio.DatagramProtocol):
    def __init__(self, queue):
        self.queue = queue
        self.readers = {}

    def connection_made(self, transport):
        pass

    def datagram_received(self, data, addr):
        ip, port = addr
        if ip not in self.readers:
            self.readers[ip] = WITMotionPacketReader()
        
        packets = self.readers[ip].feed(data)
        for packet in packets:
            if packet.get("name") == "angle":
                roll = packet["values"].get("roll_deg")
                if roll is not None:
                    elevation = get_elevation_angle(roll)
                    self.queue.put_nowait((ip, elevation))
            elif packet.get("name") == "wt55_record":
                raw_roll = packet["values"].get("angle_roll")
                if raw_roll is not None:
                    # Convert raw 16-bit to degrees
                    roll = raw_roll / 32768.0 * 180.0
                    elevation = get_elevation_angle(roll)
                    self.queue.put_nowait((ip, elevation))

async def udp_sensor_stream():
    """
    Listens for UDP packets from WitMotion sensors, decodes them, and yields the angles.
    Yields (thigh_L, shank_L, foot_L), (thigh_R, shank_R, foot_R)
    """
    config_path = os.path.join(os.path.dirname(__file__), "sensor_config.json")
    with open(config_path, "r") as f:
        config = json.load(f)
        
    udp_ip = config.get("udp", {}).get("ip", "0.0.0.0")
    udp_port = config.get("udp", {}).get("port", 1399)
    segment_by_ip = config.get("witmotion", {}).get("segment_by_ip", {})
    duplicate = config.get("pipeline", {}).get("duplicate_left_to_right", True)
    
    queue = asyncio.Queue()
    loop = asyncio.get_running_loop()
    transport, protocol = await loop.create_datagram_endpoint(
        lambda: SensorProtocol(queue),
        local_addr=(udp_ip, udp_port)
    )
    
    angles = {"thigh": 0.0, "shank": 0.0, "foot": 0.0}
    
    try:
        while True:
            ip, pitch = await queue.get()
            segment = segment_by_ip.get(ip)
            if segment in angles:
                angles[segment] = pitch
            
            left_data = (angles["thigh"], angles["shank"], angles["foot"])
            right_data = left_data if duplicate else (0.0, 0.0, 0.0)
            
            yield left_data, right_data
    finally:
        transport.close()