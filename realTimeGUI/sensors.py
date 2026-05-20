import asyncio
import math

async def mock_sensor_stream():
    """Simulates real-time gait data."""
    t = 0
    while True:
        await asyncio.sleep(0.01) # 100Hz
        t += 0.05
        # Left Leg
        t_L, s_L, f_L = 40*math.sin(t), 45*math.sin(t - 1.5), 20*math.sin(t - 2.5)
        # Right Leg
        t_R, s_R, f_R = 40*math.sin(t + math.pi), 45*math.sin(t + math.pi - 1.5), 20*math.sin(t + math.pi - 2.5)
        
        yield (t_L, s_L, f_L), (t_R, s_R, f_R)

async def udp_sensor_stream(ip="0.0.0.0", port=8080):
    """
    TODO: Insert your WitMotion API / UDP socket parsing here.
    Must yield exactly like the mock stream:
    yield (thigh_L, shank_L, foot_L), (thigh_R, shank_R, foot_R)
    """
    while True:
        await asyncio.sleep(0.01)
        # Dummy yield until implemented
        yield (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)