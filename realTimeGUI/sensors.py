import asyncio
import math
import random

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

async def udp_sensor_stream(ip="0.0.0.0", port=8080):
    """
    TODO: Insert your WitMotion API / UDP socket parsing here.
    Must yield exactly like the mock stream:
    yield (thigh_L, shank_L, foot_L), (thigh_R, shank_R, foot_R)
    """
    while True:
        await asyncio.sleep(0.01)
        yield (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)