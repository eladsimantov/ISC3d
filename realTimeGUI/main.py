import time
# pyrefly: ignore [missing-import]
from fastapi import FastAPI, WebSocket
# pyrefly: ignore [missing-import]
from fastapi.responses import HTMLResponse
from processor import GaitProcessor, ExponentialGaitProcessor
from sensors import mock_sensor_stream, udp_sensor_stream
from sensors import eval_fourier, THIGH_P, SHANK_P, FOOT_P
# pyrefly: ignore [missing-import]
import numpy as np

# --- CONFIGURATION ---
USE_MOCK_DATA = False 
WINDOW_SIZE = 50 
calib_resolution = 100
calib_data = np.array([[eval_fourier(p, THIGH_P), eval_fourier(p, SHANK_P), eval_fourier(p, FOOT_P)] for p in np.linspace(0, 1, calib_resolution)])
app = FastAPI()
proc_L = ExponentialGaitProcessor(WINDOW_SIZE, calib_data, lambda_prior=0.5)
proc_R = ExponentialGaitProcessor(WINDOW_SIZE, calib_data, lambda_prior=0.5)
# proc_L = GaitProcessor(WINDOW_SIZE)
# proc_R = GaitProcessor(WINDOW_SIZE)

@app.get("/")
async def get_index():
    with open("index.html", "r") as f:
        return HTMLResponse(f.read())

@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    
    # Choose data source based on config
    stream = mock_sensor_stream() if USE_MOCK_DATA else udp_sensor_stream()
    
    last_send = 0
    async for left_data, right_data in stream:
        res_L = proc_L.update(*left_data)
        res_R = proc_R.update(*right_data)
        
        now = time.time()
        if now - last_send >= 0.03:  # Limit transmission rate to ~33Hz to prevent browser overload
            await websocket.send_json({"left": res_L, "right": res_R})
            last_send = now

if __name__ == "__main__":
    # Make sure to run in the iscRT folder 
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # pyrefly: ignore [missing-import]
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)