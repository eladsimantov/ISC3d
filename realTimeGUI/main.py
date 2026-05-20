import time
from fastapi import FastAPI, WebSocket
from fastapi.responses import HTMLResponse
from processor import GaitProcessor
from sensors import mock_sensor_stream, udp_sensor_stream

# --- CONFIGURATION ---
USE_MOCK_DATA = True 
WINDOW_SIZE = 300 

app = FastAPI()
proc_L = GaitProcessor(WINDOW_SIZE)
proc_R = GaitProcessor(WINDOW_SIZE)

@app.get("/")
async def get_index():
    with open("index.html", "r") as f:
        return HTMLResponse(f.read())

@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    
    # Choose data source based on config
    stream = mock_sensor_stream() if USE_MOCK_DATA else udp_sensor_stream()
    
    async for left_data, right_data in stream:
        res_L = proc_L.update(*left_data)
        res_R = proc_R.update(*right_data)
        
        # Send data at ~20Hz to prevent browser lag
        if int(time.time() * 100) % 5 == 0: 
            await websocket.send_json({"left": res_L, "right": res_R})

if __name__ == "__main__":
    # Make sure to run in the iscRT folder 
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)