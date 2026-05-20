# Real-Time ISC 3D (WitMotion IMU)

This runs a local web server to stream real-time elevation angles and compute sliding-window PCA for Inter-Segmental Coordination (ISC).

## 1. Prerequisites
Install the required Python packages in your environment:
```bash
pip install fastapi uvicorn numpy websockets
```

## 2. How to Run the Server
Open your terminal (Command Prompt or Anaconda Prompt), navigate to this folder, and run:
```bash
uvicorn main:app --host 0.0.0.0 --port 8000
```
*(Leave this terminal window open while you want the server to run).*

Alternatively, you can run the `main.py` script directly after activating your environment:
```bash
python main.py
```

## 3. View the Data on your PC
Open any web browser and go to:
[http://localhost:8000](http://localhost:8000)

## 4. View the Data on your Phone (Same WiFi Network)
To view the GUI on your phone, both devices must be on the same WiFi network.

1. **Find your PC's IP Address:**
   * Open a new Command Prompt on your PC and type `ipconfig`.
   * Look for the **IPv4 Address** under your WiFi or Ethernet adapter (e.g., `192.168.1.15`).
2. **Open on Phone:**
   * Open your phone's web browser.
   * Type `http://` followed by your PC's IP address and the port `:8000`. 
   * Example: `http://192.168.1.15:8000`

## 5. Switching to Real Sensors (Lab Mode)
Currently, the script runs a mock sine-wave generator for testing at home. 
When you are ready to use the real WitMotion UDP stream:
1. Open `main.py`.
2. Change line 13 to: `USE_MOCK_DATA = False`.
3. Restart the server.

## 6. Editing the sensors and processor modules
The `sensors.py` and `processor.py` files are modularized for better organization.
- `sensors.py` contains functions for both mock data generation and real UDP streaming.
- `processor.py` contains the logic for processing the sensor data, including PCA computation.
You can edit these files to modify the data processing or sensor input methods without affecting the main server code.
