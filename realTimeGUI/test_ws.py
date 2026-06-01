import asyncio
from sensors import udp_sensor_stream, mock_sensor_stream
import main

async def test():
    print(f"USE_MOCK_DATA = {main.USE_MOCK_DATA}")
    
    stream = mock_sensor_stream() if main.USE_MOCK_DATA else udp_sensor_stream()
    
    print("Waiting for data from stream...")
    try:
        async for left, right in stream:
            print(f"Received: {left}, {right}")
            break
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    asyncio.run(test())
