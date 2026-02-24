import os
from google import genai
from dotenv import load_dotenv

load_dotenv()
api_key = os.getenv("GEMINI_API_KEY")

if api_key:
    # Initialize the NEW client
    client = genai.Client(api_key=api_key)
    
    # In Dec 2025, the standard free model is 'gemini-2.5-flash-lite'
    # It is the most reliable for free tier users
    try:
        print("Connecting to Gemini 2.5 Flash-Lite...")
        response = client.models.generate_content(
            model="gemini-2.5-flash-lite",
            contents="Connection successful! Ready for genetic reporting."
        )
        print(f"🚀 AI Response: {response.text}")
        
    except Exception as e:
        print(f"❌ Connection Failed: {e}")
        print("\n--- Urgent Quota Check ---")
        print("Google drastically reduced free tier limits this month.")
        print("Check https://aistudio.google.com/ to ensure your 'Requests Per Day' isn't 0.")
else:
    print("❌ No API Key found.")
