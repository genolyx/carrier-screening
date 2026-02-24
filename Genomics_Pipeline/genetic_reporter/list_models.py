import os
from google import genai
from dotenv import load_dotenv

load_dotenv()
client = genai.Client(api_key=os.getenv("GEMINI_API_KEY"))

print("Checking available models for your API key...")

try:
    # This fetches the list of models available to YOUR specific key
    for model in client.models.list():
        print(f"Name: {model.name} | Capabilities: {model.supported_generation_methods}")
except Exception as e:
    print(f"❌ Error fetching models: {e}")
