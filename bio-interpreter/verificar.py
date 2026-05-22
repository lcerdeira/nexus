from google import genai
from dotenv import load_dotenv

load_dotenv()
client = genai.Client()

print("Buscando modelos disponíveis para a sua chave...")
try:
    for model in client.models.list():
        print(model.name)
except Exception as e:
    print(f"Erro ao buscar modelos: {e}")
