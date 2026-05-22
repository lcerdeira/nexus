import os
import sys
import pandas as pd
from google import genai
from dotenv import load_dotenv

# Carrega as variáveis de ambiente (onde está nossa GEMINI_API_KEY)
load_dotenv()

def parse_quast_report(tsv_path):
    """Lê o relatório TSV do QUAST e extrai as métricas de comparação."""
    try:
        df = pd.read_csv(tsv_path, sep='\t', index_col=0)
        
        metricas_chave = [
            "# contigs", 
            "Total length", 
            "N50", 
            "L50", 
            "N's per 100 kbp", 
            "mismatches per 100 kbp",
            "indels per 100 kbp"
        ]
        
        df_filtered = df.loc[df.index.intersection(metricas_chave)]
        return df_filtered.to_string()
    
    except Exception as e:
        print(f"Erro ao ler o arquivo QUAST: {e}")
        sys.exit(1)

def gerar_relatorio_llm(dados_quast):
    """Envia os dados extraídos para o LLM gerar um relatório executivo."""
    
    prompt = f"""
    Você é um bioinformata experiente avaliando o desempenho de algoritmos de montagem de genomas.
    Abaixo estão os resultados do relatório do QUAST comparando nosso novo montador híbrido ("nexus_assembly") 
    com o estado da arte ("unicycler_assembly").

    DADOS BRUTOS DO QUAST:
    {dados_quast}

    TAREFA:
    Escreva um relatório executivo de no máximo 3 parágrafos focando em três pontos:
    1. Contiguidade (resolução estrutural baseada no N50 e número de contigs).
    2. Fidelidade e Correção (taxa de mismatches e indels).
    3. Conclusão direta indicando qual montador apresentou a melhor estrutura genômica.

    Use linguagem profissional e direta.
    """

    print("Enviando dados para análise de Inteligência Artificial...\n")
    
    # Inicializa o cliente da nova biblioteca (ele puxa a chave do .env automaticamente)
    client = genai.Client()
    
    # Faz a chamada usando a sintaxe moderna do modelo Flash
    response = client.models.generate_content(
        model='gemini-flash-latest',
        contents=prompt
    )
    
    return response.text

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python app.py <caminho_para_report.tsv>")
        print("Rodando simulação com dados fictícios...\n")
        
        dados_simulados = """
                             nexus_assembly  unicycler_assembly
        # contigs                 3.0             15.0
        Total length        5200000.0        5100000.0
        N50                 4800000.0        1200000.0
        L50                       1.0              3.0
        N's per 100 kbp           0.0              0.0
        mismatches per 100 kbp    0.5             12.5
        indels per 100 kbp        0.1              5.0
        """
        relatorio = gerar_relatorio_llm(dados_simulados)
        print(relatorio)
    else:
        tsv_path = sys.argv[1]
        dados_quast = parse_quast_report(tsv_path)
        relatorio = gerar_relatorio_llm(dados_quast)
        
        with open("relatorio_final_nexus.txt", "w") as f:
            f.write(relatorio)
            
        print("Análise concluída com sucesso. Leia 'relatorio_final_nexus.txt'.")
        print("\n--- PRÉVIA DO RELATÓRIO ---\n")
        print(relatorio)