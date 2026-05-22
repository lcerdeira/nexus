#!/usr/bin/env python3
import json
import sys
import argparse

def evaluate_quality(fastp_json):
    """
    Simula a análise do Agente LLM verificando métricas rígidas de qualidade.
    """
    try:
        with open(fastp_json, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"[GATEKEEPER FATAL] Erro ao ler JSON: {e}")
        sys.exit(1)

    # Extraindo métricas pós-filtragem
    summary = data.get("summary", {})
    after_filtering = summary.get("after_filtering", {})
    
    q20_rate = after_filtering.get("q20_rate", 0)
    total_reads = after_filtering.get("total_reads", 0)
    
    print("==================================================")
    print("NEXUS LLM GATEKEEPER - AVALIAÇÃO DE QUALIDADE")
    print("==================================================")
    print(f"Taxa de Bases >= Q20: {q20_rate * 100:.2f}%")
    print(f"Total de Leituras Limpas: {total_reads}")

    # Regra 1: Phred Quality Cutoff (Q20 mínimo)
    # Se menos de 85% das bases atingirem Q20, reprovar a amostra
    if q20_rate < 0.85:
        print("\n❌ [VEREDITO: NO-GO] Qualidade sub-ótima detectada.")
        print("Motivo: A taxa de bases Q20 está abaixo do limite aceitável de 85%.")
        print("Ação: Abortando montagem para economizar recursos do HPC.")
        sys.exit(1) # Isso mata o processo no Nextflow

    # Regra 2: Quantidade mínima de reads (Ex: abortar se sobrar muito pouco)
    if total_reads < 100000:
        print("\n❌ [VEREDITO: NO-GO] Profundidade de sequenciamento insuficiente.")
        print("Motivo: Sobreviveram menos de 100k leituras após o trimming.")
        print("Ação: Abortando montagem.")
        sys.exit(1)

    print("\n✅ [VEREDITO: GO] Dados aprovados. Autorizando montagem...")
    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Nexus Quality Gatekeeper")
    parser.add_argument("--json", required=True, help="Caminho para o fastp.json")
    args = parser.parse_args()
    
    evaluate_quality(args.json)
