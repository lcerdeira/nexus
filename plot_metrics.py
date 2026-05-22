import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Dados de performance (extraídos da orquestração do HPC)
dados = pd.DataFrame({
    'Montador': ['Unicycler', 'Nexus', 'Unicycler', 'Nexus'],
    'Métrica': ['Contigs', 'Contigs', 'N50 (Mb)', 'N50 (Mb)'],
    'Valor': [15, 2, 1.20, 2.21]
})

# Configuração de estilo para publicação acadêmica
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Gráfico 1: Número de Contigs (Menor é melhor)
sns.barplot(
    data=dados[dados['Métrica'] == 'Contigs'],
    x='Montador', y='Valor', hue='Montador',
    palette=['#B0BEC5', '#2E86C1'], ax=ax[0]
)
ax[0].set_title('Fragmentação do Genoma\n(Menor é Melhor)', fontweight='bold')
ax[0].set_ylabel('Número de Contigs')

# Gráfico 2: N50 (Maior é melhor)
sns.barplot(
    data=dados[dados['Métrica'] == 'N50 (Mb)'],
    x='Montador', y='Valor', hue='Montador',
    palette=['#B0BEC5', '#2E86C1'], ax=ax[1]
)
ax[1].set_title('Contiguidade Estrutural (N50)\n(Maior é Melhor)', fontweight='bold')
ax[1].set_ylabel('Tamanho (Megabases)')

plt.suptitle('Comparação de Desempenho: Montagem Híbrida De Novo', fontsize=14, fontweight='bold', y=1.05)
sns.despine()
plt.tight_layout()

# Salvar em alta resolução
plt.savefig('performance_nexus_vs_unicycler.pdf', dpi=300, bbox_inches='tight')
print("Gráfico gerado com sucesso: 'performance_nexus_vs_unicycler.pdf'")