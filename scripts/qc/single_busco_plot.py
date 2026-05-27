import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore

plt.style.use('default')
sns.set_theme(style="white")

# Load BUSCO table
df_busco = pd.read_csv('results/data_circo/full_table_busco_format.tsv', sep='\t', comment='#', 
                       names=['Busco_id', 'Status', 'Sequence', 'Gene_Start', 'Gene_End', 'Strand', 'Score', 'Length'])


df_busco = df_busco.drop_duplicates(subset=['Busco_id'])

counts = df_busco['Status'].value_counts(normalize=True) * 100

color_map = {
    'Complete': '#2ca02c',
    'Duplicated': '#ff7f0e',
    'Fragmented': '#d62728',
    'Missing': '#7f7f7f'
}

fig, ax = plt.subplots(figsize=(10, 1.2))

left_position = 0
for status in ['Complete', 'Duplicated', 'Fragmented', 'Missing']:
    if status in counts:
        pct = counts[status]
        ax.barh('BUSCO', pct, left=left_position, color=color_map[status], 
                label=f'{status} ({pct:.1f}%)', height=0.4, edgecolor='white', linewidth=1.5)
        left_position += pct

ax.set_xlim(0, 100)
ax.set_xlabel('Percentage of BUSCOs (%)', fontsize=12, labelpad=10)
ax.set_yticks([]) 

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.9), ncol=4, frameon=False)
sns.despine(left=True, top=True, right=True)

plt.savefig('clean_busco_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.close()