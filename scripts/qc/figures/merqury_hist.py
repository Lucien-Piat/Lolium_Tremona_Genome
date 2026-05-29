import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore

plt.style.use('default')
sns.set_theme(style="ticks", palette="muted")

df = pd.read_csv('results/data_circo/tremona.lmultiflorum.tremona.primary.spectra-cn.hist', 
                 sep=r'\s+', 
                 names=['Copies', 'kmer_multiplicity', 'Count'])

df['kmer_multiplicity'] = pd.to_numeric(df['kmer_multiplicity'], errors='coerce')
df['Count'] = pd.to_numeric(df['Count'], errors='coerce')
df = df.dropna()

pivot_df = df.pivot(index='kmer_multiplicity', columns='Copies', values='Count').fillna(0)

order = ['read-only', '1', '2']
existing_cols = [c for c in order if c in pivot_df.columns]
pivot_df = pivot_df[existing_cols]

merqury_colors = {
    'read-only': '#808080',  
    '1': '#e41a1c',          
    '2': '#377eb8'
}
colors = [merqury_colors[c] for c in existing_cols]

fig, ax = plt.subplots(figsize=(8, 6))

ax.stackplot(pivot_df.index, pivot_df.T, labels=pivot_df.columns, colors=colors, alpha=0.6, edgecolor='black', linewidth=0.5)

ax.set_xlim(0, 85) 
ax.set_ylim(0, 3.5e7) 
ax.set_xlabel('kmer_multiplicity', fontsize=12, fontweight='bold')
ax.set_ylabel('Count', fontsize=12, fontweight='bold')
ax.legend(title='k-mer', frameon=True)
sns.despine()
plt.tight_layout()
plt.savefig('clean_kmer_spectra.pdf', format='pdf', dpi=300)
plt.close()