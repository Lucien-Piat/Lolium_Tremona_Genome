import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.gridspec as gridspec # type: ignore
import numpy as np # type: ignore
from pycirclize import Circos # type: ignore
import matplotlib # type: ignore

# 0. Global Text & Size Styling
plt.style.use("seaborn-v0_8-white")
matplotlib.rcParams.update({
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "font.size": 25,
    "axes.titlesize": 24,
    "axes.labelsize": 25,
    "xtick.labelsize": 25,
    "ytick.labelsize": 25,
    "legend.fontsize": 25,
    "figure.titlesize": 26
})

# 1. Configuration and Input Files
DIR = "results/before_after/"

FAI_BEFORE = DIR + "lmultiflorum.tremona_before.fa.fai"
FAI_AFTER = DIR + "lmultiflorum.tremona_after.fa.fai"

BUSCO_BEFORE = DIR + "lmultiflorum.tremona_full_table_busco_format_before.tsv"
BUSCO_AFTER = DIR + "lmultiflorum.tremona_full_table_busco_format_after.tsv"

SYNTENY_BEFORE = DIR + "synteny_ks_before.tsv"
SYNTENY_AFTER = DIR + "synteny_ks_after.tsv"

BUSCO_COLORS = {
    "Complete":   "#2e7d32", 
    "Duplicated": "#ef6c00", 
    "Fragmented": "#fdd835", 
    "Missing":    "#9e9e9e", 
}

SYNTENY_INTER = "#1976d240"  
SYNTENY_INTRA = "#9c27b040"  

SYNTENY_INTER_SOLID = "#1976d2"
SYNTENY_INTRA_SOLID = "#9c27b0"

CHR_COLORS = {
    "chr1": "#1f77b4", "chr2": "#ff7f0e", "chr3": "#2ca02c", 
    "chr4": "#9467bd", "chr5": "#8c564b", "chr6": "#e377c2", "chr7": "#17becf",
}

# 2. Helper Functions
def parse_fai(fai_path):
    df = pd.read_csv(fai_path, sep='\t', header=None, names=['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH'])
    df = df[df['NAME'].str.startswith('chr')].copy()
    return {row['NAME']: row['LENGTH'] for _, row in df.iterrows()}

def parse_busco(busco_path):
    cols = ["busco_id", "status", "sequence", "gene_start", "gene_end", "strand", "score", "length"]
    return pd.read_csv(busco_path, sep='\t', comment='#', names=cols, on_bad_lines='skip')

def parse_synteny(synteny_path):
    df = pd.read_csv(synteny_path, sep='\t')
    df['link_type'] = np.where(df['q_chr'] == df['t_chr'], 'Intra', 'Inter')
    return df

def get_busco_percentages(busco_df):
    unique_buscos = busco_df.drop_duplicates(subset=['busco_id'])
    counts = unique_buscos['status'].value_counts()
    total = counts.sum()
    pcts = (counts / total * 100).to_dict() if total > 0 else {}
    return pcts

def build_circos_plot(ax, fai_path, synteny_path, busco_path, is_after=False, fai_before_path=None):
    sectors_dict = parse_fai(fai_path)
    if not sectors_dict: return
    
    if is_after and fai_before_path:
        sectors_before = parse_fai(fai_before_path)
        delta_bp = sum(sectors_before.values()) - sum(sectors_dict.values())
        if delta_bp > 0:
            sectors_dict['Purged'] = delta_bp
            
    circos = Circos(sectors_dict, space=2)
    synteny_df = parse_synteny(synteny_path)
    busco_df = parse_busco(busco_path)
    
    track_buscos = busco_df[busco_df['status'].isin(['Complete', 'Duplicated'])].dropna(subset=['gene_start', 'gene_end'])
    
    for sector in circos.sectors:
        chr_track = sector.add_track((92, 100))
        chr_color = CHR_COLORS.get(sector.name, "#757575") 
        chr_track.axis(fc=chr_color, ec="black", lw=1.0)
        
        chr_track.text(sector.name, r=112, size=25, fontweight="bold")
        
        if sector.name != 'Purged':
            band_track = sector.add_track((78, 90))
            chr_genes = track_buscos[track_buscos['sequence'] == sector.name]
            for _, row in chr_genes.iterrows():
                color = BUSCO_COLORS.get(row['status'], "#757575")
                band_track.rect(row['gene_start'], row['gene_end'], fc=color, ec=color, lw=0.1)

    for _, row in synteny_df.iterrows():
        q_name, t_name = row['q_chr'], row['t_chr']
        if q_name in sectors_dict and t_name in sectors_dict:
            circos.link(
                (q_name, row['q_start'], row['q_end']),
                (t_name, row['t_start'], row['t_end']),
                color=SYNTENY_INTRA if row['link_type'] == 'Intra' else SYNTENY_INTER,
                direction=0
            )
    circos.plotfig(ax=ax)

def plot_compact_busco(ax, busco_path, title_letter, show_legend=False):
    df_busco = parse_busco(busco_path)
    pcts = get_busco_percentages(df_busco)
    
    statuses = ["Complete", "Duplicated", "Fragmented", "Missing"]
    left = 0
    
    bar_thickness = 0.8 
    
    for status in statuses:
        val = pcts.get(status, 0)
        if val > 0:
            color = BUSCO_COLORS.get(status, "#757575")
            ax.barh(0, val, left=left, height=bar_thickness, color=color, edgecolor='white', lw=1.5, label=status)
            if status not in ["Fragmented"] and val >= 3.0:
                ax.text(left + (val / 2), 0, f"{val:.1f}%", va='center', ha='center', color='white', fontweight='bold', fontsize=22)
            left += val

    ax.set_yticks([])
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("BUSCOs (%)", fontweight="medium", labelpad=10)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    
    ax.set_title(title_letter, loc="left", fontweight="bold", pad=45)
    
    if show_legend:
        ax.legend(loc='upper center', bbox_to_anchor=(1.05, -1.2), ncol=4, frameon=False, handlelength=1.5)

def plot_intra_inter_hist(ax, df, title_letter, bins, x_col, xlabel_text, show_legend=False):
    df_intra = df[df['link_type'] == 'Intra'][x_col]
    df_inter = df[df['link_type'] == 'Inter'][x_col]
    
    ax.hist([df_intra, df_inter], bins=bins, stacked=True, 
            color=[SYNTENY_INTRA_SOLID, SYNTENY_INTER_SOLID], 
            label=['Intra-chromosomal', 'Inter-chromosomal'],
            edgecolor='white', linewidth=1.0)
    
    ax.set_xlabel(xlabel_text, fontweight="medium", labelpad=10)
    ax.set_ylabel("Blocks", fontweight="medium", labelpad=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title(title_letter, loc="left", fontweight="bold", pad=20)
    
    if show_legend:
        ax.legend(loc='upper right', frameon=False, handlelength=1.5)

# 3. Layout 
def main():
    fig = plt.figure(figsize=(22, 28))
    
    gs = gridspec.GridSpec(4, 2, height_ratios=[4.5, 0.6, 0.4, 1.5], hspace=0.35, wspace=0.1, left=0.02, right=0.98)
    
    ax_circos_before = fig.add_subplot(gs[0, 0], polar=True)
    ax_circos_after  = fig.add_subplot(gs[0, 1], polar=True)
    
    ax_busco_before  = fig.add_subplot(gs[1, 0])
    ax_busco_after   = fig.add_subplot(gs[1, 1])
    
    ax_hist_before   = fig.add_subplot(gs[3, 0])
    ax_hist_after    = fig.add_subplot(gs[3, 1])
    
    build_circos_plot(ax_circos_before, FAI_BEFORE, SYNTENY_BEFORE, BUSCO_BEFORE)
    ax_circos_before.set_title("E. Before ks-curation", loc="left", fontweight="bold", pad=30)
    
    build_circos_plot(ax_circos_after, FAI_AFTER, SYNTENY_AFTER, BUSCO_AFTER, is_after=True, fai_before_path=FAI_BEFORE)
    ax_circos_after.set_title("F. After ks-curation", loc="left", fontweight="bold", pad=30)
    
    plot_compact_busco(ax_busco_before, BUSCO_BEFORE, "G.", show_legend=True) 
    plot_compact_busco(ax_busco_after, BUSCO_AFTER, "H.", show_legend=False)
    
    df_b = parse_synteny(SYNTENY_BEFORE)
    df_a = parse_synteny(SYNTENY_AFTER)
    x_col = 'age_Mya' if 'age_Mya' in df_b.columns else 'median_ks'
    xlabel = "Age (Mya)" if x_col == 'age_Mya' else "Median Ks"
    
    min_val, max_val = min(df_b[x_col].min(), df_a[x_col].min()), max(df_b[x_col].max(), df_a[x_col].max())
    bins = np.linspace(min_val, max_val, 40)
    
    plot_intra_inter_hist(ax_hist_before, df_b, "I.", bins, x_col, xlabel, show_legend=True) 
    plot_intra_inter_hist(ax_hist_after, df_a, "J.", bins, x_col, xlabel, show_legend=False)
    
    # Sync Histogram Y-axes
    max_y = max(ax_hist_before.get_ylim()[1], ax_hist_after.get_ylim()[1])
    ax_hist_before.set_ylim(0, max_y)
    ax_hist_after.set_ylim(0, max_y)

    plt.savefig("deduplication_summary.pdf", format="pdf", dpi=300, bbox_inches='tight')
    print("Figure 'deduplication_summary.pdf' generated successfully.")

if __name__ == "__main__":
    main()