import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pycirclize import Circos

# ==========================================
# 1. Configuration and Input Files
# ==========================================
DIR = "results/before_after/"

FAI_BEFORE = DIR + "lmultiflorum.tremona_before.fa.fai"
FAI_AFTER = DIR + "lmultiflorum.tremona_after.fa.fai"

BUSCO_BEFORE = DIR + "lmultiflorum.tremona_full_table_busco_format_before.tsv"
BUSCO_AFTER = DIR + "lmultiflorum.tremona_full_table_busco_format_after.tsv"

SYNTENY_BEFORE = DIR + "synteny_ks_before.tsv"
SYNTENY_AFTER = DIR + "synteny_ks_after.tsv"

# Updated Styling
BUSCO_COLORS = {
    "Complete":   "#2e7d32", # Green
    "Duplicated": "#ef6c00", # Orange
    "Fragmented": "#fdd835", # Yellow
    "Missing":    "#9e9e9e", # Gray
}

SYNTENY_INTER = "#1976d240"  # Blue with hex alpha (40) for transparency
SYNTENY_INTRA = "#9c27b040"  # Purple with hex alpha (40)

# Matched to the image outer ring (matplotlib tab10 equivalents)
CHR_COLORS = {
    "chr1": "#1f77b4", # Blue
    "chr2": "#ff7f0e", # Orange
    "chr3": "#2ca02c", # Green
    "chr4": "#9467bd", # Purple
    "chr5": "#8c564b", # Brown
    "chr6": "#e377c2", # Pink
    "chr7": "#17becf", # Cyan
}

# ==========================================
# 2. Helper Functions
# ==========================================
def parse_fai(fai_path):
    df = pd.read_csv(fai_path, sep='\t', header=None, names=['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH'])
    df = df[df['NAME'].str.startswith('chr')].copy()
    return {row['NAME']: row['LENGTH'] for _, row in df.iterrows()}

def parse_busco(busco_path):
    cols = ["busco_id", "status", "sequence", "gene_start", "gene_end", "strand", "score", "length"]
    df = pd.read_csv(busco_path, sep='\t', comment='#', names=cols, on_bad_lines='skip')
    return df

def parse_synteny(synteny_path):
    return pd.read_csv(synteny_path, sep='\t')

def get_busco_percentages(busco_df):
    """Counts unique BUSCOs to avoid duplicating totals, returns percentages."""
    unique_buscos = busco_df.drop_duplicates(subset=['busco_id'])
    counts = unique_buscos['status'].value_counts()
    total = counts.sum()
    pcts = (counts / total * 100).to_dict() if total > 0 else {}
    return pcts, total

def build_circos_plot(ax, fai_path, synteny_path, busco_path):
    sectors_dict = parse_fai(fai_path)
    if not sectors_dict:
        return
    
    circos = Circos(sectors_dict, space=2)
    synteny_df = parse_synteny(synteny_path)
    busco_df = parse_busco(busco_path)
    
    target_statuses = ['Complete', 'Duplicated']
    track_buscos = busco_df[busco_df['status'].isin(target_statuses)].dropna(subset=['gene_start', 'gene_end'])
    
    for sector in circos.sectors:
        # Outer Track: Chromosome Ideogram (Widened: 90 to 100)
        chr_track = sector.add_track((90, 100))
        chr_color = CHR_COLORS.get(sector.name, "#cccccc") 
        chr_track.axis(fc=chr_color, ec="black", lw=0.8)
        chr_track.text(sector.name, r=105, size=9)
        
        # Inner Track: BUSCO Bands (Widened: 75 to 88)
        band_track = sector.add_track((75, 88))
        chr_genes = track_buscos[track_buscos['sequence'] == sector.name]
        
        for _, row in chr_genes.iterrows():
            color = BUSCO_COLORS.get(row['status'], "#757575")
            band_track.rect(row['gene_start'], row['gene_end'], fc=color, ec=color, lw=0.3)

    # Plot Synteny Links
    for _, row in synteny_df.iterrows():
        q_name, t_name = row['q_chr'], row['t_chr']
        if q_name in sectors_dict and t_name in sectors_dict:
            link_color = SYNTENY_INTRA if q_name == t_name else SYNTENY_INTER
            circos.link(
                (q_name, row['q_start'], row['q_end']),
                (t_name, row['t_start'], row['t_end']),
                color=link_color,
                direction=0
            )
            
    circos.plotfig(ax=ax)

# ==========================================
# 3. Figure Layout and Rendering
# ==========================================
def main():
    fig = plt.figure(figsize=(16, 14))
    
    ax_circos_before = fig.add_subplot(2, 2, 1, polar=True)
    ax_circos_after = fig.add_subplot(2, 2, 2, polar=True)
    ax_busco = fig.add_subplot(2, 2, 3)
    ax_hist = fig.add_subplot(2, 2, 4)
    
    # --- Panel a: Circos Before ---
    build_circos_plot(ax_circos_before, FAI_BEFORE, SYNTENY_BEFORE, BUSCO_BEFORE)
    ax_circos_before.set_title("a.", loc="left", fontweight="bold", fontsize=16, pad=20)
    
    # --- Panel b: Circos After ---
    build_circos_plot(ax_circos_after, FAI_AFTER, SYNTENY_AFTER, BUSCO_AFTER)
    ax_circos_after.set_title("b.", loc="left", fontweight="bold", fontsize=16, pad=20)
    
    # --- Panel c: BUSCO Percentages (Horizontal Stacked Bars) ---
    df_busco_before = parse_busco(BUSCO_BEFORE)
    df_busco_after = parse_busco(BUSCO_AFTER)
    
    pcts_before, _ = get_busco_percentages(df_busco_before)
    pcts_after, _ = get_busco_percentages(df_busco_after)
    
    statuses = ["Complete", "Duplicated", "Fragmented", "Missing"]
    
    y_pos = [1, 0]
    bar_width = 0.6
    
    left_before = 0
    left_after = 0
    
    for status in statuses:
        color = BUSCO_COLORS.get(status, "#757575")
        
        val_b = pcts_before.get(status, 0)
        val_a = pcts_after.get(status, 0)
        
        # Before bar
        if val_b > 0:
            ax_busco.barh(y_pos[0], val_b, left=left_before, height=bar_width, color=color, edgecolor='white')
            if status not in ["Fragmented", "Missing"] and val_b >= 2.0:
                ax_busco.text(left_before + (val_b / 2), y_pos[0], f"{val_b:.1f}%", 
                              va='center', ha='center', color='white', fontweight='bold', fontsize=11)
            left_before += val_b
            
        # After bar
        if val_a > 0:
            ax_busco.barh(y_pos[1], val_a, left=left_after, height=bar_width, color=color, edgecolor='white', label=status)
            if status not in ["Fragmented", "Missing"] and val_a >= 2.0:
                ax_busco.text(left_after + (val_a / 2), y_pos[1], f"{val_a:.1f}%", 
                              va='center', ha='center', color='white', fontweight='bold', fontsize=11)
            left_after += val_a

    ax_busco.set_yticks(y_pos)
    ax_busco.set_yticklabels(["Before Purging", "After Purging"], fontsize=12)
    ax_busco.set_xlabel("Percentage of Total BUSCO Groups (%)", fontsize=12)
    ax_busco.set_xlim(0, 100)
    ax_busco.spines['top'].set_visible(False)
    ax_busco.spines['right'].set_visible(False)
    ax_busco.set_title("c.", loc="left", fontweight="bold", fontsize=16, pad=20)
    ax_busco.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, frameon=False, fontsize=11)
    
    # --- Panel d: Synteny Block Age/Ks Distribution (Stacked Histogram) ---
    df_synteny_b = parse_synteny(SYNTENY_BEFORE)
    df_synteny_a = parse_synteny(SYNTENY_AFTER)
    
    x_col = 'age_Mya' if 'age_Mya' in df_synteny_b.columns else 'median_ks'
    xlabel_text = "Age (Mya)" if x_col == 'age_Mya' else "Median Ks"
    
    min_val = min(df_synteny_b[x_col].min(), df_synteny_a[x_col].min())
    max_val = max(df_synteny_b[x_col].max(), df_synteny_a[x_col].max())
    bins = np.linspace(min_val, max_val, 40)
    
    # Plot as a true stacked histogram
    ax_hist.hist([df_synteny_b[x_col], df_synteny_a[x_col]], 
                 bins=bins, 
                 stacked=True, 
                 color=['#757575', '#d81b60'], 
                 label=['Before Purging', 'After Purging'],
                 edgecolor='white',
                 linewidth=0.5)
    
    ax_hist.set_xlabel(xlabel_text, fontsize=12)
    ax_hist.set_ylabel("Number of Syntenic Blocks", fontsize=12)
    ax_hist.spines['top'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)
    ax_hist.set_title("d.", loc="left", fontweight="bold", fontsize=16, pad=20)
    ax_hist.legend(loc='upper right', frameon=False, fontsize=11)

    # Save as PDF and show
    plt.tight_layout()
    plt.savefig("deduplication_summary.pdf", format="pdf", dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()