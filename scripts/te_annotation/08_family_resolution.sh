#!/bin/bash
# scripts/te_annotation/08_family_resolution.sh


set -euo pipefail
OUTDIR=results/te_hite
RMOUT="${OUTDIR}/tremona_TE.out.gz" 
GENOME=reference_data/lmultiflorum.tremona.fa
FAI="${GENOME}.fai"
FAMTAB="${OUTDIR}/tremona_TE.family_table.tsv"
TMP="${OUTDIR}/.fam_raw.tmp"

GSIZE=$(awk '{s+=$2} END{print s}' "${FAI}")
read_out() { if [[ "${RMOUT}" == *.gz ]]; then zcat "${RMOUT}"; else cat "${RMOUT}"; fi; }

read_out | awk -v g="${GSIZE}" '
    function code(c){
        if(c=="LTR/Gypsy")return "RLG"; if(c=="LTR/Copia")return "RLC"
        if(c=="LTR/ERV")return "RLE";   if(c=="LTR/Pao")return "RLB"
        if(c=="LINE/L1")return "RIL";   if(c=="LINE/RTE-RTE")return "RIT"
        if(c=="SINE/tRNA")return "RST"
        if(c=="DNA/PIF-Harbinger")return "DTH"; if(c=="DNA/CMC-EnSpm")return "DTC"
        if(c=="DNA/MULE")return "DTM";  if(c=="DNA/hAT")return "DTA"
        if(c=="DNA/TcMar")return "DTT"; if(c=="RC/Helitron")return "DHH"
        if(c=="Unknown")return "RLX";   if(c=="Simple_repeat")return "SSR"
        if(c=="Low_complexity")return "LCX"; return "XXX"
    }
    NR>3 && NF>=11 {
        len=$7-$6+1
        fam=$10; sub(/-(int|INT|LTR|I)$/, "", fam)
        cls=$11
        bp[fam]+=len; n[fam]++
        if(cls!="Unknown") classOf[fam]=cls
        else if(!(fam in classOf)) classOf[fam]="Unknown"
    }
    END{ for(f in bp){ c=classOf[f]
        printf "%s\t%s\t%s\t%d\t%d\t%.4f\n", f, c, code(c), n[f], bp[f], 100*bp[f]/g } }
    ' > "${TMP}"

# label = code + rank within code (by bp)
sort -t$'\t' -k3,3 -k5,5rn "${TMP}" \
 | awk -F'\t' 'BEGIN{OFS="\t"}{ if($3!=prev){r=0;prev=$3} r++
       printf "%s_fam%02d\t%s\n", $3, r, $0 }' \
 | sort -t$'\t' -k6,6rn \
 | cat <(printf "label\tfamily_id\tclass\tcode\tfragments\tbp\tpct_genome\n") - \
 > "${FAMTAB}"
rm -f "${TMP}"
