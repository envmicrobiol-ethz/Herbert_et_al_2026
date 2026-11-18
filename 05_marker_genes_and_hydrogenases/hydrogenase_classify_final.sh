#!/bin/bash
#SBATCH --job-name=hyd_classify_final
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=hyd_classify_final_%j.out
#SBATCH --error=hyd_classify_final_%j.err

# ============================================================
# Definitive hydrogenase classification for 159 MAGs + assemblies
# Per Søndergaard et al. 2016 / HydDB approach:
#   1. HMMsearch (HMMer v3.4) against HydDB HMMs, e-value <1e-5
#   2. DIAMOND blastp against HydDB sequences, max-target-seqs 1
#   3. Identity cutoffs:
#      NiFe Group 4: >50%
#      NiFe other:   >30%
#      FeFe:         >45%
#      Fe:           >50%
#
# Output protein_id format matches TranscriptM gene_id for joining
# ============================================================

set -eo pipefail

source [conda_base]/etc/profile.d/conda.sh
module load stack/.2024-05-silent gcc/13.2.0 prodigal/2.6.3
module load stack/.2024-06-silent gcc/12.2.0 openmpi/4.1.6
module load hmmer/3.4-pronkqe

MAG_FAA="21_coassembly/MAGs/hyd_coxl_classification/prodigal_prefixed/all_159_MAGs_prefixed.faa"
ASSEMBLY_DIR="11_metabolic_assemblies/input_genomes"
HYD_HMM="HMMs/hyd/HydDB_all_MM2022.hmm"
HYD_DMND="HMMs/HydDB_all_hydrogenases.dmnd"
OUTDIR="27_hydrogenase_coxl_final"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}/hmmsearch" "${OUTDIR}/diamond" "${OUTDIR}/results"

conda activate diamond_2.1.10

# ════════════════════════════════════════════════════════════
# Function: HMM search + DIAMOND classify hydrogenases
# ════════════════════════════════════════════════════════════
classify_hydrogenases() {
    local FAA="$1"
    local LABEL="$2"
    local DATA_TYPE="$3"
    local HMM_OUT="${OUTDIR}/hmmsearch/${LABEL}_hyd_tblout.txt"
    local HMM_SUMMARY="${OUTDIR}/hmmsearch/${LABEL}_hyd_hmm_hits.tsv"
    local DIAMOND_OUT="${OUTDIR}/diamond/${LABEL}_hyd_diamond.tsv"
    local RESULTS="${OUTDIR}/results/${LABEL}_hydrogenase_classification.tsv"

    if [[ ! -f "${FAA}" ]]; then
        echo "  [SKIP] ${LABEL} — FAA not found"
        return
    fi
    if [[ -f "${RESULTS}" ]]; then
        echo "  [SKIP] ${LABEL} — already done"
        return
    fi

    # ── Step 1: HMMsearch ────────────────────────────────────
    echo "  [$(date '+%H:%M:%S')] HMM search: ${LABEL}..."
    hmmsearch \
        --tblout "${HMM_OUT}" \
        -E 1e-5 \
        --cpu "${THREADS}" \
        "${HYD_HMM}" "${FAA}"

    echo -e "protein_id\thmm_profile\thmm_evalue\thmm_score" > "${HMM_SUMMARY}"
    grep -v "^#" "${HMM_OUT}" | \
        awk 'BEGIN{OFS="\t"} {print $1,$3,$5,$6}' >> "${HMM_SUMMARY}"

    N_HMM=$(tail -n +2 "${HMM_SUMMARY}" | wc -l)
    echo "  HMM hits: ${N_HMM}"

    if [[ "${N_HMM}" -eq 0 ]]; then
        echo -e "data_type\tsource_label\tprotein_id\tcontig_id\thmm_profile\thmm_evalue\thmm_score\tsubject_id\tpct_identity\talign_length\tevalue\tbitscore\thydrogenase_group\tclass_flag" \
            > "${RESULTS}"
        echo "  No HMM hits for ${LABEL}"
        return
    fi

    # Extract hydrogenase protein sequences
    tail -n +2 "${HMM_SUMMARY}" | awk '{print $1}' | sort -u \
        > "${OUTDIR}/hmmsearch/${LABEL}_hyd_ids.txt"

    awk 'NR==FNR{ids[$1]=1; next}
         /^>/{id=substr($0,2); gsub(/ .*/,"",id); found=(id in ids)}
         found{print}' \
        "${OUTDIR}/hmmsearch/${LABEL}_hyd_ids.txt" "${FAA}" \
        > "${OUTDIR}/hmmsearch/${LABEL}_hyd_seqs.faa"

    N_SEQS=$(grep -c "^>" "${OUTDIR}/hmmsearch/${LABEL}_hyd_seqs.faa" || true)
    echo "  Extracted sequences: ${N_SEQS}"

    # ── Step 2: DIAMOND classification ───────────────────────
    echo "  [$(date '+%H:%M:%S')] DIAMOND classification: ${LABEL}..."
    diamond blastp \
        -q "${OUTDIR}/hmmsearch/${LABEL}_hyd_seqs.faa" \
        -d "${HYD_DMND}" \
        -o "${DIAMOND_OUT}" \
        --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --header simple \
        --threads "${THREADS}" \
        --quiet

    # ── Step 3: Apply identity cutoffs + annotate ─────────────
    # Output includes all fields needed for TranscriptM join
    echo -e "data_type\tsource_label\tprotein_id\tcontig_id\thmm_profile\thmm_evalue\thmm_score\tsubject_id\tpct_identity\talign_length\tevalue\tbitscore\thydrogenase_group\tclass_flag" \
        > "${RESULTS}"

    # Load HMM info for join
    declare -A HMM_PROF HMM_EVAL HMM_SCORE
    while IFS=$'\t' read -r pid prof he hs; do
        HMM_PROF["${pid}"]="${prof}"
        HMM_EVAL["${pid}"]="${he}"
        HMM_SCORE["${pid}"]="${hs}"
    done < <(tail -n +2 "${HMM_SUMMARY}")

    awk -v dtype="${DATA_TYPE}" \
        -v label="${LABEL}" \
        'NR>1{
            # Parse hydrogenase group from subject ID
            # Format: accession|organism|[NiFe]_Group_1h or [FeFe] etc
            split($2, arr, /\|/)
            group = (length(arr) >= 3) ? arr[3] : arr[length(arr)]
            pct = $3 + 0
            flag = "PASS"
            # Apply identity cutoffs per class
            if (group ~ /Group_4/)               { if (pct < 50) flag = "LOW_ID" }
            else if ($2 ~ /\[NiFe\]|NiFe/)       { if (pct < 30) flag = "LOW_ID" }
            else if ($2 ~ /\[FeFe\]|FeFe/)       { if (pct < 45) flag = "LOW_ID" }
            else if ($2 ~ /\[Fe\]|Fe_only/)      { if (pct < 50) flag = "LOW_ID" }
            # contig_id = protein_id with last _N stripped
            contig = $1; sub(/_[0-9]+$/, "", contig)
            print dtype"\t"label"\t"$1"\t"contig"\t.\t.\t.\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"group"\t"flag
        }' "${DIAMOND_OUT}" >> "${RESULTS}"

    # Fill in HMM info via python
    source [conda_base]/etc/profile.d/conda.sh
conda activate transcriptm_v0.3.1
python3 << PYEOF
import pandas as pd, re
res = pd.read_csv("${RESULTS}", sep="\t")
hmm = pd.read_csv("${HMM_SUMMARY}", sep="\t")
hmm_dict = hmm.set_index("protein_id")[["hmm_profile","hmm_evalue","hmm_score"]].to_dict("index")
for col in ["hmm_profile","hmm_evalue","hmm_score"]:
    res[col] = res["protein_id"].map(lambda x: hmm_dict.get(x, {}).get(col, "."))
res.to_csv("${RESULTS}", sep="\t", index=False)
PYEOF
conda activate diamond_2.1.10

    N_PASS=$(grep -c "PASS" "${RESULTS}" 2>/dev/null || echo 0)
    N_TOTAL=$(tail -n +2 "${RESULTS}" | wc -l)
    echo "  ${LABEL}: ${N_TOTAL} classified, ${N_PASS} PASS"
    echo "  Group distribution (PASS):"
    grep "PASS" "${RESULTS}" | awk -F'\t' '{print $13}' | sort | uniq -c | sort -rn | \
        awk '{print "    "$0}'
}

# ════════════════════════════════════════════════════════════
# 1. 159 MAGs
# ════════════════════════════════════════════════════════════
echo "[$(date)] === 159 MAGs ==="
classify_hydrogenases "${MAG_FAA}" "MAGs_159" "MAG"

# ════════════════════════════════════════════════════════════
# 2. Assemblies
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Assemblies ==="
for ASM in WIL1_before WIL2_before WIL2_after WILB_D1 WILB_D3 WILB_D10 WILD_D0; do
    classify_hydrogenases "${ASSEMBLY_DIR}/${ASM}.faa" "${ASM}" "assembly"
done

# ════════════════════════════════════════════════════════════
# 3. Combined master table
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Building master table ==="

MASTER="${OUTDIR}/results/all_hydrogenase_classifications.tsv"
echo -e "data_type\tsource_label\tprotein_id\tcontig_id\thmm_profile\thmm_evalue\thmm_score\tsubject_id\tpct_identity\talign_length\tevalue\tbitscore\thydrogenase_group\tclass_flag" \
    > "${MASTER}"

for TSV in "${OUTDIR}/results/"*_hydrogenase_classification.tsv; do
    [[ "${TSV}" == "${MASTER}" ]] && continue
    tail -n +2 "${TSV}" >> "${MASTER}"
done

echo "Total classified: $(tail -n +2 ${MASTER} | wc -l)"
echo "PASS:             $(grep 'PASS' ${MASTER} | wc -l)"
echo ""
echo "=== PASS hits per group ==="
grep "PASS" "${MASTER}" | awk -F'\t' '{print $13}' | sort | uniq -c | sort -rn
echo ""
echo "=== PASS hits per source ==="
grep "PASS" "${MASTER}" | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn

echo "[$(date)] Done."
echo "Master table: ${MASTER}"
echo ""
echo "Join keys:"
echo "  TranscriptM gene_id → protein_id"
echo "  CoverM contig       → contig_id"
echo "  GTDB-Tk taxonomy    → source_label (MAGs)"
echo "  Taxonkit taxonomy   → contig_id + source_label (assemblies)"

# ════════════════════════════════════════════════════════════
# 3. Unbinned dereplicated contigs
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Unbinned contigs ==="
classify_hydrogenases \
    "30_unbinned_contigs/seqs/unbinned_dereplicated.faa" \
    "unbinned" \
    "unbinned"
