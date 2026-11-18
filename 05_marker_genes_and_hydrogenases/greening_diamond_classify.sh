#!/bin/bash
#SBATCH --job-name=greening_dmnd_all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=20G
#SBATCH --time=144:00:00
#SBATCH --output=greening_dmnd_all_%j.out
#SBATCH --error=greening_dmnd_all_%j.err

# ============================================================
# DIAMOND search of all GreeningLab profiles (2020 + 2024)
# against 159 MAGs and 7 assemblies
#
# Filtering per Ricci et al. 2025:
# Default: qcov>=80%, identity>=50%
# Exceptions (identity): rho=40%, nuoF/mmoX/FeFe/coxL/amoA/nxrA/rbcL=60%
#                        psaA=80%, psbA/isoA/atpA/aro/ygfK=70%, hbsT=75%
# MAG overrides: atpA=60%, psbA=60%, rdhA=45%, cyc2=35%, rho=30%
# MAG length filter: len>=40aa OR qcov>=80% OR scov>=80%
# ============================================================

set -eo pipefail

source [conda_base]/etc/profile.d/conda.sh
conda activate diamond_2.1.10

MAG_FAA="21_coassembly/MAGs/hyd_coxl_classification/prodigal_prefixed/all_159_MAGs_prefixed.faa"
ASSEMBLY_DIR="11_metabolic_assemblies/input_genomes"
DMND_DIR="HMMs/greening_dmnd"
OUTDIR="26_greening_diamond"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}/results"

# ── Identity cutoffs (default = 50%) ─────────────────────────
declare -A ID_DEFAULT
# rho
ID_DEFAULT["Microbial_rhodopsin_RHO"]=40
# 60% group
ID_DEFAULT["NADH-ubiquinone_oxidoreductase_NuoF"]=60
ID_DEFAULT["Soluble_methane_monooxygenase_MmoA"]=60
ID_DEFAULT["FeFe_hydrogenase"]=60
ID_DEFAULT["Fe_hydrogenase"]=60
ID_DEFAULT["NiFe_hydrogenase"]=60
ID_DEFAULT["Carbon_monoxide_dehydrogenase_CoxL"]=60
ID_DEFAULT["Ammonia_monooxygenase_amoA"]=60
ID_DEFAULT["Ammonia_monooxygenase_AmoA_Feb2024"]=60
ID_DEFAULT["Nitrite_oxidoreductase_NxrA"]=60
ID_DEFAULT["Rubisco_RbcL_95_Filtered"]=60
ID_DEFAULT["Particulate_methane_monooxygenase_PmoA"]=60
ID_DEFAULT["Particulate_methane_monooxygenase_PmoA_Feb2024"]=60
ID_DEFAULT["pMMO-like_CUMMO_PxmA_Feb2024"]=60
ID_DEFAULT["Soluble_diiron_monooxygenase_Feb2024"]=60
ID_DEFAULT["FeFe_hydrogenase"]=60
# 70% group
ID_DEFAULT["Photosystem_II_reaction_centre_PsbA_sequences"]=70
ID_DEFAULT["Isoprene_monooxygenase_subunit_A"]=70
ID_DEFAULT["F-type_ATP_synthase_AtpA"]=70
ID_DEFAULT["Arsenite_oxidase_ARO"]=70
ID_DEFAULT["Selenate_reductase_YgfK"]=70
# 75%
ID_DEFAULT["Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences"]=75
# 80%
ID_DEFAULT["Photosystem_I_reaction_centre_PsaA_sequences"]=80

# ── MAG-specific overrides ────────────────────────────────────
declare -A ID_MAG
ID_MAG["F-type_ATP_synthase_AtpA"]=60
ID_MAG["Photosystem_II_reaction_centre_PsbA_sequences"]=60
ID_MAG["Reductive_dehalogenase_RdhA"]=45
ID_MAG["Iron_oxidising_cytochrome_Cyc2"]=35
ID_MAG["Microbial_rhodopsin_RHO"]=30

# ════════════════════════════════════════════════════════════
# Function: search all profiles against one FAA
# ════════════════════════════════════════════════════════════
run_search() {
    local FAA="$1"
    local LABEL="$2"
    local DATA_TYPE="$3"
    local RESULTS="${OUTDIR}/results/${LABEL}_classifications.tsv"

    if [[ ! -f "${FAA}" ]]; then
        echo "  [SKIP] ${LABEL} — FAA not found"
        return
    fi
    if [[ -f "${RESULTS}" ]]; then
        echo "  [SKIP] ${LABEL} — already done ($(grep -c 'PASS' ${RESULTS} 2>/dev/null || echo 0) PASS hits)"
        return
    fi

    echo "  [$(date '+%H:%M:%S')] Processing ${LABEL} (${DATA_TYPE})..."

    echo -e "data_type\tsource_label\tprotein_id\tcontig_id\tprofile\t"\
"subject_id\tpct_identity\tquery_coverage\tsubject_coverage\talign_length\t"\
"evalue\tbitscore\tclass_flag" \
        > "${RESULTS}"

    for DMND in "${DMND_DIR}"/*.dmnd; do
        PROFILE=$(basename "${DMND}" .dmnd)

        # Get identity cutoff
        if [[ "${DATA_TYPE}" == "MAG" && -n "${ID_MAG[${PROFILE}]+x}" ]]; then
            ID_CUT="${ID_MAG[${PROFILE}]}"
        else
            ID_CUT="${ID_DEFAULT[${PROFILE}]:-50}"
        fi

        diamond blastp \
            -q "${FAA}" \
            -d "${DMND}" \
            -o "/tmp/${LABEL}_${PROFILE}_hits.tsv" \
            --max-target-seqs 1 \
            --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length evalue bitscore \
            --query-cover 0 \
            --id 0 \
            --sensitive \
            --threads "${THREADS}" \
            --quiet 2>/dev/null || true

        [[ ! -s "/tmp/${LABEL}_${PROFILE}_hits.tsv" ]] && \
            rm -f "/tmp/${LABEL}_${PROFILE}_hits.tsv" && continue

        while IFS=$'\t' read -r qid sid pct qcov scov alen eval bscore; do
            contig="${qid%_*}"
            flag="LOW"

            if [[ "${DATA_TYPE}" == "MAG" ]]; then
                id_ok=0; cov_ok=0
                (( $(echo "${pct} >= ${ID_CUT}" | bc -l) )) && id_ok=1
                ( (( $(echo "${alen} >= 40" | bc -l) )) || \
                  (( $(echo "${qcov} >= 80" | bc -l) )) || \
                  (( $(echo "${scov} >= 80" | bc -l) )) ) && cov_ok=1
                [[ "${id_ok}" -eq 1 && "${cov_ok}" -eq 1 ]] && flag="PASS"
            else
                ( (( $(echo "${pct} >= ${ID_CUT}" | bc -l) )) && \
                  (( $(echo "${qcov} >= 80" | bc -l) )) ) && flag="PASS"
            fi

            echo -e "${DATA_TYPE}\t${LABEL}\t${qid}\t${contig}\t${PROFILE}\t"\
"${sid}\t${pct}\t${qcov}\t${scov}\t${alen}\t${eval}\t${bscore}\t${flag}"
        done < "/tmp/${LABEL}_${PROFILE}_hits.tsv" >> "${RESULTS}"

        rm -f "/tmp/${LABEL}_${PROFILE}_hits.tsv"
    done

    N_PASS=$(grep -c "PASS" "${RESULTS}" 2>/dev/null || echo 0)
    echo "  ${LABEL}: ${N_PASS} PASS hits across all profiles"
}

# ════════════════════════════════════════════════════════════
# 1. 159 MAGs
# ════════════════════════════════════════════════════════════
echo "[$(date)] === 159 MAGs ==="
run_search "${MAG_FAA}" "MAGs_159" "MAG"

# ════════════════════════════════════════════════════════════
# 2. Assemblies
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Assemblies ==="
for ASM in WIL1_before WIL2_before WIL2_after WILB_D1 WILB_D3 WILB_D10 WILD_D0; do
    run_search "${ASSEMBLY_DIR}/${ASM}.faa" "${ASM}" "assembly"
done

# ════════════════════════════════════════════════════════════
# 3. Combined master table
# ════════════════════════════════════════════════════════════
echo ""

# ════════════════════════════════════════════════════════════
# 3. vOTUs and pOTUs
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === vOTUs ==="
run_search "09_proteins/vOTUs.faa" "vOTUs" "assembly"

echo ""
echo "[$(date)] === pOTUs ==="
run_search "09_proteins/pOTUs.faa" "pOTUs" "assembly"

echo "[$(date)] === Building master table ==="

MASTER="${OUTDIR}/results/all_greening_classifications.tsv"
echo -e "data_type\tsource_label\tprotein_id\tcontig_id\tprofile\t"\
"subject_id\tpct_identity\tquery_coverage\tsubject_coverage\talign_length\t"\
"evalue\tbitscore\tclass_flag" \
    > "${MASTER}"

for TSV in "${OUTDIR}/results/"*_classifications.tsv; do
    [[ "${TSV}" == "${MASTER}" ]] && continue
    tail -n +2 "${TSV}" >> "${MASTER}"
done

echo "Total hits: $(tail -n +2 ${MASTER} | wc -l)"
echo "PASS hits:  $(grep 'PASS' ${MASTER} | wc -l)"
echo ""
echo "=== PASS hits per profile ==="
grep "PASS" "${MASTER}" | awk -F'\t' '{print $5}' | sort | uniq -c | sort -rn
echo ""
echo "=== PASS hits per source ==="
grep "PASS" "${MASTER}" | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn

echo "[$(date)] Done."
echo "Master table: ${MASTER}"
echo ""
echo "Join keys for cross-dataset matching:"
echo "  TranscriptM: protein_id"
echo "  CoverM:      contig_id"
echo "  Taxonkit:    contig_id + source_label (assemblies)"
echo "  GTDB-Tk:     source_label (MAGs)"
