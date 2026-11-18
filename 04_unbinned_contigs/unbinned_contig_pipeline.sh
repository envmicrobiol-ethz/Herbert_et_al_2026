#!/bin/bash
#SBATCH --job-name=unbinned_contigs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=16G
#SBATCH --time=96:00:00
#SBATCH --output=unbinned_contigs_%j.out
#SBATCH --error=unbinned_contigs_%j.err

# ============================================================
# Unbinned contig dereplication and analysis pipeline
#
# Steps:
#   1. Extract unbinned contigs from all assemblies + co-assembly
#   2. Remove contigs that map to MAGs (binned fraction leakage)
#   3. Dereplicate at 99% ANI using MMseqs2
#   4. Create lookup table: representative → all duplicates
#      (so original taxonkit taxonomy can be reused)
#   5. CoverM metaG (read coverage, 7 samples)
#   6. CoverM metaT (TPM, 3 metaT samples)
#   7. GreeningLab DIAMOND gene search
#
# Output: 30_unbinned_contigs/
# ============================================================

set -eo pipefail

source [conda_base]/etc/profile.d/conda.sh
module load stack/.2024-06-silent gcc/12.2.0
module load samtools/1.16.1

ASSEMBLIES_DIR="../Wildstrubel_assemblies"
COASM_FA="21_coassembly/assembly/assembly.fasta"
C2M_ORIG="assembly_stats/contig_to_mag_95.tsv"
C2M_WILBD0="assembly_stats/contig_to_mag_wilbd0.tsv"
COASM_C2M="21_coassembly/MAGs/coasm_contig_to_mag.tsv"
MAG_DIR="21_coassembly/MAGs/drep_95/dereplicated_genomes"
GREENING_DMND_DIR="HMMs/greening_dmnd"
OUTDIR="30_unbinned_contigs"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "${OUTDIR}/seqs" "${OUTDIR}/mmseqs" "${OUTDIR}/coverm" \
         "${OUTDIR}/diamond" "${OUTDIR}/logs"

# ════════════════════════════════════════════════════════════
# 1. Extract unbinned contigs per assembly
# ════════════════════════════════════════════════════════════
COMBINED_UNBINNED="${OUTDIR}/seqs/all_unbinned_raw.fasta"
if [[ -f "${COMBINED_UNBINNED}" ]]; then
    echo "[$(date)] [SKIP] Step 1 — all_unbinned_raw.fasta already exists ($(grep -c '^>' ${COMBINED_UNBINNED}) contigs)"
else
echo "[$(date)] === Step 1: Extracting unbinned contigs ==="

# Build set of all binned contig names (bare contig IDs)
ALL_BINNED="${OUTDIR}/logs/all_binned_contig_ids.txt"
cut -f1 "${C2M_ORIG}" | sort -u > "${ALL_BINNED}"
cut -f1 "${C2M_WILBD0}" | sort -u >> "${ALL_BINNED}"
sort -u "${ALL_BINNED}" -o "${ALL_BINNED}"

# Build set of co-assembly binned contigs
COASM_BINNED="${OUTDIR}/logs/coasm_binned_contig_ids.txt"
cut -f1 "${COASM_C2M}" | sort -u > "${COASM_BINNED}"

COMBINED_UNBINNED="${OUTDIR}/seqs/all_unbinned_raw.fasta"
> "${COMBINED_UNBINNED}"
RENAME_TABLE="${OUTDIR}/logs/contig_rename_table.tsv"
echo -e "original_id\trenamed_id\tsource_assembly" > "${RENAME_TABLE}"

# Per-sample assemblies
for ASM in WIL1_before WIL2_before WIL2_after WILB_D1 WILB_D3 WILB_D10; do
    FA="${ASSEMBLIES_DIR}/${ASM}.fasta"
    [[ ! -f "${FA}" ]] && continue
    echo "  Extracting unbinned: ${ASM}..."

    python3 << PYEOF
import gzip, re
binned = set()
with open("${ALL_BINNED}") as f:
    binned = set(l.strip() for l in f)

out = open("${COMBINED_UNBINNED}", "a")
rename = open("${RENAME_TABLE}", "a")
seq_id = None
seq_lines = []
n_kept = 0

def write_seq():
    global n_kept
    if seq_id and seq_id not in binned:
        renamed = "${ASM}__" + seq_id
        out.write(f">{renamed}\n")
        out.write("".join(seq_lines))
        rename.write(f"{seq_id}\t{renamed}\t${ASM}\n")
        n_kept += 1

with open("${FA}") as f:
    for line in f:
        if line.startswith(">"):
            write_seq()
            seq_id = line[1:].strip().split()[0]
            seq_lines = []
        else:
            seq_lines.append(line)
    write_seq()

out.close()
rename.close()
print(f"  ${ASM}: {n_kept} unbinned contigs extracted")
PYEOF
done

# WILD_D0 (different c2m file)
echo "  Extracting unbinned: WILD_D0..."
python3 << PYEOF
binned = set()
with open("${ALL_BINNED}") as f:
    binned = set(l.strip() for l in f)
with open("${C2M_WILBD0}") as f:
    for line in f:
        binned.add(line.strip().split()[0])

out = open("${COMBINED_UNBINNED}", "a")
rename = open("${RENAME_TABLE}", "a")
seq_id = None; seq_lines = []; n_kept = 0

def write_seq():
    global n_kept
    if seq_id and seq_id not in binned:
        renamed = "WILD_D0__" + seq_id
        out.write(f">{renamed}\n")
        out.write("".join(seq_lines))
        rename.write(f"{seq_id}\t{renamed}\tWILD_D0\n")
        n_kept += 1

with open("${ASSEMBLIES_DIR}/WILD_D0.fasta") as f:
    for line in f:
        if line.startswith(">"):
            write_seq()
            seq_id = line[1:].strip().split()[0]; seq_lines = []
        else:
            seq_lines.append(line)
    write_seq()
out.close(); rename.close()
print(f"  WILD_D0: {n_kept} unbinned contigs extracted")
PYEOF

# Co-assembly unbinned
echo "  Extracting unbinned: co-assembly..."
python3 << PYEOF
binned = set()
with open("${COASM_BINNED}") as f:
    binned = set(l.strip() for l in f)

out = open("${COMBINED_UNBINNED}", "a")
rename = open("${RENAME_TABLE}", "a")
seq_id = None; seq_lines = []; n_kept = 0

def write_seq():
    global n_kept
    if seq_id and seq_id not in binned:
        renamed = "COASM__" + seq_id
        out.write(f">{renamed}\n")
        out.write("".join(seq_lines))
        rename.write(f"{seq_id}\t{renamed}\tco-assembly\n")
        n_kept += 1

with open("${COASM_FA}") as f:
    for line in f:
        if line.startswith(">"):
            write_seq()
            seq_id = line[1:].strip().split()[0]; seq_lines = []
        else:
            seq_lines.append(line)
    write_seq()
out.close(); rename.close()
print(f"  co-assembly: {n_kept} unbinned contigs extracted")
PYEOF

TOTAL_UNBINNED=$(grep -c "^>" "${COMBINED_UNBINNED}")
echo "  Total unbinned contigs before dereplication: ${TOTAL_UNBINNED}"
fi  # end step 1 skip block

# ════════════════════════════════════════════════════════════
# 2. Remove contigs mapping to MAGs (≥80% cov)
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 2: Removing contigs mapping to MAGs ==="
conda activate minimap2_V2.22

# Build prefixed MAG combined FASTA (prefix each contig with MAG name)
# This avoids duplicate contig name errors in samtools/minimap2
MAG_COMBINED="${OUTDIR}/seqs/all_mags_combined_prefixed.fasta"
if [[ ! -f "${MAG_COMBINED}" ]]; then
    echo "  Building prefixed MAG reference FASTA..."
    > "${MAG_COMBINED}"
    for FA in "${MAG_DIR}"/*.fa; do
        MAG=$(basename "${FA}" .fa)
        awk -v mag="${MAG}" '/^>/{print ">"mag"__"substr($0,2); next} {print}' "${FA}"
    done > "${MAG_COMBINED}"
    samtools faidx "${MAG_COMBINED}"
    echo "  MAG contigs: $(grep -c '^>' ${MAG_COMBINED})"
else
    echo "  [SKIP] Prefixed MAG FASTA already exists"
    [[ ! -f "${MAG_COMBINED}.fai" ]] && samtools faidx "${MAG_COMBINED}"
fi

FILTERED="${OUTDIR}/seqs/unbinned_filtered.fasta"
if [[ -f "${FILTERED}" ]]; then
    echo "  [SKIP] Step 2 — unbinned_filtered.fasta already exists ($(grep -c '^>' ${FILTERED}) contigs)"
else
# Map unbinned contigs against MAGs
echo "  Mapping unbinned contigs against MAGs..."
minimap2 -ax asm20 -t "${THREADS}" \
    "${MAG_COMBINED}" \
    "${COMBINED_UNBINNED}" | \
    samtools view -F 4 -q 10 | \
    awk '$3!="*" && $6!="*" {
        # Parse alignment length from CIGAR
        cig=$6; al=0
        while(match(cig,/[0-9]+[MID]/)) {
            s=substr(cig,RSTART,RLENGTH)
            al+=substr(s,1,length(s)-1)
            cig=substr(cig,RSTART+RLENGTH)
        }
        qlen=length($10)
        if(qlen>0 && al/qlen>=0.80) print $1
    }' | sort -u \
    > "${OUTDIR}/logs/contigs_mapping_to_mags.txt"

N_MAG_MAPPED=$(wc -l < "${OUTDIR}/logs/contigs_mapping_to_mags.txt")
echo "  Contigs mapping to MAGs (≥80% cov): ${N_MAG_MAPPED}"

# Filter out MAG-mapping contigs
python3 << PYEOF

exclude = set()
with open("${OUTDIR}/logs/contigs_mapping_to_mags.txt") as f:
    exclude = set(l.strip() for l in f)

out = open("${OUTDIR}/seqs/unbinned_filtered.fasta", "w")
seq_id = None; seq_lines = []; kept = 0; excl = 0

def write_seq():
    global kept, excl
    if seq_id is None: return
    if seq_id in exclude: excl += 1; return
    out.write(f">{seq_id}\n")
    out.write("".join(seq_lines))
    kept += 1

with open("${COMBINED_UNBINNED}") as f:
    for line in f:
        if line.startswith(">"):
            write_seq()
            seq_id = line[1:].strip().split()[0]; seq_lines = []
        else:
            seq_lines.append(line)
    write_seq()
out.close()
print(f"  Kept: {kept}, Excluded (map to MAGs): {excl}")
PYEOF

fi  # end step 2 skip block
conda deactivate

# ════════════════════════════════════════════════════════════
# 3. Dereplicate at 99% ANI with MMseqs2
# ════════════════════════════════════════════════════════════
DEDUP="${OUTDIR}/seqs/unbinned_dereplicated.fasta"
MMSEQS_OUT="${OUTDIR}/mmseqs"

if [[ -f "${DEDUP}" ]]; then
    echo "[$(date)] [SKIP] Step 3 — unbinned_dereplicated.fasta already exists ($(grep -c '^>' ${DEDUP}) contigs)"
else
echo ""
echo "[$(date)] === Step 3: Dereplicating at 99% ANI ==="
conda activate mmseqs2

echo "  Input contigs: $(grep -c '^>' ${FILTERED})"

# Create MMseqs2 database
mmseqs createdb "${FILTERED}" "${MMSEQS_OUT}/contigs_db" \
    --threads "${THREADS}"

# Cluster at 99% ANI, 80% coverage
mmseqs linclust \
    "${MMSEQS_OUT}/contigs_db" \
    "${MMSEQS_OUT}/cluster_result" \
    "${MMSEQS_OUT}/tmp" \
    --min-seq-id 0.99 \
    --cov-mode 1 \
    -c 0.80 \
    --threads "${THREADS}"

# Get cluster TSV (representative → member)
mmseqs createtsv \
    "${MMSEQS_OUT}/contigs_db" \
    "${MMSEQS_OUT}/contigs_db" \
    "${MMSEQS_OUT}/cluster_result" \
    "${MMSEQS_OUT}/cluster_table.tsv" \
    --threads "${THREADS}"

# Extract representative sequences
mmseqs result2repseq \
    "${MMSEQS_OUT}/contigs_db" \
    "${MMSEQS_OUT}/cluster_result" \
    "${MMSEQS_OUT}/rep_db" \
    --threads "${THREADS}"

mmseqs convert2fasta \
    "${MMSEQS_OUT}/rep_db" \
    "${OUTDIR}/seqs/unbinned_dereplicated.fasta"

N_REPS=$(grep -c "^>" "${OUTDIR}/seqs/unbinned_dereplicated.fasta")
N_INPUT=$(grep -c "^>" "${FILTERED}")
echo "  Representatives after 99% dereplication: ${N_REPS} / ${N_INPUT}"

# Build duplicate lookup table
# Format: representative_id TAB member_id TAB source_assembly
echo -e "representative_id\tmember_id\tsource_assembly" \
    > "${OUTDIR}/logs/dedup_cluster_table.tsv"
awk 'BEGIN{OFS="\t"} {
    # Extract source from prefix (before __)
    split($1, a, "__"); rep_src=a[1]
    split($2, b, "__"); mem_src=b[1]
    print $1, $2, mem_src
}' "${MMSEQS_OUT}/cluster_table.tsv" >> "${OUTDIR}/logs/dedup_cluster_table.tsv"

echo "  Cluster table: ${OUTDIR}/logs/dedup_cluster_table.tsv"
fi  # end step 3 skip block

conda deactivate

# ════════════════════════════════════════════════════════════
# 4. CoverM metaG — all 7 samples
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 4: CoverM metaG ==="
if [[ -f "${OUTDIR}/coverm/unbinned_metaG_coverage.tsv" ]]; then
    echo "[SKIP] Step 4 — metaG coverage already done"
else
conda activate coverm_0.6.1

DEDUP="${OUTDIR}/seqs/unbinned_dereplicated.fasta"

coverm contig \
    --reference "${DEDUP}" \
    --single \
        "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz" \
        "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode09_filtered.fastq.gz" \
        "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode10_filtered.fastq.gz" \
        "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz" \
        "[reads_dir_run2]/fastq_trimmed/trimmed_filtered/barcode17_filtered.fastq.gz" \
        "[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode11_filtered.fastq.gz" \
        "[reads_dir_run1]/fastq_trimmed/trimmed_filtered/barcode12_filtered.fastq.gz" \
    --methods mean covered_bases rpkm \
    --min-covered-fraction 0 \
    --min-read-percent-identity 85 \
    --mapper minimap2-ont \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/coverm/unbinned_metaG_coverage.tsv"

echo "  metaG done: $(tail -n +2 ${OUTDIR}/coverm/unbinned_metaG_coverage.tsv | wc -l) contigs"
fi  # end step 4 skip

# ════════════════════════════════════════════════════════════
# 5. CoverM metaT — 3 metaT samples
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 5: CoverM metaT ==="
if [[ -f "${OUTDIR}/coverm/unbinned_metaT_coverage.tsv" ]]; then
    echo "[SKIP] Step 5 — metaT coverage already done"
else

READS_DIR="[metaT_reads_dir]"

coverm contig \
    --reference "${DEDUP}" \
    --coupled \
        "${READS_DIR}/WAf15_8_R1.trimmed.fastq.gz" \
        "${READS_DIR}/WAf15_8_R2.trimmed.fastq.gz" \
    --coupled \
        "${READS_DIR}/WBe15_8_R1.trimmed.fastq.gz" \
        "${READS_DIR}/WBe15_8_R2.trimmed.fastq.gz" \
    --coupled \
        "${READS_DIR}/WBe218_8_R1.trimmed.fastq.gz" \
        "${READS_DIR}/WBe218_8_R2.trimmed.fastq.gz" \
    --methods mean covered_bases rpkm \
    --min-covered-fraction 0 \
    --min-read-percent-identity 90 \
    --mapper minimap2-sr \
    --threads "${THREADS}" \
    --output-file "${OUTDIR}/coverm/unbinned_metaT_coverage.tsv"
fi  # end step 5 skip

echo "  metaT done: $(tail -n +2 ${OUTDIR}/coverm/unbinned_metaT_coverage.tsv | wc -l) contigs"
conda deactivate

# ════════════════════════════════════════════════════════════
# 6. GreeningLab DIAMOND gene search
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 6: GreeningLab DIAMOND gene search ==="
module load stack/.2024-05-silent gcc/13.2.0 prodigal/2.6.3
conda activate diamond_2.1.10

# Predict proteins
FAA="${OUTDIR}/seqs/unbinned_dereplicated.faa"
if [[ ! -f "${FAA}" ]]; then
    echo "  Running Prodigal..."
    prodigal \
        -i "${DEDUP}" \
        -a "${FAA}" \
        -p meta -f gff \
        -o "${OUTDIR}/seqs/unbinned_dereplicated.gff" -q
    echo "  Proteins: $(grep -c '^>' ${FAA})"
fi

# DIAMOND against all GreeningLab profiles
mkdir -p "${OUTDIR}/diamond/results"
echo -e "contig_id\tprotein_id\tprofile\tsubject_id\tpct_identity\tquery_coverage\tsubject_coverage\talign_length\tevalue\tbitscore\tclass_flag" \
    > "${OUTDIR}/diamond/results/unbinned_greening_classifications.tsv"

# Per-profile identity cutoffs (same as MAG/assembly analysis)
declare -A ID_CUT
ID_CUT["Microbial_rhodopsin_RHO"]=40
ID_CUT["NADH-ubiquinone_oxidoreductase_NuoF"]=60
ID_CUT["Soluble_methane_monooxygenase_MmoA"]=60
ID_CUT["Carbon_monoxide_dehydrogenase_CoxL"]=60
ID_CUT["Ammonia_monooxygenase_amoA"]=60
ID_CUT["Ammonia_monooxygenase_AmoA_Feb2024"]=60
ID_CUT["Nitrite_oxidoreductase_NxrA"]=60
ID_CUT["Rubisco_RbcL_95_Filtered"]=60
ID_CUT["Particulate_methane_monooxygenase_PmoA"]=60
ID_CUT["Particulate_methane_monooxygenase_PmoA_Feb2024"]=60
ID_CUT["pMMO-like_CUMMO_PxmA_Feb2024"]=60
ID_CUT["Photosystem_I_reaction_centre_PsaA_sequences"]=80
ID_CUT["Photosystem_II_reaction_centre_PsbA_sequences"]=70
ID_CUT["Isoprene_monooxygenase_subunit_A"]=70
ID_CUT["ATP-citrate_lyase_AclB_sequences"]=70
ID_CUT["Arsenite_oxidase_ARO"]=70
ID_CUT["Reductive_dehalogenase_RdhA"]=70
ID_CUT["Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences"]=75
ID_CUT["F-type_ATP_synthase_AtpA"]=70
ID_CUT["Selenate_reductase_YgfK"]=70

for DMND in "${GREENING_DMND_DIR}"/*.dmnd; do
    PROFILE=$(basename "${DMND}" .dmnd)
    CUT="${ID_CUT[${PROFILE}]:-50}"

    diamond blastp \
        -q "${FAA}" \
        -d "${DMND}" \
        -o "/tmp/unbinned_${PROFILE}_hits.tsv" \
        --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length evalue bitscore \
        --query-cover 0 --id 0 \
        --threads "${THREADS}" \
        --quiet 2>/dev/null || true

    [[ ! -s "/tmp/unbinned_${PROFILE}_hits.tsv" ]] && continue

    # Assembly filter: identity >= cutoff AND (len>=40 OR qcov>=80 OR scov>=80)
    awk -v profile="${PROFILE}" -v cut="${CUT}" 'BEGIN{OFS="\t"} {
        contig=$1; sub(/_[0-9]+$/, "", contig)
        flag=(($3>=cut) && ($6>=40 || $4>=80 || $5>=80)) ? "PASS" : "LOW"
        print contig,$1,profile,$2,$3,$4,$5,$6,$7,$8,flag
    }' "/tmp/unbinned_${PROFILE}_hits.tsv" \
        >> "${OUTDIR}/diamond/results/unbinned_greening_classifications.tsv"

    rm -f "/tmp/unbinned_${PROFILE}_hits.tsv"
done

N_PASS=$(grep -c "PASS" "${OUTDIR}/diamond/results/unbinned_greening_classifications.tsv" || echo 0)
echo "  PASS hits: ${N_PASS}"
echo ""
echo "  Top profiles:"
grep "PASS" "${OUTDIR}/diamond/results/unbinned_greening_classifications.tsv" | \
    awk -F'\t' '{print $3}' | sort | uniq -c | sort -rn | head -15

# ════════════════════════════════════════════════════════════
# 7. Taxonomy lookup table
# ════════════════════════════════════════════════════════════
echo ""
echo "[$(date)] === Step 7: Building taxonomy lookup ==="

# Create lookup pointing to original taxonkit output
# Key: renamed_id → original_contig_id → taxonkit file
conda activate transcriptm_v0.3.1
python3 << PYEOF
import pandas as pd, os

rename = pd.read_csv("${OUTDIR}/logs/contig_rename_table.tsv", sep="\t")
cluster = pd.read_csv("${OUTDIR}/logs/dedup_cluster_table.tsv", sep="\t")

# Map representative → original contig → taxonkit file
taxonkit_files = {
    "WIL1_before": "24_taxonkit_assemblies/WIL1_before_read_GTDB_r220_taxonomy.txt",
    "WIL2_before": "24_taxonkit_assemblies/WIL2_before_read_GTDB_r220_taxonomy.txt",
    "WIL2_after":  "24_taxonkit_assemblies/WIL2_after_read_GTDB_r220_taxonomy.txt",
    "WILB_D1":     "24_taxonkit_assemblies/WILB_D1_read_GTDB_r220_taxonomy.txt",
    "WILB_D3":     "24_taxonkit_assemblies/WILB_D3_read_GTDB_r220_taxonomy.txt",
    "WILB_D10":    "24_taxonkit_assemblies/WILB_D10_read_GTDB_r220_taxonomy.txt",
    "WILD_D0":     "24_taxonkit_assemblies/WILD_D0_read_GTDB_r220_taxonomy.txt",
    "co-assembly": "24_taxonkit_assemblies/COASM_read_GTDB_r220_taxonomy.txt",
}

# Build lookup: representative_id → original_contig_id, taxonkit_file
rename["source"] = rename["source_assembly"]
rename_dict = rename.set_index("renamed_id").to_dict("index")

rows = []
for _, row in cluster.iterrows():
    rep = row["representative_id"]
    mem = row["member_id"]
    src = row["source_assembly"]
    orig = rename_dict.get(mem, {}).get("original_id", mem)
    taxf = taxonkit_files.get(src, "N/A")
    rows.append({"representative_id": rep, "member_id": mem,
                 "original_contig_id": orig, "source_assembly": src,
                 "taxonkit_file": taxf})

df = pd.DataFrame(rows)
df.to_csv("${OUTDIR}/logs/taxonomy_lookup.tsv", sep="\t", index=False)
print(f"Taxonomy lookup: {len(df)} entries")
print(f"Representatives with taxonomy available: {(df['taxonkit_file']!='N/A').sum()}")
PYEOF

echo ""
echo "[$(date)] === Pipeline complete ==="
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Key files:"
echo "  Dereplicated contigs:    ${OUTDIR}/seqs/unbinned_dereplicated.fasta"
echo "  Duplicate lookup:        ${OUTDIR}/logs/dedup_cluster_table.tsv"
echo "  Taxonomy lookup:         ${OUTDIR}/logs/taxonomy_lookup.tsv"
echo "  metaG coverage:          ${OUTDIR}/coverm/unbinned_metaG_coverage.tsv"
echo "  metaT coverage:          ${OUTDIR}/coverm/unbinned_metaT_coverage.tsv"
echo "  Gene classifications:    ${OUTDIR}/diamond/results/unbinned_greening_classifications.tsv"
echo ""
echo "Join keys for cross-dataset matching:"
echo "  Gene results → contig_id (matches dereplicated contig names)"
echo "  Coverage     → contig name (same)"
echo "  Taxonomy     → use taxonomy_lookup.tsv to map rep → original_contig → taxonkit output"
