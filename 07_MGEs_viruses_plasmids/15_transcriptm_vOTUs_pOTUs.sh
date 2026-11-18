#!/bin/bash
#SBATCH --job-name=transcriptm_vp
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=144:00:00
#SBATCH --output=transcriptm_vp_%j.out
#SBATCH --error=transcriptm_vp_%j.err

source [conda_base]/etc/profile.d/conda.sh
conda activate transcriptm_v0.3.1

READS_DIR="[metaT_reads_dir]"
VOTUS_DIR="07_viruses/vOTUs"
POTUS_DIR="08_plasmids/pOTUs/dereplicated_genomes"
OUTDIR="15_transcriptm"

# ── Prepare vOTU split directory ──────────────────────────────
VOTU_SPLIT="${OUTDIR}/vOTUs_split"
if [[ ! -d "${VOTU_SPLIT}" ]] || [[ $(ls "${VOTU_SPLIT}"/*.fasta 2>/dev/null | wc -l) -eq 0 ]]; then
    mkdir -p "${VOTU_SPLIT}"
    awk -v outdir="${VOTU_SPLIT}" \
        '/^>/{if(file) close(file);
              name=substr($0,2); gsub(/[^A-Za-z0-9._-]/,"_",name);
              file=outdir"/"name".fasta"}
         {print > file}' "${VOTUS_DIR}/vOTUs.fasta"
    echo "vOTUs split: $(ls ${VOTU_SPLIT}/*.fasta | wc -l) files"
fi

# ── Prepare pOTU directory with .fasta extension ──────────────
POTU_FASTA="${OUTDIR}/pOTUs_fasta"
if [[ ! -d "${POTU_FASTA}" ]] || [[ $(ls "${POTU_FASTA}"/*.fasta 2>/dev/null | wc -l) -eq 0 ]]; then
    mkdir -p "${POTU_FASTA}"
    for FA in "${POTUS_DIR}"/*.fa; do
        [[ -f "${FA}" ]] || continue
        cp "${FA}" "${POTU_FASTA}/$(basename ${FA} .fa).fasta"
    done
    echo "pOTUs ready: $(ls ${POTU_FASTA}/*.fasta | wc -l) files"
fi

# ── vOTUs ─────────────────────────────────────────────────────
for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do
    R1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    OUT="${OUTDIR}/${SAMPLE}_vOTUs"
    echo "Running TranscriptM for sample: ${SAMPLE} vs vOTUs"
    transcriptm count \
        -1 "$R1" \
        -2 "$R2" \
        -g "$VOTU_SPLIT" \
        -x fasta \
        -o "$OUT" \
        --skip-qc
done

# ── pOTUs ─────────────────────────────────────────────────────
for SAMPLE in WAf15_8 WBe15_8 WBe218_8; do
    R1="${READS_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    OUT="${OUTDIR}/${SAMPLE}_pOTUs"
    echo "Running TranscriptM for sample: ${SAMPLE} vs pOTUs"
    transcriptm count \
        -1 "$R1" \
        -2 "$R2" \
        -g "$POTU_FASTA" \
        -x fasta \
        -o "$OUT" \
        --skip-qc
done

echo "=== ALL TranscriptM vOTU/pOTU runs completed ==="
