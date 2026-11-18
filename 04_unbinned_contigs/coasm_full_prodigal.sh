#!/bin/bash
#SBATCH --job-name=coasm_prodigal_full
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=50G
#SBATCH --time=24:00:00
#SBATCH --output=coasm_prodigal_full_%j.out
#SBATCH --error=coasm_prodigal_full_%j.err

# ============================================================
# Run Prodigal on full co-assembly to generate a single
# consistent GFF + FAA that can be used for both:
#   1. GreeningLab/HydDB DIAMOND classification
#   2. TranscriptM --ref/--gff input
#
# This ensures gene IDs match exactly between classification
# and expression analysis.
#
# Output: 21_coassembly/prodigal/
#   COASM_full.gff  — gene coordinates (use with TranscriptM --gff)
#   COASM_full.faa  — protein sequences (use for DIAMOND)
#   COASM_full.fna  — nucleotide sequences
#   COASM_full_renamed.fna — with COASM__ prefix for TranscriptM --ref
# ============================================================

set -eo pipefail

module load stack/.2024-05-silent gcc/13.2.0 prodigal/2.6.3

COASM_FA="21_coassembly/assembly/assembly.fasta"
OUTDIR="21_coassembly/prodigal"
mkdir -p "${OUTDIR}"

if [[ -f "${OUTDIR}/COASM_full.gff" ]]; then
    echo "[$(date)] [SKIP] Prodigal already done"
    echo "  GFF: $(grep -c $'\tCDS\t' ${OUTDIR}/COASM_full.gff) CDS"
    echo "  FAA: $(grep -c '^>' ${OUTDIR}/COASM_full.faa) proteins"
    exit 0
fi

echo "[$(date)] Running Prodigal on full co-assembly (211,503 contigs)..."
echo "  Input: ${COASM_FA}"

prodigal \
    -i "${COASM_FA}" \
    -f gff \
    -o "${OUTDIR}/COASM_full.gff" \
    -a "${OUTDIR}/COASM_full.faa" \
    -d "${OUTDIR}/COASM_full.fna" \
    -p meta \
    -q

echo "[$(date)] Done."
echo "  CDS: $(grep -c $'\tCDS\t' ${OUTDIR}/COASM_full.gff)"
echo "  FAA: $(grep -c '^>' ${OUTDIR}/COASM_full.faa) proteins"

# Build renamed ref FASTA for TranscriptM (COASM__ prefix)
echo "[$(date)] Building renamed reference FASTA for TranscriptM..."
awk '/^>/{print ">COASM__"substr($0,2); next}{print}' \
    "${COASM_FA}" > "${OUTDIR}/COASM_full_renamed.fna"

# Build renamed GFF to match
sed 's/^contig_/COASM__contig_/' \
    "${OUTDIR}/COASM_full.gff" > "${OUTDIR}/COASM_full_renamed.gff"

# Build renamed FAA to match
sed 's/^>contig_/>COASM__contig_/' \
    "${OUTDIR}/COASM_full.faa" > "${OUTDIR}/COASM_full_renamed.faa"

echo "[$(date)] Renamed files ready for TranscriptM."
echo ""
echo "Next steps:"
echo "  1. Run DIAMOND GreeningLab + HydDB on ${OUTDIR}/COASM_full_renamed.faa"
echo "  2. Filter PASS hits to unbinned contigs"
echo "  3. Extract marker contigs from co-assembly"
echo "  4. Run TranscriptM with:"
echo "       --ref ${OUTDIR}/COASM_full_renamed.fna (or marker subset)"
echo "       --gff ${OUTDIR}/COASM_full_renamed.gff (or marker subset)"
echo "  5. Join: TranscriptM gene_id → DIAMOND protein_id via seqnum in GFF"
