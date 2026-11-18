#!/bin/bash
# ============================================================
# Rename assembly contigs to <=20 characters for Prokka
# Keeps a conversion table to map back to original names
# ============================================================
# Run from MAGs directory before submitting TranscriptM
# ============================================================

set -euo pipefail

ASSEMBLIES_DIR="15_transcriptm_assemblies/renamed_assemblies"
OUTDIR="15_transcriptm_assemblies/prokka_ready_assemblies"
CONV_TABLE="${OUTDIR}/assembly_contig_conversion_table.tsv"

mkdir -p "${OUTDIR}"
echo -e "assembly\toriginal_contig\trenamed_contig" > "${CONV_TABLE}"

for FASTA in "${ASSEMBLIES_DIR}"/*.fasta; do
    [[ -f "${FASTA}" ]] || continue
    ASM=$(basename "${FASTA}" .fasta)
    OUT="${OUTDIR}/${ASM}.fasta"

    if [[ -f "${OUT}" ]]; then
        echo "[SKIP] ${ASM} — already renamed"
        continue
    fi

    echo "[$(date '+%H:%M:%S')] Renaming: ${ASM}..."
    > "${OUT}"
    seq_count=1

    while read -r line; do
        if [[ "${line}" == ">"* ]]; then
            ORIG="${line#>}"
            # Short unique name: first 8 chars of asm + seq number
            ASM_SHORT=$(echo "${ASM}" | sed 's/[^A-Za-z0-9]//g' | cut -c1-8)
            NEW_NAME="${ASM_SHORT}_${seq_count}"
            # Ensure <=20 chars
            NEW_NAME=$(echo "${NEW_NAME}" | cut -c1-20)
            echo ">${NEW_NAME}" >> "${OUT}"
            echo -e "${ASM}\t${ORIG}\t${NEW_NAME}" >> "${CONV_TABLE}"
            ((seq_count++))
        else
            echo "${line}" >> "${OUT}"
        fi
    done < "${FASTA}"

    N=$(grep -c "^>" "${OUT}")
    echo "  ${ASM}: ${N} contigs renamed"
done

echo ""
echo "Done. Prokka-ready assemblies: ${OUTDIR}/"
echo "Conversion table: ${CONV_TABLE}"

# Check max name length
echo ""
echo "Max contig name length:"
grep "^>" "${OUTDIR}"/*.fasta | awk -F'>' '{print length($2)}' | sort -n | tail -1
