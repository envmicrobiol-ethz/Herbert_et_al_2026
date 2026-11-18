#!/usr/bin/env python3
"""
Extract CoxL sequences for phylogenetic tree.
- MAGs: Prokka FAA, filtered to 159 final MAGs
- Unbinned: Prodigal FAA from co-assembly, filtered to 30_unbinned
- Plasmid: pOTU WILB_D3__contig_507_14 (manually validated)
- References: KEGG K03520 only (no GreeningLab)
"""

import pandas as pd, re, os, glob

OUTDIR = "34_phylogenetics/coxL"
os.makedirs(OUTDIR, exist_ok=True)

# ── 159 MAG set ───────────────────────────────────────────
mag_159 = set(os.path.basename(f).replace(".fa","")
              for f in glob.glob("21_coassembly/MAGs/drep_95/dereplicated_genomes/*.fa"))
print(f"159 MAGs: {len(mag_159)}")

# ── GTDB taxonomy ──────────────────────────────────────────
gtdb = pd.concat([
    pd.read_csv("23_phylogenetics/gtdbtk_full/gtdbtk.bac120.summary.tsv",
                sep="\t", usecols=["user_genome","classification"]),
    pd.read_csv("23_phylogenetics/gtdbtk_full/gtdbtk.ar53.summary.tsv",
                sep="\t", usecols=["user_genome","classification"])
], ignore_index=True)
gtdb_dict = gtdb.set_index("user_genome")["classification"].to_dict()

# ── MAG name lookups ──────────────────────────────────────
mag_lookup   = pd.read_csv("21_coassembly/MAGs/transcriptm/coasm_mag_id_lookup.tsv", sep="\t")
orig_lookup  = pd.read_csv("prokka_ready_mags/original_mag_id_lookup.tsv", sep="\t")
wilbd0_lookup= pd.read_csv("15_transcriptm_wilbd0/wilbd0_mag_id_lookup.tsv", sep="\t")

# ── Build gene_id → MAG name + contig via feature counts ──
gene_to_mag = {}
gene_to_contig_fc = {}
for pattern, lkdf in [
    ("21_coassembly/MAGs/transcriptm/WAf15_8_MAGs/feature_count/feature_count_table", mag_lookup),
    ("15_transcriptm/WAf15_8_MAGs/feature_count/feature_count_table", orig_lookup),
    ("15_transcriptm_wilbd0/WAf15_8_MAGs/feature_count/feature_count_table", wilbd0_lookup),
]:
    lk = lkdf.set_index("mag_id")["mag_name"].to_dict()
    try:
        fc = pd.read_csv(pattern, sep="\t")
        ncols = len(fc.columns)
        fc.columns = ["gene_id","contig","start","end","strand","len"] + \
                     [f"a{i}" for i in range(ncols-7)] + ["count"]
        for _, r in fc.iterrows():
            contig = str(r["contig"])
            if "_seq_" in contig:
                raw_id = contig.split("_seq_")[0]
                m = re.match(r"coasm(\d+)$", raw_id)
                if m: mag_id = f"coasm_mag_{m.group(1)}"
                else:
                    m2 = re.match(r"mag(\d+)$", raw_id)
                    if m2: mag_id = "mag_" + m2.group(1)
                    else:
                        m3 = re.match(r"wilbd0_mag(\d+)$", raw_id)
                        mag_id = "wilbd0_mag_" + m3.group(1) if m3 else raw_id
                mag_name = lk.get(mag_id,"")
                if mag_name:
                    gene_to_mag[r["gene_id"]] = mag_name
                    gene_to_contig_fc[r["gene_id"]] = contig
    except Exception as e:
        print(f"  Warning: {e}")

print(f"Gene→MAG lookup: {len(gene_to_mag)}")

# ── Prokka coord → gene_id lookup ─────────────────────────
prodigal_coords = {}
faa = "21_coassembly/MAGs/hyd_coxl_classification/prodigal_prefixed/all_159_MAGs_prefixed.faa"
with open(faa) as f:
    for line in f:
        if not line.startswith(">"): continue
        parts = line[1:].strip().split(" # ")
        if len(parts) < 3: continue
        prodigal_coords[parts[0].strip()] = (int(parts[1].strip()), int(parts[2].strip()))

# Conversion tables
def load_conv(fpath):
    df = pd.read_csv(fpath, sep="\t")
    df["bare"] = df["Old_Header"].str.replace(">","")
    df["renamed"] = df["New_Header"].str.replace(">","")
    return df.set_index("bare")["renamed"].to_dict()

conv_coasm  = load_conv("21_coassembly/MAGs/transcriptm/coasm_mag_conversion_table.tsv")
conv_orig   = load_conv("prokka_ready_mags/mag_conversion_table.tsv")
conv_wilbd0 = load_conv("15_transcriptm_wilbd0/wilbd0_mag_conversion_table.tsv")

def build_coord_lookup(fpath):
    lookup = {}
    try:
        fc = pd.read_csv(fpath, sep="\t")
        ncols = len(fc.columns)
        fc.columns = ["gene_id","contig","start","end","strand","len"] + \
                     [f"a{i}" for i in range(ncols-7)] + ["count"]
        for _, r in fc.iterrows():
            lookup[(r["contig"], int(r["start"]), int(r["end"]))] = r["gene_id"]
    except: pass
    return lookup

coord_coasm  = build_coord_lookup("21_coassembly/MAGs/transcriptm/WAf15_8_MAGs/feature_count/feature_count_table")
coord_orig   = build_coord_lookup("15_transcriptm/WAf15_8_MAGs/feature_count/feature_count_table")
coord_wilbd0 = build_coord_lookup("15_transcriptm_wilbd0/WAf15_8_MAGs/feature_count/feature_count_table")

def get_prokka_id(prodigal_pid):
    coords = prodigal_coords.get(prodigal_pid)
    if not coords: return None
    start, end = coords
    bare = re.search(r"contig_(\d+)_\d+$", prodigal_pid)
    if not bare: return None
    bare_contig = f"contig_{bare.group(1)}"
    if prodigal_pid.startswith("COASM_"):
        renamed = conv_coasm.get(bare_contig)
        coord_lk = coord_coasm
    elif prodigal_pid.startswith("WILBD0_"):
        renamed = conv_wilbd0.get(bare_contig)
        coord_lk = coord_wilbd0
    else:
        renamed = conv_orig.get(bare_contig)
        coord_lk = coord_orig
    return coord_lk.get((renamed, start, end)) if renamed else None

# ── Load Prokka CoxL PASS hits ────────────────────────────
prokka = pd.read_csv(
    "26_greening_diamond/prokka/results/all_prokka_greening_classifications.tsv", sep="\t")
coxl = prokka[(prokka["class_flag"]=="PASS") &
              (prokka["profile"]=="Carbon_monoxide_dehydrogenase_CoxL")].copy()
coxl = coxl.sort_values("bitscore",ascending=False).drop_duplicates("gene_id",keep="first")

# Filter to 159 MAGs only
rows = []
for _, r in coxl[coxl["source_label"].isin(["coasm_mags","original_mags","wilbd0_mags"])].iterrows():
    mag_name = gene_to_mag.get(r["gene_id"],"")
    if not mag_name or mag_name not in mag_159: continue
    prokka_id = r["gene_id"]
    # Get Prodigal protein ID for sequence extraction
    # Find via coord reverse lookup
    rows.append({
        "prokka_id": prokka_id,
        "mag_name": mag_name,
        "source": "MAG",
        "taxonomy": gtdb_dict.get(mag_name,"unknown"),
        "pct_identity": r["pct_identity"],
        "subject_id": r.get("subject_id",""),
    })

print(f"MAG CoxL (159 only): {len(rows)}")

# ── Co-assembly unbinned CoxL ─────────────────────────────
lookup = pd.read_csv("30_unbinned_contigs/logs/taxonomy_lookup_complete.tsv", sep="\t")
unbinned_reps = set(lookup["representative_id"])
coasm_tax = pd.read_csv(
    "24_taxonkit_assemblies/COASM_unbinned_proteins_read_GTDB_r220_taxonomy.txt",
    sep="\t", header=None, names=["contig","taxid","taxonomy"])
coasm_tax_dict = coasm_tax.set_index("contig")["taxonomy"].to_dict()

coasm_ub = pd.read_csv(
    "21_coassembly/prodigal/diamond/coasm_unbinned_greening_classifications.tsv", sep="\t")
coxl_ub = coasm_ub[(coasm_ub["class_flag"]=="PASS") &
                   (coasm_ub["profile"]=="Carbon_monoxide_dehydrogenase_CoxL")].copy()
coxl_ub = coxl_ub.sort_values("bitscore",ascending=False).drop_duplicates("protein_id",keep="first")

for _, r in coxl_ub.iterrows():
    cid = r["contig_id"]
    if cid not in unbinned_reps: continue
    bare = cid.replace("COASM__","")
    rows.append({
        "prokka_id": r["protein_id"],
        "mag_name": cid,
        "source": "unbinned",
        "taxonomy": coasm_tax_dict.get(bare,"unknown"),
        "pct_identity": r["pct_identity"],
        "subject_id": r.get("subject_id",""),
    })

n_ub = sum(1 for r in rows if r["source"]=="unbinned")
print(f"Unbinned CoxL (30_unbinned): {n_ub}")

# Add plasmid CoxL
rows.append({
    "prokka_id": "WILB_D3__contig_507_14",
    "mag_name": "WILB_D3__contig_507",
    "source": "plasmid",
    "taxonomy": "unknown",
    "pct_identity": 37.0,
    "subject_id": "MANUAL_HMM_VALIDATED",
})
print(f"Total: {len(rows)}")

# ── Extract sequences ─────────────────────────────────────
all_ids = set(r["prokka_id"] for r in rows)

FAA_SOURCES = [
    "26_greening_diamond/prokka/faa/coasm_mags_prokka.faa",
    "26_greening_diamond/prokka/faa/original_mags_prokka.faa",
    "26_greening_diamond/prokka/faa/wilbd0_mags_prokka.faa",
    "21_coassembly/prodigal/COASM_full_renamed.faa",
    "09_proteins/pOTUs.faa",
]

seqs = {}
for faa_path in FAA_SOURCES:
    if not os.path.exists(faa_path): continue
    seq_id = None; seq_lines = []
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                if seq_id and seq_id in all_ids:
                    seqs[seq_id] = "".join(seq_lines)
                seq_id = line[1:].strip().split()[0]; seq_lines = []
            else:
                seq_lines.append(line.strip())
        if seq_id and seq_id in all_ids:
            seqs[seq_id] = "".join(seq_lines)

print(f"Sequences found: {len(seqs)} / {len(all_ids)}")

# ── Write FASTA with informative labels ───────────────────
ann_dict = {r["prokka_id"]: r for r in rows}
out_faa = f"{OUTDIR}/wildstrubel_coxL.faa"
n = 0
with open(out_faa, "w") as f:
    for pid, seq in seqs.items():
        a = ann_dict.get(pid, {})
        mag = a.get("mag_name","")
        src = a.get("source","")
        tax = str(a.get("taxonomy","unknown"))
        # Shorten taxonomy to class level
        tax_short = ";".join(tax.split(";")[:4]).replace(" ","_")
        # Clean label
        label = f"{pid}|{src}|{mag}|{tax_short}".replace(" ","_")
        f.write(f">{label}\n{seq}\n")
        n += 1

print(f"Wildstrubel CoxL FASTA: {n} sequences → {out_faa}")

# Save annotation table
pd.DataFrame(rows).to_csv(f"{OUTDIR}/coxL_annotations.tsv", sep="\t", index=False)
print(f"Annotations: {OUTDIR}/coxL_annotations.tsv")
