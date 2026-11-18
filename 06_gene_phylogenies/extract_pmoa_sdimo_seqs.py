#!/usr/bin/env python3
"""
Extract protein sequences for pmoA and soluble monooxygenase PASS hits
from MAGs and unbinned contigs, with taxonomy annotations.
"""

import pandas as pd
import os

OUTDIR = "32_pmoA_sdimo"
os.makedirs(OUTDIR, exist_ok=True)

profiles_of_interest = [
    "Particulate_methane_monooxygenase_PmoA",
    "Particulate_methane_monooxygenase_PmoA_Feb2024",
    "Soluble_methane_monooxygenase_MmoA",
    "Soluble_diiron_monooxygenase_Feb2024",
    "pMMO-like_CUMMO_PxmA_Feb2024"
]

# ── Load classification results ────────────────────────────────
greening_mag = pd.read_csv("26_greening_diamond/results/all_greening_classifications.tsv", sep="\t")
mag_pass = greening_mag[(greening_mag["class_flag"]=="PASS") &
               (greening_mag["data_type"]=="MAG") &
               (greening_mag["profile"].isin(profiles_of_interest))].copy()
mag_pass = mag_pass.sort_values("bitscore", ascending=False).drop_duplicates(subset="protein_id", keep="first")

greening_ub = pd.read_csv("30_unbinned_contigs/diamond/results/unbinned_greening_classifications.tsv", sep="\t")
ub_pass = greening_ub[(greening_ub["class_flag"]=="PASS") &
             (greening_ub["profile"].isin(profiles_of_interest))].copy()
ub_pass = ub_pass.sort_values("bitscore", ascending=False).drop_duplicates(subset="protein_id", keep="first")

print(f"MAG PASS hits: {len(mag_pass)}")
print(f"Unbinned PASS hits: {len(ub_pass)}")

# ── Load GTDB taxonomy for MAGs ────────────────────────────────
gtdb_bac = pd.read_csv("23_phylogenetics/gtdbtk_full/gtdbtk.bac120.summary.tsv", sep="\t",
                       usecols=["user_genome","classification"])
gtdb_arc = pd.read_csv("23_phylogenetics/gtdbtk_full/gtdbtk.ar53.summary.tsv", sep="\t",
                       usecols=["user_genome","classification"])
gtdb = pd.concat([gtdb_bac, gtdb_arc], ignore_index=True)
gtdb_dict = gtdb.set_index("user_genome")["classification"].to_dict()

def get_mag_name(protein_id):
    parts = protein_id.split("__contig")
    return parts[0] if len(parts) > 1 else "unknown"

mag_pass["mag_name"] = mag_pass["protein_id"].apply(get_mag_name)
mag_pass["taxonomy"] = mag_pass["mag_name"].map(gtdb_dict).fillna("unknown")
mag_pass["sample"] = mag_pass["mag_name"].apply(
    lambda x: x.split("__")[0] if "__" in x else x)
mag_pass["source"] = "MAG"

# ── Load taxonkit taxonomy for unbinned contigs ────────────────
# Load taxonomy_lookup to map representative contig → original contig
lookup = pd.read_csv("30_unbinned_contigs/logs/taxonomy_lookup.tsv", sep="\t")
lookup_dedup = lookup.drop_duplicates(subset="representative_id", keep="first")
lookup_dict = lookup_dedup.set_index("representative_id")[["original_contig_id","source_assembly"]].to_dict("index")
print(f"Taxonomy lookup entries: {len(lookup_dict)}")

# Load all taxonkit files into contig → taxonomy dict
taxonkit_files = {
    "WIL1_before": "24_taxonkit_assemblies/WIL1_before_read_GTDB_r220_taxonomy.txt",
    "WIL2_before": "24_taxonkit_assemblies/WIL2_before_read_GTDB_r220_taxonomy.txt",
    "WIL2_after":  "24_taxonkit_assemblies/WIL2_after_read_GTDB_r220_taxonomy.txt",
    "WILB_D1":     "24_taxonkit_assemblies/WILB_D1_read_GTDB_r220_taxonomy.txt",
    "WILB_D3":     "24_taxonkit_assemblies/WILB_D3_read_GTDB_r220_taxonomy.txt",
    "WILB_D10":    "24_taxonkit_assemblies/WILB_D10_read_GTDB_r220_taxonomy.txt",
    "WILD_D0":     "24_taxonkit_assemblies/WILD_D0_read_GTDB_r220_taxonomy.txt",
}

contig_tax = {}
for sample, fpath in taxonkit_files.items():
    if not os.path.exists(fpath):
        continue
    with open(fpath) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                contig_tax[parts[0]] = parts[2]
print(f"Taxonkit entries loaded: {len(contig_tax)}")

def get_unbinned_tax(contig_id):
    # contig_id is representative ID e.g. COASM__contig_15042
    # look up original contig via lookup_dict
    if contig_id in lookup_dict:
        orig = lookup_dict[contig_id]["original_contig_id"]
        return contig_tax.get(orig, "unknown")
    # Try stripping prefix
    bare = contig_id.split("__")[-1] if "__" in contig_id else contig_id
    return contig_tax.get(bare, "unknown")

ub_pass["taxonomy"] = ub_pass["contig_id"].apply(get_unbinned_tax)
ub_pass["sample"] = ub_pass["contig_id"].apply(
    lambda x: x.split("__")[0] if "__" in x else x)
ub_pass["mag_name"] = ub_pass["contig_id"]
ub_pass["source"] = "unbinned"

# ── Combine and save annotation table ────────────────────────
mag_ann = mag_pass[["protein_id","profile","mag_name","sample","source","taxonomy",
                     "pct_identity","query_coverage","evalue","bitscore"]].copy()
ub_ann = ub_pass[["protein_id","profile","mag_name","sample","source","taxonomy",
                   "pct_identity","query_coverage","evalue","bitscore"]].copy()

combined_ann = pd.concat([mag_ann, ub_ann], ignore_index=True)
combined_ann.to_csv(f"{OUTDIR}/pmoA_sdimo_annotations.tsv", sep="\t", index=False)
print(f"\nAnnotation table: {OUTDIR}/pmoA_sdimo_annotations.tsv")

# ── Extract protein sequences ─────────────────────────────────
all_ids = set(combined_ann["protein_id"].tolist())
print(f"\nExtracting {len(all_ids)} protein sequences...")

mag_faa = "21_coassembly/MAGs/hyd_coxl_classification/prodigal_prefixed/all_159_MAGs_prefixed.faa"
ub_faa = "30_unbinned_contigs/seqs/unbinned_dereplicated.faa"

seqs = {}
for faa in [mag_faa, ub_faa]:
    if not os.path.exists(faa):
        print(f"  WARNING: {faa} not found")
        continue
    seq_id = None; seq_lines = []
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                if seq_id and seq_id in all_ids:
                    seqs[seq_id] = "".join(seq_lines)
                seq_id = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if seq_id and seq_id in all_ids:
            seqs[seq_id] = "".join(seq_lines)

print(f"Sequences found: {len(seqs)} / {len(all_ids)}")

ann_dict = combined_ann.set_index("protein_id").to_dict("index")

# ── Write per-profile FASTAs ──────────────────────────────────
for profile_short, profile_names in [
    ("pmoA", ["Particulate_methane_monooxygenase_PmoA",
              "Particulate_methane_monooxygenase_PmoA_Feb2024",
              "pMMO-like_CUMMO_PxmA_Feb2024"]),
    ("sdimo_mmoA", ["Soluble_methane_monooxygenase_MmoA",
                    "Soluble_diiron_monooxygenase_Feb2024"])
]:
    profile_ids = combined_ann[combined_ann["profile"].isin(profile_names)]["protein_id"].tolist()
    out_faa = f"{OUTDIR}/{profile_short}_sequences.faa"
    n = 0
    with open(out_faa, "w") as f:
        for pid in profile_ids:
            if pid not in seqs:
                continue
            ann = ann_dict.get(pid, {})
            mag = ann.get("mag_name", "unknown")
            sample = ann.get("sample", "unknown")
            profile = ann.get("profile", "unknown")
            tax = ann.get("taxonomy", "unknown")
            tax_short = ";".join(str(tax).split(";")[:6]) if tax != "unknown" else "unknown"
            f.write(f">{pid} | {profile} | {mag} | {sample} | {tax_short}\n")
            f.write(seqs[pid] + "\n")
            n += 1
    print(f"{profile_short}: {n} sequences → {out_faa}")

# ── Write combined FASTA ──────────────────────────────────────
out_combined = f"{OUTDIR}/all_pmoA_sdimo_sequences.faa"
with open(out_combined, "w") as f:
    for pid, seq in seqs.items():
        ann = ann_dict.get(pid, {})
        mag = ann.get("mag_name", "unknown")
        sample = ann.get("sample", "unknown")
        profile = ann.get("profile", "unknown")
        tax = ann.get("taxonomy", "unknown")
        tax_short = ";".join(str(tax).split(";")[:6]) if tax != "unknown" else "unknown"
        f.write(f">{pid} | {profile} | {mag} | {sample} | {tax_short}\n")
        f.write(seq + "\n")
print(f"Combined: {len(seqs)} sequences → {out_combined}")

print(f"\nSummary by profile:")
print(combined_ann.groupby("profile").size().to_string())
print(f"\nSummary by sample:")
print(combined_ann.groupby("sample").size().sort_values(ascending=False).to_string())
print(f"\nSummary by taxonomy (top 10):")
combined_ann["phylum"] = combined_ann["taxonomy"].apply(
    lambda x: str(x).split(";")[1] if len(str(x).split(";")) > 1 else "unknown")
print(combined_ann.groupby("phylum").size().sort_values(ascending=False).head(10).to_string())

# ════════════════════════════════════════════════════════════
# Extract Group 1 and 2a hydrogenases + amoA
# ════════════════════════════════════════════════════════════

print("\n" + "="*60)
print("Extracting hydrogenases Group 1 + 2a and amoA")
print("="*60)

hyd_groups = ["[NiFe]_Group_1a","[NiFe]_Group_1c","[NiFe]_Group_1h",
              "[NiFe]_Group_1l","[NiFe]_Group_2a"]
amoa_profiles = ["Ammonia_monooxygenase_amoA","Ammonia_monooxygenase_AmoA_Feb2024"]

# Load hydrogenase results
hyd = pd.read_csv("27_hydrogenase_coxl_final/results/all_hydrogenase_classifications.tsv", sep="\t")

# MAG hydrogenases
hyd_mag = hyd[(hyd["class_flag"]=="PASS") &
              (hyd["data_type"]=="MAG") &
              (hyd["hydrogenase_group"].isin(hyd_groups))].copy()
hyd_mag = hyd_mag.drop_duplicates(subset="protein_id", keep="first")
hyd_mag["mag_name"] = hyd_mag["protein_id"].apply(get_mag_name)
hyd_mag["taxonomy"] = hyd_mag["mag_name"].map(gtdb_dict).fillna("unknown")
hyd_mag["sample"] = hyd_mag["mag_name"].apply(lambda x: x.split("__")[0] if "__" in x else x)
hyd_mag["source"] = "MAG"
hyd_mag = hyd_mag.rename(columns={"hydrogenase_group":"profile"})

# Unbinned hydrogenases
hyd_ub = hyd[(hyd["class_flag"]=="PASS") &
             (hyd["data_type"]=="unbinned") &
             (hyd["hydrogenase_group"].isin(hyd_groups))].copy()
hyd_ub = hyd_ub.drop_duplicates(subset="protein_id", keep="first")
hyd_ub["taxonomy"] = hyd_ub["contig_id"].apply(get_unbinned_tax) if "contig_id" in hyd_ub.columns else "unknown"
hyd_ub["sample"] = hyd_ub["protein_id"].apply(lambda x: x.split("__")[0] if "__" in x else x)
hyd_ub["mag_name"] = hyd_ub["protein_id"].apply(lambda x: "__".join(x.split("__")[:-1]) if "__" in x else x)
hyd_ub["source"] = "unbinned"
hyd_ub = hyd_ub.rename(columns={"hydrogenase_group":"profile"})

print(f"Hydrogenase Group 1+2a MAG: {len(hyd_mag)}")
print(f"Hydrogenase Group 1+2a unbinned: {len(hyd_ub)}")

# amoA from GreeningLab
amoa_mag = greening_mag[(greening_mag["class_flag"]=="PASS") &
               (greening_mag["data_type"]=="MAG") &
               (greening_mag["profile"].isin(amoa_profiles))].copy()
amoa_mag = amoa_mag.drop_duplicates(subset="protein_id", keep="first")
amoa_mag["mag_name"] = amoa_mag["protein_id"].apply(get_mag_name)
amoa_mag["taxonomy"] = amoa_mag["mag_name"].map(gtdb_dict).fillna("unknown")
amoa_mag["sample"] = amoa_mag["mag_name"].apply(lambda x: x.split("__")[0] if "__" in x else x)
amoa_mag["source"] = "MAG"

amoa_ub = greening_ub[(greening_ub["class_flag"]=="PASS") &
             (greening_ub["profile"].isin(amoa_profiles))].copy()
amoa_ub = amoa_ub.drop_duplicates(subset="protein_id", keep="first")
amoa_ub["taxonomy"] = amoa_ub["contig_id"].apply(get_unbinned_tax)
amoa_ub["sample"] = amoa_ub["contig_id"].apply(lambda x: x.split("__")[0] if "__" in x else x)
amoa_ub["mag_name"] = amoa_ub["contig_id"]
amoa_ub["source"] = "unbinned"

print(f"amoA MAG: {len(amoa_mag)}")
print(f"amoA unbinned: {len(amoa_ub)}")

# Collect all new IDs and extract sequences
cols = ["protein_id","profile","mag_name","sample","source","taxonomy","pct_identity","query_coverage","evalue","bitscore"]
new_ann = pd.concat([
    hyd_mag[["protein_id","profile","mag_name","sample","source","taxonomy","pct_identity","evalue","bitscore"]], hyd_ub[["protein_id","profile","mag_name","sample","source","taxonomy","pct_identity","evalue","bitscore"]],
    amoa_mag[cols], amoa_ub[cols]
], ignore_index=True)
new_ann.to_csv(f"{OUTDIR}/hyd_group1_2a_amoa_annotations.tsv", sep="\t", index=False)

new_ids = set(new_ann["protein_id"].tolist())
new_seqs = {}
for faa in [mag_faa, ub_faa]:
    if not os.path.exists(faa): continue
    seq_id = None; seq_lines = []
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                if seq_id and seq_id in new_ids:
                    new_seqs[seq_id] = "".join(seq_lines)
                seq_id = line[1:].strip().split()[0]; seq_lines = []
            else:
                seq_lines.append(line.strip())
        if seq_id and seq_id in new_ids:
            new_seqs[seq_id] = "".join(seq_lines)

print(f"Sequences found: {len(new_seqs)} / {len(new_ids)}")

new_ann_dict = new_ann.set_index("protein_id").to_dict("index")

for label, profiles in [
    ("hyd_group1_2a", hyd_groups),
    ("amoA", amoa_profiles)
]:
    pids = new_ann[new_ann["profile"].isin(profiles)]["protein_id"].tolist()
    out = f"{OUTDIR}/{label}_sequences.faa"
    n = 0
    with open(out, "w") as f:
        for pid in pids:
            if pid not in new_seqs: continue
            a = new_ann_dict.get(pid, {})
            tax = str(a.get("taxonomy","unknown"))
            tax_short = ";".join(tax.split(";")[:6]) if tax != "unknown" else "unknown"
            f.write(f">{pid} | {a.get('profile','?')} | {a.get('mag_name','?')} | {a.get('sample','?')} | {tax_short}\n")
            f.write(new_seqs[pid] + "\n")
            n += 1
    print(f"{label}: {n} sequences → {out}")

print(f"\nGroup breakdown:")
print(new_ann.groupby(["profile","source"]).size().to_string())
