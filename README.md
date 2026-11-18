# Depth-partitioned oxidation of atmospheric hydrogen, carbon monoxide and methane in a glacier forefield soil analysis scripts


Key scripts used for the manuscript (Herbert et al., 2026 "Depth-partitioned oxidation of atmospheric hydrogen, carbon monoxide and methane in a glacier forefield soil", under submission).
These are the SLURM / shell / Python scripts that were run on a SLURM HPC cluster.
Scripts are grouped by the section of the Methods they correspond to, and numbered in the order they were run.
Per-sample re-runs, test runs, logs and merged intermediate tables are not included.

Raw reads: NCBI BioProject PRJNA1559487.

| Folder | Methods section | Scripts | Tools |
|---|---|---|---|
| `01_read_QC_and_assembly` | Metagenomic data processing | `01_barcoding_fastq_processing.sh`, `02_qc_trim_filter.sh`, `03_flye_assembly.sh`, `coassembly_flye.sh`, `assembly_stats.sh` | Chopper, Porechop, Flye (per-sample and co-assembly) |
| `02_taxonomic_profiling_reads` | Metagenomic data processing | `04_kraken2_bracken_gtdb.sh`, `07_bracken.sh` | Kraken2, Bracken (GTDB r220) |
| `03_MAG_recovery` | MAG recovery | `run_pipeline.sh` (master), `01_read_mapping.sh`, `02_binning.sh`, `03_dastool.sh`, `04_drep_95.sh`, `05_gtdbtk_coverm.sh`, `21_coassembly_MAG_pipeline.sh`, `coassembly_checkm2_gtdb_drep.sh`, `coverm_159_MAGs_metaG.sh` | minimap2, MetaBAT2, CONCOCT, MaxBin2, DAS Tool, CheckM2, dRep (95% ANI), GTDB-Tk, CoverM |
| `04_unbinned_contigs` | Unbinned contigs | `unbinned_contig_pipeline.sh`, `coasm_full_prodigal.sh`, `coverm_metaG_assemblies.sh` | minimap2, MMseqs2, Prodigal, CoverM |
| `05_marker_genes_and_hydrogenases` | Functional annotation | `build_greening_dmnd.sh`, `greening_diamond_classify.sh`, `hydrogenase_classify_final.sh` | DIAMOND, HMMsearch (GreeningLab database, HydDB) |
| `06_gene_phylogenies` | Functional annotation (trees) | `extract_coxL_tree.py`, `extract_pmoa_sdimo_seqs.py` | Extraction of coxL and pmoA/amoA/SDIMO sequences for the trees |
| `07_MGEs_viruses_plasmids` | MGE identification and activity | `06_genomad.sh`, `07_checkv_votu.sh`, `08_plasmid_potu.sh`, `09_predict_proteins.sh`, `iphop_add_mags_v141.sh`, `17_coverm_metaT_vOTU_pOTU.sh`, `15_transcriptm_vOTUs_pOTUs.sh` | geNomad, CheckV, dRep, Prodigal, iPHoP, CoverM, TranscriptM |
| `08_metatranscriptomics` | Metatranscriptomic processing | `15_transcriptm_MAGs.sh`, `15b_transcriptm_assemblies.sh`, `rename_assembly_contigs_transcriptm.sh`, `coasm_transcriptm.sh`, `coasm_transcriptm_markers.sh`, `17_coverm_MAGs.sh`, `17_coverm_metaT_assemblies.sh`, `coverm_metaT_coasm_mags.sh` | TranscriptM, CoverM |

## Placeholders
Cluster-specific paths have been replaced with bracketed placeholders. Replace these with your own paths before running.

| Placeholder | Meaning |
|---|---|
| `[reads_dir]` | Directory of the sequencing run (trimmed and filtered ONT metagenomic reads under `fastq_trimmed/trimmed_filtered`) |
| `[reads_dir_run1]`, `[reads_dir_run2]`, `[reads_dir_run3]` | Directories of the separate ONT sequencing runs, used where one script reads samples from several runs |
| `[metaT_reads_dir]` | Directory of the Illumina metatranscriptome reads |
| `[assemblies_dir]` | Directory holding the per-sample Flye assemblies |
| `[MAG_dir]` | Working directory for the MAG-recovery steps |
| `[conda_base]` | Base directory of your conda installation |
| `[software_dir]` | Directory of locally installed software and virtual environments |
| `[kraken2_gtdb_db]` | Kraken2/Bracken GTDB r220 database |
| `[gtdbtk_db]` | GTDB-Tk reference data (GTDB r220) |
| `[genomad_db]`, `[genomad_venv]` | geNomad database and virtual environment |
| `[iphop_db]` | iPHoP database (v1.4.1, June 2025 release) |
| `[transcriptm_db]` | TranscriptM database |
| `[contamination_db]` | ONT lambda DCS control sequence used for read decontamination |
| `[slurm_account]` | Your SLURM account |

## Notes
- The MAG set was built in stages (per-sample bins, additional bins from one sample, then a co-assembly of all samples), and results were compiled afterwards.
  All 159 MAGs are listed in Supplementary Table S5.
- Conda environment names inside the scripts are those used on the original system and may need recreating.
- Databases and software versions are given in the Methods of the manuscript.
- Phylogenetic trees (IQ-TREE v3.0.1, `-MFP -bb 1000`) are available on Zenodo along with the sequences used to build them.
