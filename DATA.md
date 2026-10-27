# DATA: sources, versions, checksums, licences

Everything under `data/` and `starting_data/` is git-ignored and must not be committed (the PIE data
are CC BY-NC-SA 4.0, and derived caches inherit that licence).
To rebuild everything except the manually provided inputs, run `scripts/data/run_all.sh` (about 3 min plus downloads).
`scripts/data/checksum_inputs.py` regenerates the full checksum table `docs/data_checksums.tsv`
(sha256 + md5 + size for every input file). The tables below are copied from that file on 2026-10-27.
Numbers describing the data are in `docs/DATA_PROFILE.md`.

## 1. DepMap Public 26Q3: `starting_data/` (primary CN + expression + metadata)

| | |
|---|---|
| Source | DepMap portal, https://depmap.org/portal/data_page/?tab=allData (files provided by the PI in `starting_data/`) |
| Release | **DepMap Public 26Q3** (portal announcement 2026-10-05, https://forum.depmap.org/t/announcing-the-26q3-release/4692). This is *inferred*, not md5-verified; see below |
| Licence | DepMap data are released under CC BY 4.0 (the licence stated on DepMap's figshare+ deposits); the DepMap portal terms of use also apply |
| Obtained | Manual browser download (by the PI) |

How the release was identified (`scripts/data/identify_release.py` → `data/processed/release_evidence.json`):
- The authoritative check would match each md5 against the portal listing `https://depmap.org/portal/api/download/files`. That endpoint returns a
  human-verification page and asks not to be scraped, so it was **not** used. Figshare+ deposits
  (md5-comparable) stop at 24Q4, and the 24Q4 md5s of Model/ModelCondition/OmicsProfiles do **not** match.
- Markers in the files:
  - The omics matrices carry `SequencingID, ModelConditionID, ModelID, IsDefaultEntryForMC, IsDefaultEntryForModel` (introduced in 25Q3).
  - `SubtypeTree.csv` exists (25Q3+).
  - SubtypeTree has PAM50 breast subtype nodes ("new breast cancer subtypes", 26Q3 release note).
  - The latest RNA sequencing date is 2025-10-15.
  - 290/291 NextGen-manuscript models with WGS have a WGS profile here (these models were integrated in 26Q3).
- **Action for the PI (optional):** in a browser, compare these md5s with the 26Q3 listing on the portal.

| file | bytes | sha256 | md5 |
|---|---:|---|---|
| `starting_data/Model.csv` | 789,805 | `65eb1895800e71a963a3ab56e5d08c8d006cb0d9e7aa68f7e34a3810a2a3cbc7` | `a7423dbf345fcf5a822ab5fc7e854e5c` |
| `starting_data/ModelCondition.csv` | 314,820 | `a1edda78cd60718e5c99690ce8ab80315b65ddc9346643338d2c1c774fc51326` | `840237c503a47dd32859c4fbf3af2f95` |
| `starting_data/OmicsCNGeneWGS.csv` | 661,789,195 | `31cace66f202adaace306d975b5661c8b77ea82f07d2ed34045eec9c7958fcd6` | `8e5158c0b8c929cf7574f0ae7a7d922a` |
| `starting_data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 394,997,160 | `252d1fedc7014953b19ffe5d68151b4df973c8078b12c9e5ccd97e292548d177` | `b45585365eb9580bcf4618791ea87f76` |
| `starting_data/OmicsProfiles.csv` | 882,568 | `4241e537f0b23139c0e525b4d78fb53ee65c5507771a33e2cd5debd78f7179ea` | `6a025f0c96c66e1cc49086ea9b2f7687` |
| `starting_data/SubtypeTree.csv` | 138,880 | `67cbe65fe3f7cf08efb85b065d93eb98a80d8173839493d961277aa9fb8fedf2` | `a8283200302dbba0292fe9af4dda461f` |

### 1b. 26Q3 covariate files that are NOT available programmatically

Covariates such as ploidy, WGD, MSI score, CIN and LoH (`OmicsGlobalSignatures.csv`), and arm-level CNA calls (if a
26Q3 file exists), can only be downloaded through the bot-protected portal. If the PI downloads them
in a browser into `data/external/depmap_26Q3/`, `download_depmap_extra.py` verifies them and records their checksums.
**Note:** since 25Q2, DepMap has frozen PureCN-derived ploidy, WGD, LoH, absolute CN and arm-level calls
"for samples released till 2024 Q4" (25Q2 release notes), and OmicsAbsoluteCNGene was discontinued after 24Q4 (DepMap
forum, 2026-08). So even the 26Q3 portal files would not cover newer models; §2 holds the same
information for older models.

## 2. DepMap Public 24Q4 (figshare+): covariates frozen since 24Q4

| | |
|---|---|
| Source | https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248 (DOI 10.25452/figshare.plus.27993248.v1) |
| Licence | CC BY 4.0 |
| Obtained | `scripts/data/download_depmap_extra.py` (figshare API; each file's md5 is checked against figshare `supplied_md5`) |
| Use | Model- and profile-level Ploidy, WGD, CIN, LoHFraction, Aneuploidy (PureCN) and MSIScore (MSIsensor2); PureCN absolute CN (used only to verify the 26Q3 CN scale) |

| file | bytes | sha256 | md5 |
|---|---:|---|---|
| `data/external/depmap_24Q4/README.txt` | 43,103 | `f744a5cf5c112e01276e772679157aac8e85cdc4a8f6d3b86d808e3bfd332891` | `54628c15a10b9ff6536db245cfed5231` |
| `data/external/depmap_24Q4/OmicsSignatures.csv` | 131,912 | `2eb6d83065db5b591d9904339a1eaeb769a4db8593fe80a135399d7708b39bef` | `14bdb13656a584dab52e2959a1e11c5a` |
| `data/external/depmap_24Q4/OmicsSignaturesProfile.csv` | 194,252 | `7412ec6b0c607c7174a29e449a8f5b13369782c4fd8a74ddb490cfe0d7119abd` | `94f8afe88b5a896126b4f122d836bb44` |
| `data/external/depmap_24Q4/OmicsDefaultModelProfiles.csv` | 90,080 | `096a96c39d88374cb7c37058816e6cde29d6617b820b8a14d33ae9dd821c6dc5` | `4d342a03701c0a4c0742fedb56c1de40` |
| `data/external/depmap_24Q4/OmicsProfiles.csv` | 254,733 | `0bb19ee57935d4b213d378a4db9caacaa1db7318a63ae81485d84f91c20a3bf2` | `7195088a5fc41eddb3f896c6e984a992` |
| `data/external/depmap_24Q4/Model.csv` | 645,696 | `b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa` | `675210d17675f3517b0ce39a3c274f16` |
| `data/external/depmap_24Q4/OmicsAbsoluteCNSegmentsProfile.csv` | 33,904,573 | `5bb6b52550c372a929f9bd8f09f13344833472ef61967814022cac37c856fced` | `e5d3a4425a3ce6d77c42ca75cd13b71b` |
| `data/external/depmap_24Q4/OmicsAbsoluteCNGene.csv` | 238,808,692 | `5f530c797bf12bce68cccc3e763b29491c368bbbced8971733bed807a13b8c13` | `16afffb33230de2fd13e15e753e28a9b` |

Linkage to 26Q3: a signature is taken from `OmicsSignaturesProfile` when its ProfileID equals the 26Q3 default
WGS ProfileID (`sig24_source = same_profile`), otherwise from model-level `OmicsSignatures`
(`model_level_24Q4`, which may come from WES). See `data/processed/depmap_models.tsv`.

**Not available anywhere checked:**
- **Purity:** DepMap ships none. For cell lines it is ~1 by construction, and PureCN purity is not in the public files. A possible
  alternative is CCLE 2019 ABSOLUTE purity (Ghandi et al. 2019 supplement, also portal-gated; it covers only older lines).
- **Arm-level CNA for 26Q3 WGS profiles:** derived here from gene-level CN (`data/processed/depmap_arm_cn_derived.parquet`).

## 3. Sanger Cell Model Passports model list (independent MSI / ploidy)

| | |
|---|---|
| Source | https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz (the gzip header mtime is 2021-04-27, so the file is old) |
| Licence | Cell Model Passports terms (https://cellmodelpassports.sanger.ac.uk); licence text not verified, check before redistribution |
| Obtained | `scripts/data/download_depmap_extra.py`; linked to DepMap through `Model.csv` `SangerModelID` |
| File | `data/external/sanger_cmp/model_list_latest.csv.gz`, 145,029 B, sha256 `333affaa690c9f7e517467a6f78b00ffb418cb8449e2eb1cc5d93d85f7b0862e`, md5 `28f53085bf4db7b855bd8b251ca7d841` |

## 4. PIE X-Atlas/Orion DE tables (benchmark)

| | |
|---|---|
| Source | https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion |
| Version | HF commit `abfa3da0f3bb98c9df9c90082386ed0c7c9282d3` |
| Underlying data | X-Atlas/Orion genome-wide CRISPRi Perturb-seq (Huang et al. 2025), HCT116 (CVCL_0291) and HEK293T (CVCL_0063). DE was computed by PIE `pie process` (gpudge 0.7.0, Mann-Whitney U vs non-targeting, BH; we verified BH is applied within each target) |
| Licence | **CC BY-NC-SA 4.0** (`DATA_LICENSE.md`): non-commercial; derived files share the licence; never commit |
| Obtained | Downloaded by the PI in M0 with `huggingface_hub` (local-dir download pinned to the commit; `.cache/huggingface/download/*.metadata` records commit + etag). Present: `de/*.parquet`, `preprocessed/{contexts.yaml,meta.json}`, `README.md`, `DATA_LICENSE.md`; the large `preprocessed/*.npy` arrays were deliberately not downloaded. Only the pre-computed DE tables are used |
| Verification | sha256 of both parquet files equals the HF LFS etag (`checksum_inputs.py`) |

| file | bytes | sha256 | md5 |
|---|---:|---|---|
| `de/hct116.parquet` | 3,395,023,058 | `0858718bc0051835874920a64db804c80cba30c7b70a24b0852643da153d867f` | `1baad538111bf78dc17d375f5f7b9737` |
| `de/hek293t.parquet` | 1,820,568,723 | `02924904206e943ac616c35544145cabc7b62c2ec3de1b74c8255bb3de8c61c0` | `f941022efd177e7f3a7624da709aa05b` |
| `preprocessed/meta.json` | 473,417 | `16b83a176823c2e87b5ccdfb6fd2cc8283e60f860e1797acdc9d01e6aa77398c` | `022a3ccebc665615c4c66c3e0f094bf0` |
| `preprocessed/contexts.yaml` | 129 | `13b724c3de36980ddc200f34fbd424acbd1e4ac55bfa89177ee395927d40bbc3` | `5b0922e0e039177c1a2fd5c6c3bdf9a4` |
| `README.md` | 4,988 | `22608de3f8928593a526aa5e935cfa379a77ffbbc7c64544a2c9bb879d1eb182` | `fd8c36507de9cd73cca1efe13a65db05` |
| `DATA_LICENSE.md` | 18,859 | `8dad297919c8c508c5cae28aac6fa059e34daf66071ef0f30f3d1528b262c86a` | `684c5b22fe35557a27fff38eca6d21c3` |

## 5. Gene annotation resources

| dataset | source / version | licence | file, sha256 |
|---|---|---|---|
| HGNC protein-coding set (`starting_data/gene_with_protein_product.txt`) | HGNC `gene_with_protein_product.txt`, provided by the PI; latest `date_modified` 2026-10-23 | CC0 | 9,498,389 B, `6a274291fd9029a2ca41b8b00cd114b039629db3a6d68a52eedd3e6b422a3304` (md5 `b0d895d293f9ae43c4018eed73b2ac61`) |
| HGNC complete set (all locus types; needed for non-coding Perturb-seq features) | https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-10-06.tsv (monthly archive, pinned) | CC0 | 16,973,125 B, `8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d`; md5 matches the GCS object md5 |
| Ensembl GTF, GRCh38 | https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz (**Ensembl 116**, current at download) | Ensembl: no restrictions on use | 141,121,632 B, `ed992f0eac7197d9627bda618f8f831ba355c95bd5d0796af785387d462828b6`; BSD `sum` matches Ensembl CHECKSUMS |
| UCSC hg38 cytoBand | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz (Last-Modified 2022-10-29) | UCSC Genome Browser data, free use | 12,118 B, `e514e48b0a12a516fd5231b44241a20a282f290bdd6131cadcddc6276ee26082` |

All of these are obtained by `scripts/data/download_annotation.py`, except the PI-provided HGNC protein-coding file.

## 6. NextGen 3D model metadata (used only for release identification)

https://figshare.com/articles/dataset/29472362 (`model_metadata.csv`, figshare file 64494438), CC BY 4.0; downloaded by
`identify_release.py` to `data/external/nextgen_2026/`; 272,338 B, sha256
`a000a5ee1ebe2b2488a146b2a413ec4b93472488e6ad31b56126f522d652557c`, md5 `97b26ec52b8a05e52df7abd507b8139f` (asserted).

## Derived files (`data/processed/`, `data/interim/`; not committed)

| file | script | content |
|---|---|---|
| `data/interim/depmap_{cn_wgs,expr_tpm_logp1}.parquet` | `convert_depmap.py` | all profile rows, float32 (max abs rounding error 7e-6) |
| `gene_annotation.tsv` | `build_gene_annotation.py` | one row per DepMap gene column: HGNC id/symbol, Entrez, Ensembl, chrom, arm, cytoband, GRCh38 start/end (Ensembl 116), chrom_class, PAR flag |
| `pie_gene_map.tsv` | `build_gene_annotation.py` | one row per PIE target/feature identifier: mapping method, HGNC, location, and the DepMap CN/expression column it links to (one-to-one) |
| `depmap_models.tsv` | `build_depmap_tables.py` | one row per model with a default CN and/or RNA profile: profile ids, Model/ModelCondition annotations, 24Q4 signatures (with provenance), CMP MSI/ploidy |
| `depmap_{cn,expr}_matched.parquet` | `build_depmap_tables.py` | 1,724 matched models, identical row order (sorted by ModelID) |
| `depmap_arm_cn_derived.parquet` | `build_depmap_tables.py` | per model × arm, median/mean log2 relative CN (floor 2^-5) |
| `depmap_cn_variability.tsv`, `depmap_cis_spearman.tsv`, `depmap_profile.json` | `profile_depmap.py` | per-gene CN variability statistics, cis CN–expression correlation, profile numbers |
| `pie/tested_<cell>.parquet` | `cache_pie_de.py` | **all** tested (target, feature) pairs: uint16 codes, float32 lfc, −log10 p, −log10 p_adj |
| `pie/loose_<cell>.parquet` | `cache_pie_de.py` | full-precision rows with p < 1e-3 or p_adj < 0.1 |
| `pie/{vocab_*,targets_*,features_*,target_summary_*}.tsv` | `cache_pie_de.py`, `profile_benchmark.py` | code vocabularies; per-target cell counts, on-target LFC and hit counts; per-feature control mean |
| `benchmark_profile.json`, `release_evidence.json` | `profile_benchmark.py`, `identify_release.py` | profile numbers |
