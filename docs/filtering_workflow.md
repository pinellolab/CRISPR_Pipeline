# Filtering and MuData assembly order

The Gary Hon samplesheet contains three measurement-set batches:
`IGVFDS6244NAXC`, `IGVFDS8721BKRO`, and `IGVFDS9613DDRB`. Its `lane` column is
empty, so `measurement_sets` is the key that pairs RNA, guide, and hashing data.
Each batch has one RNA FASTQ pair, two guide FASTQ pairs, and one hashing FASTQ
pair. The two guide pairs are pooled into one mapping call for that batch.

The `lane` field is not the grouping key in this checkout. It is empty for all
Gary Hon rows, and `meta.id` is built as
`sample_${file_modality}_${measurement_sets}`. Consequently, this run has one
mapping tuple per modality and measurement set, not one tuple per lane.

```mermaid
flowchart TD
    INPUT["Gary Hon samplesheet"]
    GROUP["Group by file_modality + measurement_sets"]

    subgraph BATCHES["Three matched measurement-set batches"]
        B1["IGVFDS6244NAXC<br/>RNA: 1 pair<br/>guide: 2 pairs pooled<br/>hash: 1 pair"]
        B2["IGVFDS8721BKRO<br/>RNA: 1 pair<br/>guide: 2 pairs pooled<br/>hash: 1 pair"]
        B3["IGVFDS9613DDRB<br/>RNA: 1 pair<br/>guide: 2 pairs pooled<br/>hash: 1 pair"]
    end

    INPUT --> GROUP
    GROUP --> B1
    GROUP --> B2
    GROUP --> B3

    subgraph RNA["scRNA processing"]
        R1["Map RNA 6244"]
        R2["Map RNA 8721"]
        R3["Map RNA 9613"]
        RC["Concatenate RNA AnnData<br/>barcodes qualified by measurement set"]
        RK{"QC_barcode_filter"}
        RUMI["Knee or knee2 barcode filter<br/>keep cells above total RNA UMI threshold"]
        RGENE["Minimum genes per cell<br/>QC_min_genes_per_cell"]
        RGB["Standard gene prefilter<br/>gene detected in at least 10 cells"]
        RMITO["Mitochondrial filter<br/>percent_mito less than QC_pct_mito"]
        RF["filtered_anndata.h5ad"]

        R1 --> RC
        R2 --> RC
        R3 --> RC
        RC --> RK
        RK -- "knee / knee2" --> RUMI
        RK -- "none" --> RGENE
        RUMI --> RGB
        RGENE --> RGB
        RGB --> RMITO --> RF
    end

    subgraph GUIDE["Guide-count processing"]
        G1["Map pooled guide FASTQs 6244"]
        G2["Map pooled guide FASTQs 8721"]
        G3["Map pooled guide FASTQs 9613"]
        GC["Concatenate guide AnnData<br/>barcodes qualified by measurement set"]

        G1 --> GC
        G2 --> GC
        G3 --> GC
    end

    subgraph HASH["Hashing processing when enabled"]
        H1["Map hash 6244"]
        H2["Map hash 8721"]
        H3["Map hash 9613"]
        HC["Concatenate hashing AnnData<br/>barcodes qualified by measurement set"]
        HI["Intersect qualified hashing keys<br/>with RNA-QC-surviving keys"]
        HS["Split hashing data by batch"]
        HD["GMM-demux each batch<br/>retain filtered demultiplexed cells"]
        HCAT["Concatenate demultiplexed hashing AnnData"]

        H1 --> HC
        H2 --> HC
        H3 --> HC
        HC --> HI --> HS --> HD --> HCAT
    end

    B1 --> R1
    B1 --> G1
    B1 --> H1
    B2 --> R2
    B2 --> G2
    B2 --> H2
    B3 --> R3
    B3 --> G3
    B3 --> H3

    RF -. "surviving RNA barcodes" .-> HI

    CREATE["CreateMuData<br/>intersect cell barcodes across gene + guide<br/>and hashing when enabled"]
    FIRST["mudata.h5mu<br/>single MuData containing all three batches"]

    RF --> CREATE
    GC --> CREATE
    HCAT --> CREATE
    CREATE --> FIRST

    BSPLIT["Split MuData by batch"]
    A1["CLEANSER assignment<br/>IGVFDS6244NAXC"]
    A2["CLEANSER assignment<br/>IGVFDS8721BKRO"]
    A3["CLEANSER assignment<br/>IGVFDS9613DDRB"]
    ACAT["Concatenate assigned MuData objects"]
    GF["Global gene prevalence filter<br/>expression in more than<br/>QC_min_cells_per_gene x total cells"]
    DUAL{"DUAL_GUIDE?"}
    COLLAPSE["Collapse assigned guides<br/>to intended-target elements"]
    FINAL["concat_mudata.h5mu<br/>input to inference"]

    FIRST --> BSPLIT
    BSPLIT --> A1
    BSPLIT --> A2
    BSPLIT --> A3
    A1 --> ACAT
    A2 --> ACAT
    A3 --> ACAT
    ACAT --> GF --> DUAL
    DUAL -- "false" --> FINAL
    DUAL -- "true" --> COLLAPSE --> FINAL
```

## Verified Nextflow channel lineage

```mermaid
flowchart TD
    CSV["samplesheet CSV"] --> CS["ch_samplesheet: groupTuple by meta.id"]
    CS --> C["ch_samples"]
    C --> CR["ch_rna: 3 tuples"]
    C --> CG["ch_guide: 3 tuples"]
    C --> CH["ch_hash: 3 tuples"]

    CR --> MR["mappingscRNA x3"]
    CG --> MG["mappingGuide x3; 4 FASTQs per task"]
    CH --> MH["mappingHashing x3"]

    MR --> RC["collect/sort -> concat_anndata_rna"]
    MG --> GC["collect/sort -> concat_anndata_guide"]
    MH --> HC["collect/sort -> concat_anndata_hashing"]

    RC --> PRE["PreprocessAnnData -> filtered_anndata_rna"]
    PRE --> FH["filter_hashing"]
    HC --> FH
    FH --> DM["three files -> demultiplex x3"]
    DM --> HCC["collect/sort -> hashing_concat"]

    PRE --> CMD["CreateMuData"]
    GC --> CMD
    HCC --> CMD
    CMD --> M1["mudata.h5mu"]
    M1 --> PA["prepare_assignment -> three batch MuData files"]
    PA --> CL["CLEANSER x3"]
    CL --> MC["collect/sort -> mudata_concat"]
    MC --> M2["concat_mudata.h5mu: 69,123 cells"]
    M2 --> INF["inference_pipeline"]
```

### Evidence from the completed May 14, 2026 run

| Stage | Observed result |
| --- | --- |
| Guide mapping | The 6244 `mappingGuide` task passed four FASTQs to one `kb count`, proving that its two guide pairs were pooled. |
| Modality concatenation | RNA, guide, and hash logs all sorted measurement sets as 6244, 8721, 9613. |
| RNA QC | `PreprocessAnnData` started with 2,131,101 cells and wrote `filtered_anndata.h5ad`. |
| RNA/hash intersection | `filter_hashing` found 128,814 shared keys, then wrote one file for each measurement set. |
| Assignment split | `prepare_assignment` wrote exactly three files: one each for 6244, 8721, and 9613. |
| Final assignment concatenation | `mudata_concat` combined three files and wrote `concat_mudata.h5mu` with 69,123 cells. |

The evidence comes from the checkout at commit `47a36ce`, its `.nextflow.log`,
and task objects under `gs://igvf-pertub-seq-pipeline-data/work`.

## Important details

- Samplesheet rows are grouped by `file_modality + measurement_sets`.
- For each Gary Hon batch, the two guide rows share the same grouping key, so
  four guide FASTQ files (two R1/R2 pairs) enter one guide mapping call.
- RNA, guide, and hashing AnnData objects are physically concatenated across
  batches, but cell keys remain batch-qualified. In the updated workspace, for example,
  `AAAC..._IGVFDS6244NAXC` cannot intersect
  `AAAC..._IGVFDS8721BKRO`.
- The completed May 14, 2026 checkout used positional suffixes `_0`, `_1`, and `_2`.
  Its intersections remained same-batch because all three modality
  concatenations used the same sorted measurement-set order. The updated code
  removes that ordering dependency by using the measurement-set ID directly.
- With `QC_barcode_filter = 'knee'` or `'knee2'`, total RNA UMI depth selects
  cells and `QC_min_genes_per_cell` is skipped. With
  `QC_barcode_filter = 'none'`, the minimum-gene filter is used instead.
- The guide count matrix is not filtered by a fixed UMI cutoff before MuData
  creation. SCEPTRE or CLEANSER assigns guides from the per-cell guide counts
  after the modalities have been intersected.
- When hashing is enabled, hashing is first restricted to cells that survived
  RNA QC. `CreateMuData` then performs the final barcode intersection across
  the gene, guide, and demultiplexed hashing modalities.
- When hashing is disabled, the final intersection uses only gene and guide
  barcodes. Optional Scrublet doublet removal runs on the resulting MuData
  before guide assignment.
- `mudata.h5mu` is one combined multimodal object containing all three batches.
  It is temporarily split by `batch` (`measurement_sets`) so guide assignment
  runs independently, then recombined as `concat_mudata.h5mu`.
- `QC_min_cells_per_gene` is a fraction in `[0, 1)` and is applied after guide
  assignment and batch concatenation. A gene must be detected in strictly more
  than `QC_min_cells_per_gene * total_cells`; `0` retains every gene detected
  in at least one cell.
- `TAPSEQ_QC_MODE = true` removes the standard 10-cell preprocessing floor so
  observed TAP-seq genes reach the final fractional filter. Pair it with a very
  small fraction (for example `0.000001`) when all observed genes should be
  retained.
