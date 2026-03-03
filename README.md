# Proteomics DIA Analysis Pipeline

An end-to-end bioinformatics pipeline for label-free quantitative (LFQ) proteomics analysis — from MaxQuant output to pathway-level biological interpretation.

This project was built to practice and demonstrate a real proteomics workflow using a publicly available neuronal dataset. The dataset comes from SH-SY5Y human neuroblastoma cells, a widely used model in neurodegeneration research, comparing three differentiation states: undifferentiated, RA-differentiated, and RA+PMA-differentiated.

---

## Why I Built This

I wanted to work through a complete proteomics analysis independently — not just run tools, but understand every step: why we impute missing values the way we do, why limma outperforms a simple t-test at low sample sizes, and how pathway enrichment connects protein-level changes to biology. The dataset was chosen specifically because the findings (cytoskeletal remodeling, synaptic pathways, axon guidance) are directly relevant to neuronal maintenance and neuromuscular research.

---

## What the Pipeline Does

```
data/proteinGroups.txt
        │
        ▼
scripts/proteomics_analysis.py
  ├── Filter contaminants and reverse hits
  ├── Log2 transformation of LFQ intensities
  ├── Perseus-style downshift imputation for missing values
  ├── PCA to check sample clustering
  ├── Welch's t-test + Benjamini-Hochberg FDR correction
  └── Volcano plots and heatmaps
        │
        ▼
scripts/limma_analysis.R
  ├── Empirical Bayes linear modeling (limma)
  ├── Two contrasts: RA vs Undiff, PMA vs Undiff
  └── Volcano plots
        │
        ▼
scripts/pathway_enrichment.py
  ├── Over-representation analysis via Enrichr (gseapy)
  ├── GO Biological Process, GO Cellular Component, KEGG, Reactome
  ├── Dotplots of top enriched terms
  └── Filtered output highlighting neuronal and cytoskeletal pathways
```

---

## Repository Structure

```
proteomics_dia_analysis_pipeline/
├── data/
│   └── proteinGroups.txt                        # MaxQuant output — not tracked by git
├── env/
│   └── environment.yml                          # Conda environment
├── scripts/
│   ├── proteomics_analysis.py                   # Preprocessing, PCA, t-test, plots
│   ├── limma_analysis.R                         # Differential expression with limma
│   └── pathway_enrichment.py                   # GO/KEGG/Reactome enrichment
├── results/
│   ├── pca.png
│   ├── volcano_undiff_vs_RA.png
│   ├── volcano_undiff_vs_PMA.png
│   ├── volcano_limma_RA_vs_Undiff.png
│   ├── volcano_limma_PMA_vs_Undiff.png
│   ├── heatmap_undiff_vs_RA.png
│   ├── heatmap_undiff_vs_PMA.png
│   ├── MA_plot_limma.png
│   ├── significant_proteins_undiff_vs_RA.csv
│   ├── significant_proteins_undiff_vs_PMA.csv
│   ├── significant_proteins_limma_RA_vs_Undiff.csv
│   ├── significant_proteins_limma_PMA_vs_Undiff.csv
│   ├── all_proteins_undiff_vs_RA.csv
│   ├── all_proteins_undiff_vs_PMA.csv
│   ├── all_proteins_limma_RA_vs_Undiff.csv
│   ├── all_proteins_limma_PMA_vs_Undiff.csv
│   └── enrichment/
│       ├── dotplot_RA_vs_Undiff.png
│       ├── dotplot_PMA_vs_Undiff.png
│       ├── enrichment_RA_vs_Undiff.csv
│       ├── enrichment_PMA_vs_Undiff.csv
│       ├── DLR_relevant_pathways_RA_vs_Undiff.csv
│       └── DLR_relevant_pathways_PMA_vs_Undiff.csv
├── .gitignore
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Getting Started

### Clone the repository

```bash
git clone https://github.com/Kin-zala/proteomics_dia_analysis_pipeline.git
cd proteomics_dia_analysis_pipeline
```

### Set up the Python environment

Using conda (recommended):

```bash
conda env create -f env/environment.yml
conda activate proteomics_env
```

Or using pip:

```bash
pip install -r requirements.txt
```

### Set up R dependencies

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("limma")
install.packages(c("dplyr", "ggplot2"))
```

### Get the data

Download `proteinGroups.txt` from [PRIDE PXD031054](https://www.ebi.ac.uk/pride/archive/projects/PXD031054) and place it in the `data/` folder. The file is not tracked by git due to its size.

---

## Running the Pipeline

Run the three scripts in order:

**Step 1 — Preprocessing and differential analysis**

```bash
python scripts/proteomics_analysis.py
```

**Step 2 — Limma differential expression**

```bash
Rscript scripts/limma_analysis.R
```

**Step 3 — Pathway enrichment**

```bash
python scripts/pathway_enrichment.py
```

All outputs are saved automatically to the `results/` folder.

---

## Results Summary

| Comparison | Significant proteins (t-test) | Significant proteins (limma) |
|---|---|---|
| Undiff vs RA | 76 | 118 |
| Undiff vs PMA | 143 | 162 |

Limma consistently finds more significant hits than the t-test because it borrows variance information across all proteins (empirical Bayes smoothing), making it more sensitive — especially with n=6 per group.

**Key pathways enriched (neuronal and cytoskeletal relevance):**

| Pathway | Comparison | FDR |
|---|---|---|
| Cytoskeleton (GO:0005856) | Both | < 0.02 |
| Actin Cytoskeleton (GO:0015629) | PMA | < 0.02 |
| Regulation of actin cytoskeleton | Both | < 0.03 |
| Axon Guidance | Both | < 0.03 |
| GABAergic / Glutamatergic / Cholinergic synapse | RA | < 0.05 |
| Muscle Contraction | RA | < 0.04 |
| ATR Replication Stress Response | RA | < 0.01 |

The stronger cytoskeletal signal in the PMA condition compared to RA alone is consistent with RA+PMA being a more complete differentiation protocol — something visible both in the protein counts and the enrichment FDR values.

---

## Methods

| Step | Method |
|---|---|
| Filtering | Remove reverse hits, contaminants, only-identified-by-site |
| Normalization | Log2 transformation of LFQ intensities |
| Missing value filter | Minimum 50% valid values per protein across all samples |
| Imputation | Perseus-style downshift (shift = 1.8 SD, width = 0.3 SD) |
| Differential expression | Welch's t-test + BH FDR (Python); limma eBayes (R) |
| Significance threshold | FDR < 0.05 and \|log2FC\| > 1 |
| Enrichment | Over-representation analysis via Enrichr (gseapy) |
| Databases | GO BP, GO CC, KEGG 2021, Reactome 2022 |

---

## Dataset

| Field | Detail |
|---|---|
| Source | PRIDE Archive |
| Accession | [PXD031054](https://www.ebi.ac.uk/pride/archive/projects/PXD031054) |
| Organism | *Homo sapiens* |
| Cell line | SH-SY5Y neuroblastoma |
| Conditions | Undifferentiated, RA, RA+PMA |
| Replicates | n = 6 per condition |
| Format | MaxQuant proteinGroups.txt (LFQ) |

---

## Author

**Kinnari Zala**  
MSc Bioinformatics  
[github.com/Kin-zala](https://github.com/Kin-zala)

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.