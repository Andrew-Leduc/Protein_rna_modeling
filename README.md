# Single-Cell Protein-mRNA Modeling

Statistical modeling of the relationship between mRNA and protein abundances at single-cell resolution using XX distributions and log-normal noise models.

## Overview

This project models protein abundance from mRNA counts in single cells by:
1. Fitting negative binomial (NB) distributions to mRNA count data
2. Simulating protein levels using a log-normal multiplicative noise model
3. Jointly optimizing mRNA and protein parameters to match observed distributions

The core model assumes:
```
P = (c * m + ε) * LogNormal(0, σ)
```
where `m` is mRNA count (NB distributed), `c` is a scaling factor, `ε` prevents log(0), and `σ` captures biological/technical noise.

## Repository Structure

```
├── models/
│   └── Model_V1.ipynb          # Main analysis notebook (Google Colab)
├── preprocessing/
│   └── preproc_and_filter.R    # R preprocessing script (Seurat-based)
├── dat/                        # Local data directory (optional)
└── requirements.txt            # Python dependencies
```

## Data Access

Raw sequencing data is available on GEO:
**[GSE244215](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244215)**

Processed data for modeling is stored on Google Drive:
**[Data Folder](https://drive.google.com/drive/folders/1gv9as_seGLJxYvc90Ogdz1t5bHj6IKIz?usp=sharing)**

### Recommended Google Drive Structure

```
Drive/
└── Protein_RNA_Modeling/
    ├── raw/           # Original, unmodified data
    ├── processed/     # Cleaned/filtered protein data for our modeling purposes
    ├── results/       # Model outputs and figures
    └── notebooks/     # Working copies of Colab notebooks
```

### Notebook Links

- **GitHub source** links (used by the Colab badge) will break if files are moved/renamed in the repo — update the badge URL if you reorganize.
- **Google Drive** links are based on file ID, not path — notebooks can be moved freely within Drive without breaking share links.

### Data Loading (Built into Notebook)

The notebook includes a data loading cell that mounts Google Drive and loads the data. Just update the `DATA_PATH` variable to point to your data folder.

## Quick Start

1. Open [Model_V1.ipynb](https://colab.research.google.com/github/Andrew-Leduc/Protein_rna_modeling/blob/main/models/Model_V1.ipynb) in Google Colab
2. Run the first cell to mount Google Drive
3. Update `DATA_PATH` to point to your data folder
4. Run cells sequentially

## Core Functions

| Function | Description |
|----------|-------------|
| `fit_nb_mle(x)` | Fit negative binomial to mRNA counts via MLE |
| `simulate_protein_log2fc_from_mrna_nb()` | Simulate protein log2FC from mRNA model |
| `fit_sigma_log_to_protein()` | Grid search for optimal noise parameter σ |
| `joint_fit_nb_and_sigma()` | Joint optimization of NB and noise parameters |

## Dependencies

**Python (Colab):**
- numpy
- scipy
- matplotlib

**R (preprocessing):**
- Seurat
- stringr
- seqinr
- dplyr

## Usage Example

```python
# Fit NB to mRNA counts
mu_hat, r_hat, _ = fit_nb_mle(x_mrna)

# Find optimal sigma to match observed protein distribution
best_sigma = fit_sigma_log_to_protein(y_obs, mu_hat, r_hat)

# Or jointly optimize all parameters
result = joint_fit_nb_and_sigma(
    x_mrna,
    y_obs,
    t_half_m_hours=2.0,   # mRNA half-life
    t_half_p_hours=24.0   # protein half-life
)
```

## License

MIT
