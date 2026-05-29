from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURE_DIR = PROJECT_ROOT / "figures"
REPORTS_DIR = PROJECT_ROOT / "reports"

# Analysis-ready outputs.
WIDE_OUTPUT_FILE = PROCESSED_DIR / "checkpoint_expression_wide.csv"

def main():
    wide = pd.read_csv(WIDE_OUTPUT_FILE)

    gene_matrix = wide.set_index("gene")
    gene_correlation = gene_matrix.T.corr()

    plt.figure(figsize=(10, 8))
    mask = np.triu(np.ones_like(gene_correlation, dtype=bool))
    sns.heatmap(gene_correlation, annot=True, cmap="coolwarm", vmin=-1, vmax=1, mask=mask, linewidths=0.5, cbar_kws={"shrink": 0.8})
    plt.title("Correlation of Gene Expression Across Cell Types")

    plt.savefig(FIGURE_DIR / "gene_gene_correlations.png")
    plt.close()

if __name__ == "__main__":
    main()