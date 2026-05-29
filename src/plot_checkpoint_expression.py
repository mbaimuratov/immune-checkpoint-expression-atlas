from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
CHECKPOINT_EXPRESSION_FILE = PROCESSED_DIR / "checkpoint_expression_wide.csv"
RANKED_EXPRESSION_FILE = PROJECT_ROOT / "reports" / "top_cell_types_per_gene.csv"
FIGURE_DIR = PROJECT_ROOT / "figures"
TOP_GENE_PLOTS = ["TIGIT", "PDCD1", "CTLA4", "LAG3", "HAVCR2"]
TOP_CELL_TYPES_TO_PLOT = 8

def save_figure(filepath):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(filepath, dpi=200, bbox_inches="tight")


def plot_checkpoint_expression(checkpoint_expression, title=None):
    plt.figure(figsize=(10, 6))
    sns.heatmap(
        checkpoint_expression.set_index("gene").T,
        cmap="viridis",
        cbar_kws={"label": "Average TPM"},
    )
    plt.title(title or "Checkpoint Gene Expression Across Immune Cell Types")
    plt.xlabel("Checkpoint Genes")
    plt.ylabel("Immune Cell Types")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    save_figure(FIGURE_DIR / "checkpoint_expression_heatmap.png")
    plt.close()


def plot_top_cell_types_per_gene(ranked_expression: pd.DataFrame) -> None:
    plot_dir = FIGURE_DIR / "top_cell_types"
    for gene in TOP_GENE_PLOTS:
        gene_data = ranked_expression.loc[ranked_expression["gene"] == gene]
        if gene_data.empty:
            continue

        subset = gene_data.nsmallest(TOP_CELL_TYPES_TO_PLOT, "rank").sort_values(
            "expression"
        )

        plt.figure(figsize=(8, 4.5))
        sns.barplot(
            data=subset,
            x="expression",
            y="cell_type",
            color="#4c72b0",
            order=subset["cell_type"].tolist(),
        )
        plt.title(f"Top immune cell types for {gene}")
        plt.xlabel("Average TPM")
        plt.ylabel("Immune cell type")
        plt.tight_layout()
        save_figure(plot_dir / f"{gene.lower()}_top_cell_types.png")
        plt.close()


def main():
    checkpoint_expression = pd.read_csv(CHECKPOINT_EXPRESSION_FILE)
    ranked_expression = pd.read_csv(RANKED_EXPRESSION_FILE)
    plot_checkpoint_expression(checkpoint_expression)
    plot_top_cell_types_per_gene(ranked_expression)


if __name__ == "__main__":
    main()
