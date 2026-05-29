from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
CHECKPOINT_EXPRESSION_FILE = PROCESSED_DIR / "checkpoint_expression_long.csv"

GENE_GROUPS = {
    "inhibitory": ["PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "TOX"],
    "ligand": ["CD274"],
    "activation": ["ICOS", "TNFRSF9"],
    "functional": ["ENTPD1"],
}

def save_figure(filepath):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(filepath, dpi=200, bbox_inches="tight")


def compare_marker_groups(expression_data: pd.DataFrame) -> None:
    expression_data["group"] = None
    for group_name, genes in GENE_GROUPS.items():
        expression_data.loc[expression_data["gene"].isin(genes), "group"] = group_name

    group_summary = (
        expression_data.groupby(["cell_type", "group"])["expression"]
        .mean()
        .reset_index()
    )

    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=group_summary,
        x="cell_type",
        y="expression",
        hue="group",
        palette="Set2",
    )
    plt.title("Average Expression of Checkpoint Gene Groups Across Immune Cell Types")
    plt.xlabel("Immune Cell Type")
    plt.ylabel("Average TPM")
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Gene Group")
    plt.tight_layout()
    save_figure(PROJECT_ROOT / "figures" / "checkpoint_gene_groups_by_cell_type.png")
    plt.close()

if __name__ == "__main__":
    expression_data = pd.read_csv(CHECKPOINT_EXPRESSION_FILE)
    compare_marker_groups(expression_data)