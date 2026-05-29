from pathlib import Path

import pandas as pd

from expression_helpers import (
    load_raw_expression,
    standardize_cell_type,
    standardize_gene_symbol,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_FILE = PROJECT_ROOT / "data" / "raw" / "rna_immune_cell.tsv"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
REPORTS_DIR = PROJECT_ROOT / "reports"

# Analysis-ready outputs.
WIDE_OUTPUT_FILE = PROCESSED_DIR / "checkpoint_expression_wide.csv"
LONG_OUTPUT_FILE = PROCESSED_DIR / "checkpoint_expression_long.csv"
RANKED_OUTPUT_FILE = REPORTS_DIR / "top_cell_types_per_gene.csv"

# Immune checkpoint genes selected for the expression atlas.
TARGET_CHECKPOINT_GENES = [
    "PDCD1",
    "CD274",
    "CTLA4",
    "TIGIT",
    "LAG3",
    "HAVCR2",
    "TOX",
    "ENTPD1",
    "ICOS",
    "TNFRSF9",
]

def clean_checkpoint_expression(raw_expression: pd.DataFrame) -> pd.DataFrame:
    required_columns = {"Gene name", "Immune cell", "TPM"}
    missing_columns = required_columns.difference(raw_expression.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Raw expression table is missing required columns: {missing}")

    expression = raw_expression.loc[:, ["Gene name", "Immune cell", "TPM"]].copy()
    expression = expression.rename(
        columns={
            "Gene name": "gene",
            "Immune cell": "cell_type",
            "TPM": "expression",
        }
    )

    expression["gene"] = expression["gene"].map(standardize_gene_symbol)
    expression["cell_type"] = expression["cell_type"].map(standardize_cell_type)
    expression["expression"] = pd.to_numeric(
        expression["expression"], errors="coerce"
    ).fillna(0.0)

    expression = expression.dropna(subset=["gene", "cell_type"])
    expression = expression[expression["gene"].isin(TARGET_CHECKPOINT_GENES)]

    found_genes = set(expression["gene"].unique())
    missing_genes = sorted(set(TARGET_CHECKPOINT_GENES).difference(found_genes))
    if missing_genes:
        print(
            "Warning: target checkpoint genes not found in raw dataset: "
            + ", ".join(missing_genes)
        )

    gene_order = {gene: index for index, gene in enumerate(TARGET_CHECKPOINT_GENES)}
    expression["gene_order"] = expression["gene"].map(gene_order)

    long_expression = (
        expression.groupby(["gene_order", "gene", "cell_type"], as_index=False)[
            "expression"
        ]
        .mean()
        .sort_values(["gene_order", "cell_type"])
        .drop(columns="gene_order")
        .reset_index(drop=True)
    )

    return long_expression


def make_wide_expression(long_expression: pd.DataFrame) -> pd.DataFrame:
    wide_expression = (
        long_expression.pivot_table(
            index="gene",
            columns="cell_type",
            values="expression",
            aggfunc="mean",
            fill_value=0.0,
        )
        .reindex(TARGET_CHECKPOINT_GENES)
        .dropna(how="all")
        .reset_index()
    )

    wide_expression.columns.name = None
    cell_type_columns = sorted(
        column for column in wide_expression.columns if column != "gene"
    )

    return wide_expression.loc[:, ["gene", *cell_type_columns]]


def make_ranked_expression(long_expression: pd.DataFrame) -> pd.DataFrame:
    ranked_expression = long_expression.sort_values(
        ["gene", "expression", "cell_type"],
        ascending=[True, False, True],
    ).copy()
    ranked_expression["rank"] = ranked_expression.groupby("gene").cumcount() + 1
    ranked_expression = ranked_expression.loc[
        :, ["gene", "rank", "cell_type", "expression"]
    ].sort_values(["gene", "rank"])
    return ranked_expression.reset_index(drop=True)


def main() -> None:
    raw_expression = load_raw_expression(RAW_FILE)
    long_expression = clean_checkpoint_expression(raw_expression)
    wide_expression = make_wide_expression(long_expression)
    ranked_expression = make_ranked_expression(long_expression)

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    wide_expression.to_csv(WIDE_OUTPUT_FILE, index=False)
    long_expression.to_csv(LONG_OUTPUT_FILE, index=False)
    ranked_expression.to_csv(RANKED_OUTPUT_FILE, index=False)

    print(f"Saved wide checkpoint expression table to {WIDE_OUTPUT_FILE}")
    print(f"Saved long checkpoint expression table to {LONG_OUTPUT_FILE}")
    print(f"Saved ranked checkpoint expression table to {RANKED_OUTPUT_FILE}")


if __name__ == "__main__":
    main()
