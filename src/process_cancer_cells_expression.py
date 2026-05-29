import argparse
from collections.abc import Iterator
from pathlib import Path

import pandas as pd

from expression_helpers import standardize_gene_symbol

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_FILE = PROJECT_ROOT / "data" / "raw" / "rna_cancer_sample.tsv.gz"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
REPORTS_DIR = PROJECT_ROOT / "reports"
DEFAULT_CHUNKSIZE = 500_000
RAW_COLUMNS = ["Gene", "Sample", "Cancer", "pTPM"]

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

# Human Ensembl IDs used by the HPA cancer sample RNA table.
CANCER_GENE_ID_TO_SYMBOL = {
    "ENSG00000188389": "PDCD1",
    "ENSG00000120217": "CD274",
    "ENSG00000163599": "CTLA4",
    "ENSG00000181847": "TIGIT",
    "ENSG00000089692": "LAG3",
    "ENSG00000135077": "HAVCR2",
    "ENSG00000198846": "TOX",
    "ENSG00000138185": "ENTPD1",
    "ENSG00000163600": "ICOS",
    "ENSG00000049249": "TNFRSF9",
}
LONG_OUTPUT_FILE = PROCESSED_DIR / "cancer_checkpoint_expression_long.csv"
WIDE_OUTPUT_FILE = PROCESSED_DIR / "cancer_checkpoint_expression_wide.csv"
RANKED_OUTPUT_FILE = REPORTS_DIR / "cancer_checkpoint_expression_ranked.csv"


def _validate_raw_columns(columns: pd.Index) -> None:
    missing_columns = set(RAW_COLUMNS).difference(columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Raw expression table is missing required columns: {missing}")


def _standardize_gene_identifiers(genes: pd.Series) -> pd.Series:
    gene_keys = (
        genes.astype("string")
        .str.strip()
        .str.upper()
        .str.replace(r"\..*$", "", regex=True)
    )
    gene_symbols = gene_keys.map(CANCER_GENE_ID_TO_SYMBOL)

    non_ensembl_mask = gene_symbols.isna() & ~gene_keys.str.startswith(
        "ENSG", na=False
    )
    if non_ensembl_mask.any():
        symbol_candidates = gene_keys.loc[non_ensembl_mask].map(
            standardize_gene_symbol
        )
        symbol_candidates = symbol_candidates.where(
            symbol_candidates.isin(TARGET_CHECKPOINT_GENES)
        )
        gene_symbols.loc[non_ensembl_mask] = symbol_candidates

    return gene_symbols


def _prepare_cancer_checkpoint_expression(
    raw_expression: pd.DataFrame,
) -> pd.DataFrame:
    _validate_raw_columns(raw_expression.columns)

    expression = raw_expression.loc[:, ["Gene", "Sample", "Cancer", "pTPM"]].copy()
    expression["gene"] = _standardize_gene_identifiers(expression["Gene"])
    expression = expression[expression["gene"].isin(TARGET_CHECKPOINT_GENES)]
    if expression.empty:
        return pd.DataFrame(columns=["gene", "sample", "cancer_type", "expression"])

    expression = expression.rename(
        columns={
            "Sample": "sample",
            "Cancer": "cancer_type",
            "pTPM": "expression",
        }
    )
    expression["sample"] = expression["sample"].astype("string").str.strip()
    expression["cancer_type"] = expression["cancer_type"].str.strip()
    expression["expression"] = pd.to_numeric(
        expression["expression"], errors="coerce"
    ).fillna(0.0)

    expression = expression.dropna(subset=["gene", "sample", "cancer_type"])
    return expression.loc[:, ["gene", "sample", "cancer_type", "expression"]]


def _sort_by_target_gene(
    expression: pd.DataFrame,
    sort_columns: list[str],
) -> pd.DataFrame:
    gene_order = {gene: index for index, gene in enumerate(TARGET_CHECKPOINT_GENES)}
    return (
        expression.assign(_gene_order=expression["gene"].map(gene_order))
        .sort_values(["_gene_order", *sort_columns])
        .drop(columns="_gene_order")
        .reset_index(drop=True)
    )


def _finalize_long_cancer_expression(expression: pd.DataFrame) -> pd.DataFrame:
    if expression.empty:
        return pd.DataFrame(columns=["gene", "sample", "cancer_type", "expression"])

    long_expression = (
        expression.groupby(["gene", "sample", "cancer_type"], as_index=False)[
            "expression"
        ]
        .mean()
    )

    return _sort_by_target_gene(long_expression, ["cancer_type", "sample"])


def clean_cancer_checkpoint_expression(raw_expression: pd.DataFrame) -> pd.DataFrame:
    expression = _prepare_cancer_checkpoint_expression(raw_expression)
    return _finalize_long_cancer_expression(expression)


def iter_cancer_checkpoint_expression_chunks(
    raw_file: Path = RAW_FILE,
    chunksize: int = DEFAULT_CHUNKSIZE,
) -> Iterator[pd.DataFrame]:
    header = pd.read_csv(raw_file, sep="\t", compression="infer", nrows=0)
    _validate_raw_columns(header.columns)

    raw_chunks = pd.read_csv(
        raw_file,
        sep="\t",
        compression="infer",
        usecols=RAW_COLUMNS,
        dtype={"Gene": "string", "Sample": "string", "Cancer": "string"},
        chunksize=chunksize,
    )
    for raw_chunk in raw_chunks:
        cleaned_chunk = _prepare_cancer_checkpoint_expression(raw_chunk)
        if not cleaned_chunk.empty:
            yield cleaned_chunk


def load_cancer_checkpoint_expression(
    raw_file: Path = RAW_FILE,
    chunksize: int = DEFAULT_CHUNKSIZE,
    save: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    chunks = list(iter_cancer_checkpoint_expression_chunks(raw_file, chunksize))
    if chunks:
        filtered_expression = pd.concat(chunks, ignore_index=True)
    else:
        filtered_expression = pd.DataFrame(
            columns=["gene", "sample", "cancer_type", "expression"]
        )

    long_expression = _finalize_long_cancer_expression(filtered_expression)
    _warn_missing_target_genes(long_expression)
    wide_expression = make_wide_cancer_expression(long_expression)
    ranked_expression = make_ranked_cancer_expression(long_expression)

    if save:
        save_cancer_checkpoint_expression(
            long_expression,
            wide_expression,
            ranked_expression,
        )

    return long_expression, wide_expression, ranked_expression


def make_wide_cancer_expression(long_expression: pd.DataFrame) -> pd.DataFrame:
    cancer_summary = (
        long_expression.groupby(["gene", "cancer_type"], as_index=False)["expression"]
        .mean()
    )
    wide_expression = (
        cancer_summary.pivot_table(
            index="gene",
            columns="cancer_type",
            values="expression",
            aggfunc="mean",
            fill_value=0.0,
        )
        .reindex(TARGET_CHECKPOINT_GENES)
        .dropna(how="all")
        .reset_index()
    )

    wide_expression.columns.name = None
    cancer_type_columns = sorted(
        column for column in wide_expression.columns if column != "gene"
    )

    return wide_expression.loc[:, ["gene", *cancer_type_columns]]


def make_ranked_cancer_expression(long_expression: pd.DataFrame) -> pd.DataFrame:
    ranked_expression = (
        long_expression.groupby(["gene", "cancer_type"], as_index=False)["expression"]
        .mean()
        .sort_values(
            ["gene", "expression", "cancer_type"],
            ascending=[True, False, True],
        )
    )
    ranked_expression["rank"] = ranked_expression.groupby("gene").cumcount() + 1
    ranked_expression = ranked_expression.loc[
        :, ["gene", "rank", "cancer_type", "expression"]
    ]
    return _sort_by_target_gene(ranked_expression, ["rank"])


def _warn_missing_target_genes(long_expression: pd.DataFrame) -> None:
    found_genes = set(long_expression["gene"].unique())
    missing_genes = sorted(set(TARGET_CHECKPOINT_GENES).difference(found_genes))
    if missing_genes:
        print(
            "Warning: target checkpoint genes not found in raw dataset: "
            + ", ".join(missing_genes)
        )


def save_cancer_checkpoint_expression(
    long_expression: pd.DataFrame,
    wide_expression: pd.DataFrame,
    ranked_expression: pd.DataFrame,
) -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    long_expression.to_csv(LONG_OUTPUT_FILE, index=False)
    wide_expression.to_csv(WIDE_OUTPUT_FILE, index=False)
    ranked_expression.to_csv(RANKED_OUTPUT_FILE, index=False)

    print(f"Saved long cancer checkpoint expression table to {LONG_OUTPUT_FILE}")
    print(f"Saved wide cancer checkpoint expression table to {WIDE_OUTPUT_FILE}")
    print(f"Saved ranked cancer checkpoint expression table to {RANKED_OUTPUT_FILE}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stream the compressed HPA cancer sample RNA table."
    )
    parser.add_argument(
        "--raw-file",
        type=Path,
        default=RAW_FILE,
        help="Path to rna_cancer_sample.tsv.gz.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=DEFAULT_CHUNKSIZE,
        help="Rows to read per pandas chunk.",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Write long, wide, and ranked CSV outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    long_expression, wide_expression, ranked_expression = (
        load_cancer_checkpoint_expression(
            raw_file=args.raw_file,
            chunksize=args.chunksize,
            save=args.save,
        )
    )

    print(
        "Loaded cancer checkpoint expression: "
        f"{len(long_expression):,} long rows, "
        f"{len(wide_expression):,} genes, "
        f"{long_expression['cancer_type'].nunique():,} cancer types."
    )
    if not args.save:
        print("CSV outputs were not saved. Re-run with --save to write them.")


if __name__ == "__main__":
    main()
