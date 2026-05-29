from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_FILE = PROJECT_ROOT / "data" / "raw" / "rna_immune_cell.tsv"

# Common aliases for checkpoint genes. Keys are normalized before lookup.
GENE_SYMBOL_ALIASES = {
    "PD-1": "PDCD1",
    "PD1": "PDCD1",
    "PDL1": "CD274",
    "PD-L1": "CD274",
    "CTLA-4": "CTLA4",
    "LAG-3": "LAG3",
    "TIM-3": "HAVCR2",
    "CD39": "ENTPD1",
    "4-1BB": "TNFRSF9",
}


def standardize_gene_symbol(symbol: object) -> str | None:
    """Return an uppercase HGNC-like gene symbol with known aliases resolved."""
    if pd.isna(symbol):
        return None

    standardized = str(symbol).strip().upper()
    if not standardized:
        return None

    return GENE_SYMBOL_ALIASES.get(standardized, standardized)


def standardize_cell_type(cell_type: object) -> str | None:
    """Return a stable, column-friendly immune cell type label."""
    if pd.isna(cell_type):
        return None

    standardized = (
        str(cell_type)
        .strip()
        .lower()
        .replace("+", "pos")
        .replace("-", "_")
        .replace(" ", "_")
        .replace("/", "_")
    )
    if not standardized:
        return None

    return standardized


def load_raw_expression(path: Path = RAW_FILE) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")
