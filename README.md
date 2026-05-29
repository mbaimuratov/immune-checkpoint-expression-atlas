# Immune Checkpoint Expression Atlas

## Question
Which immune checkpoint genes are expressed across immune cell types?

## Data
Human Protein Atlas immune-cell RNA expression data from `data/raw/rna_immune_cell.tsv`.

## Genes
PDCD1, CD274, CTLA4, TIGIT, LAG3, HAVCR2, TOX, ENTPD1, ICOS, TNFRSF9.

## Methods
Data cleaning, expression filtering, heatmap visualization, top-cell-type ranking, gene correlation analysis.

## Main findings
- TIGIT has the strongest overall checkpoint signal, especially in T reg cells, memory CD8 T cells, and memory CD4 T cells.
- CTLA4 is concentrated in T reg cells.
- PDCD1 is highest in memory CD8 T cells, memory CD4 T cells, and GDT cells.
- HAVCR2 is strongest in myeloid DC, NK cells, and intermediate monocytes, giving it a distinct innate/myeloid-leaning profile.
- ENTPD1 is strongest in eosinophils, neutrophils, and classical monocytes.
- TIGIT, TNFRSF9, CTLA4, ICOS, and TOX show related T-cell-enriched expression patterns.

## How to run
Create or activate a virtual environment, then install dependencies:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Regenerate processed tables and static figures:

```bash
.venv/bin/python src/process_checkpoint_expression.py
.venv/bin/python src/plot_checkpoint_expression.py
.venv/bin/python src/compare_marker_groups.py
.venv/bin/python src/gene_gene_relationships.py
```

Run the analysis notebook:

```bash
.venv/bin/jupyter notebook notebooks/immune_checkpoint_expression_atlas.ipynb
```

Launch the interactive dashboard:

```bash
.venv/bin/python -m streamlit run streamlit_app.py
```

The dashboard lets users select a gene, inspect expression across immune cell types, select a cell type, compare checkpoint profiles, compare two genes, and view top expressing cell types.

## Limitations
Bulk/aggregated expression, not functional validation, expression does not prove protein activity.
