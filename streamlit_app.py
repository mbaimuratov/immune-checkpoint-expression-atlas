from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


PROJECT_ROOT = Path(__file__).resolve().parent
LONG_EXPRESSION_FILE = PROJECT_ROOT / "data" / "processed" / "checkpoint_expression_long.csv"
WIDE_EXPRESSION_FILE = PROJECT_ROOT / "data" / "processed" / "checkpoint_expression_wide.csv"
RANKED_EXPRESSION_FILE = PROJECT_ROOT / "reports" / "top_cell_types_per_gene.csv"


ACRONYMS = {
    "b": "B",
    "cd4": "CD4",
    "cd8": "CD8",
    "dc": "DC",
    "gdt": "GDT",
    "mait": "MAIT",
    "nk": "NK",
    "pbmc": "PBMC",
    "t": "T",
}


def format_cell_type(cell_type: str) -> str:
    parts = str(cell_type).split("_")
    labels = [ACRONYMS.get(part, part.capitalize()) for part in parts]
    return " ".join(labels).replace(" Pos ", "+ ")


@st.cache_data
def load_expression_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    long_expression = pd.read_csv(LONG_EXPRESSION_FILE)
    wide_expression = pd.read_csv(WIDE_EXPRESSION_FILE)
    ranked_expression = pd.read_csv(RANKED_EXPRESSION_FILE)

    required_long = {"gene", "cell_type", "expression"}
    required_wide = {"gene"}
    required_ranked = {"gene", "rank", "cell_type", "expression"}

    missing_long = required_long.difference(long_expression.columns)
    missing_wide = required_wide.difference(wide_expression.columns)
    missing_ranked = required_ranked.difference(ranked_expression.columns)
    if missing_long or missing_wide or missing_ranked:
        raise ValueError(
            "Input CSV schemas do not match the dashboard contract. "
            f"Missing long={sorted(missing_long)}, "
            f"wide={sorted(missing_wide)}, ranked={sorted(missing_ranked)}."
        )

    long_expression["expression"] = pd.to_numeric(
        long_expression["expression"], errors="coerce"
    ).fillna(0.0)
    ranked_expression["expression"] = pd.to_numeric(
        ranked_expression["expression"], errors="coerce"
    ).fillna(0.0)

    return long_expression, wide_expression, ranked_expression


def add_cell_labels(expression: pd.DataFrame) -> pd.DataFrame:
    labeled = expression.copy()
    labeled["cell_label"] = labeled["cell_type"].map(format_cell_type)
    return labeled


def expression_bar(
    expression: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    color: str = "#277c78",
) -> go.Figure:
    figure = px.bar(expression, x=x, y=y, title=title, text_auto=".2g")
    figure.update_traces(marker_color=color, textposition="outside", cliponaxis=False)
    figure.update_layout(
        height=430,
        margin=dict(l=10, r=10, t=60, b=40),
        xaxis_title=None,
        yaxis_title="Average TPM",
        title_font_size=18,
        showlegend=False,
    )
    return figure


st.set_page_config(
    page_title="Immune Checkpoint Expression Atlas",
    page_icon="IC",
    layout="wide",
)

st.markdown(
    """
    <style>
    .block-container {padding-top: 2rem; padding-bottom: 3rem;}
    div[data-testid="stMetric"] {
        border: 1px solid #d7dee8;
        border-radius: 8px;
        padding: 0.85rem 1rem;
        background: #f8fafc;
    }
    div[data-testid="stMetricValue"] {font-size: 1.55rem;}
    </style>
    """,
    unsafe_allow_html=True,
)


try:
    long_expression, wide_expression, ranked_expression = load_expression_data()
except Exception as error:
    st.error(f"Dashboard data could not be loaded: {error}")
    st.stop()


genes = wide_expression["gene"].tolist()
cell_types = [column for column in wide_expression.columns if column != "gene"]

st.title("Immune Checkpoint Expression Atlas")
st.caption(
    "Interactive view of Human Protein Atlas immune-cell RNA expression for selected "
    "checkpoint and checkpoint-associated genes. Values are average TPM."
)

with st.sidebar:
    st.header("Controls")
    selected_gene = st.selectbox("Gene", genes, index=0)
    selected_cell_type = st.selectbox(
        "Cell type",
        cell_types,
        index=cell_types.index("t_reg") if "t_reg" in cell_types else 0,
        format_func=format_cell_type,
    )
    first_gene = st.selectbox("Compare gene A", genes, index=0)
    second_gene = st.selectbox(
        "Compare gene B",
        genes,
        index=genes.index("TIGIT") if "TIGIT" in genes else min(1, len(genes) - 1),
    )
    top_n = st.slider("Top expressing cell types", 3, min(15, len(cell_types)), 8)

metric_columns = st.columns(4)
metric_columns[0].metric("Genes", len(genes))
metric_columns[1].metric("Cell types", len(cell_types))
metric_columns[2].metric(
    "Max TPM", f"{long_expression['expression'].max():.1f}"
)
metric_columns[3].metric(
    "Zero values", f"{(long_expression['expression'] == 0).sum()}"
)

gene_expression = add_cell_labels(
    long_expression.loc[long_expression["gene"] == selected_gene]
).sort_values("expression", ascending=False)

cell_profile = add_cell_labels(
    long_expression.loc[long_expression["cell_type"] == selected_cell_type]
).sort_values("expression", ascending=False)

top_ranked = add_cell_labels(
    ranked_expression.loc[ranked_expression["gene"] == selected_gene]
    .nsmallest(top_n, "rank")
    .sort_values("expression", ascending=True)
)

st.subheader(f"{selected_gene} expression across immune cell types")
st.plotly_chart(
    expression_bar(
        gene_expression,
        x="cell_label",
        y="expression",
        title=f"{selected_gene} across immune cell types",
    ),
    width="stretch",
)

profile_column, top_column = st.columns(2)

with profile_column:
    st.subheader(f"Checkpoint profile in {format_cell_type(selected_cell_type)}")
    st.plotly_chart(
        expression_bar(
            cell_profile,
            x="gene",
            y="expression",
            title=f"{format_cell_type(selected_cell_type)} checkpoint profile",
            color="#4b6cb7",
        ),
        width="stretch",
    )

with top_column:
    st.subheader(f"Top {top_n} cell types for {selected_gene}")
    figure = px.bar(
        top_ranked,
        x="expression",
        y="cell_label",
        orientation="h",
        title=f"Highest {selected_gene} expression",
        text_auto=".2g",
    )
    figure.update_traces(marker_color="#c47f2c", textposition="outside", cliponaxis=False)
    figure.update_layout(
        height=430,
        margin=dict(l=10, r=10, t=60, b=40),
        xaxis_title="Average TPM",
        yaxis_title=None,
        title_font_size=18,
        showlegend=False,
    )
    st.plotly_chart(figure, width="stretch")

st.subheader("Compare two genes")
compare_data = add_cell_labels(
    long_expression.loc[long_expression["gene"].isin([first_gene, second_gene])]
)
compare_order = (
    compare_data.groupby("cell_label")["expression"]
    .max()
    .sort_values(ascending=False)
    .index.tolist()
)
compare_data["cell_label"] = pd.Categorical(
    compare_data["cell_label"], categories=compare_order, ordered=True
)
compare_data = compare_data.sort_values("cell_label")

correlation_matrix = wide_expression.set_index("gene").T.corr()
correlation = correlation_matrix.loc[first_gene, second_gene]

compare_metric, compare_chart = st.columns([1, 4])
with compare_metric:
    st.metric("Correlation", f"{correlation:.2f}")
    if first_gene == second_gene:
        st.info("Select two different genes to compare distinct expression patterns.")
with compare_chart:
    compare_figure = px.bar(
        compare_data,
        x="cell_label",
        y="expression",
        color="gene",
        barmode="group",
        title=f"{first_gene} vs {second_gene}",
        text_auto=".2g",
        color_discrete_sequence=["#277c78", "#4b6cb7"],
    )
    compare_figure.update_traces(textposition="outside", cliponaxis=False)
    compare_figure.update_layout(
        height=430,
        margin=dict(l=10, r=10, t=60, b=40),
        xaxis_title=None,
        yaxis_title="Average TPM",
        title_font_size=18,
        legend_title_text="Gene",
    )
    st.plotly_chart(compare_figure, width="stretch")

with st.expander("Full checkpoint expression heatmap and data table"):
    heatmap_matrix = wide_expression.set_index("gene")
    heatmap_figure = go.Figure(
        data=go.Heatmap(
            z=heatmap_matrix.to_numpy().T,
            x=heatmap_matrix.index.tolist(),
            y=[format_cell_type(cell_type) for cell_type in heatmap_matrix.columns],
            colorscale="Viridis",
            colorbar=dict(title="TPM"),
        )
    )
    heatmap_figure.update_layout(
        height=620,
        margin=dict(l=10, r=10, t=25, b=40),
        xaxis_title="Gene",
        yaxis_title="Immune cell type",
    )
    st.plotly_chart(heatmap_figure, width="stretch")

    table = add_cell_labels(long_expression).loc[
        :, ["gene", "cell_label", "expression"]
    ]
    table = table.rename(
        columns={
            "gene": "Gene",
            "cell_label": "Cell type",
            "expression": "Average TPM",
        }
    )
    st.dataframe(table, width="stretch", hide_index=True)
