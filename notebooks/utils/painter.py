from plotly import graph_objects as go
from typing import Literal
from pandas import DataFrame


def plot_roc_curve(
    df: DataFrame,
    title: str,
    outPath: str,
    width: int = 1200,
    height: int = 900,
    x_column: str = "fpr",
    y_column: str = "tpr",
    trace_name: list[str] = ["selection_model", "train_model"],
):
    """
    Plot the ROC curve for the given data frame
    
    Parameters
    ----------
    df : DataFrame
        The data frame containing the ROC curve data
    title : str
        The title of the plot
    outPath : str
        The path to save the plot
    width : int
        The width of the plot
    height : int
        The height of the plot
    x_column : str
        The name of the column containing the x values
    y_column : str
        The name of the column containing the y values
    trace_name : list[str]
        The list of column names to use as the trace names
    """
    fig = go.Figure()
    for _, row in df.iterrows():
        fig.add_trace(
            go.Scatter(
                x=row[x_column],
                y=row[y_column],
                mode="lines",
                name="_".join([str(x) for x in row[trace_name]]),
            )
        )

    # Add the random guess line
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(dash="dash", color="grey"),
            name="Random Guess",
        )
    )

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title="False Positive Rate (FPR)",
        yaxis_title="True Positive Rate (TPR)",
        legend_title="Model",
        template="plotly_white",
        width=width,
        height=height,
    )

    # store the plot
    fig.write_html(outPath)
    print(f"ROC curve saved to {outPath}")


# plot the accuracy difference using ploty


def create_performance_barchart(
    df: DataFrame,
    color_mapping: dict,
    metric: str,
    out_path: str,
    title: str,
    x_axis_label: str,
    y_axis_label: str,
    width: int = 800,
    height: int = 1600,
    orientation: Literal["v", "h"] = "v",
    template: str = "plotly_white",
):
    fig = go.Figure()
    for key, value in color_mapping.items():
        fig.add_trace(
            go.Bar(
                x=df[f"{key}"],
                y=df[metric],
                name=f"{key}",
                marker=dict(color=value),
                orientation=orientation,
            )
        )

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title=x_axis_label,
        yaxis_title=y_axis_label,
        template=template,
        width=width,
        height=height,
    )

    fig.write_html(out_path)
    print(f"Performance difference saved to {out_path}")

