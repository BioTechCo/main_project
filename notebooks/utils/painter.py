from plotly import graph_objects as go
from typing import Literal
from pandas import DataFrame


def plot_roc_curve(
    df: DataFrame,
    title: str,
    outPath: str,
    width: int = 1200,
    height: int = 900,
    mode: Literal["train", "test"] = "train",
):
    fig = go.Figure()
    for index, row in df.iterrows():
        fig.add_trace(
            go.Scatter(
                x=row[f"fpr_{mode}"],
                y=row[f"tpr_{mode}"],
                mode="lines",
                name=f"{row['selection_model']} {row['train_model']} (AUC = {row[f'AUC_{mode}']:.3f})",
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


def plot_performance_diff(
    df: DataFrame, colors: dict, outPath: str, width: int = 800, height: int = 1600
):
    fig = go.Figure()
    for key, value in colors.items():
        fig.add_trace(
            go.Scatter(
                x=df[f"{key}_diff"],
                y=df["selection_train_model_features"],
                mode="lines+markers",
                name=f"{key} Difference",
                line=dict(color=value, width=2),
                marker=dict(color=value, size=8),
            )
        )

    # Update layout
    fig.update_layout(
        title="Performance Difference between Training and Testing Set",
        xaxis_title="Performance Difference (Training - Testing)",
        yaxis_title="Combination of Selection Model, Train Model and Features",
        template="plotly_white",
        width=width,
        height=height,
    )

    fig.write_html(outPath)
    print(f"Performance difference saved to {outPath}")
