import pandas as pd
from .config_helper import update_nested_toml
from api.logger import logging_config
from logging.config import dictConfig
import logging


import seaborn as sns
import matplotlib.pyplot as plt


dictConfig(logging_config)

logger = logging.getLogger("dbeta_avg_helper")

def IQR(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate the Interquartile Range (IQR) for the given DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        tuple: Upper and lower fences for outlier detection.
    """
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    upper_fence = Q3 + IQR * 1.5
    lower_fence = Q1 - IQR * 1.5
    return upper_fence, lower_fence


def no_outlier(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove outliers from the given DataFrame.
    
    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: DataFrame without outliers.
    """
    upper_fence, lower_fence = IQR(df)
    ddf = df[(df > lower_fence) & (df < upper_fence)]
    return ddf

def get_dbeta_avg(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the mean of the tumor samples.
    
    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: Mean of the tumor samples.
    """
    normal_count = int((df.iloc[-1, 1:] == 0).sum())
    all_beta_normalized_normal = df.iloc[:-1, 1 : normal_count + 1 :].T

    all_beta_normalized_tumor = df.iloc[:-1, normal_count + 1 : :].T

    all_beta_normalized_normal = no_outlier(all_beta_normalized_normal)
    all_beta_normalized_tumor = no_outlier(all_beta_normalized_tumor)

    train_normal_avg = all_beta_normalized_normal.mean(skipna=True, axis=0)

    all_beta_normalized_tumor = all_beta_normalized_tumor.subtract(
        train_normal_avg, axis=1
    )

    all_beta_normalized_tumor = no_outlier(all_beta_normalized_tumor)

    train_tumor_mean = all_beta_normalized_tumor.mean(skipna=True, axis=0)

    delta_beta = pd.merge(
        df.iloc[:-1, :1],
        pd.DataFrame(train_tumor_mean, columns=["dbeta"]),
        left_index=True,
        right_index=True,
    )

    return delta_beta

def drop_dbeta_nan(delta_beta: pd.DataFrame, log_postfix:str="") -> pd.DataFrame:
    """
    Drop NaN values from the given DataFrame.
    
    Args:
        delta_beta (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: DataFrame without NaN values.
    """
    update_nested_toml("preprocess.dbeta", f"delta_beta_avg_feature_num{log_postfix}", delta_beta.shape[0])
    update_nested_toml("preprocess.dbeta", f"NaN_dbeta_feature{log_postfix}", delta_beta.loc[pd.isna(delta_beta["dbeta"]), "Unnamed: 0"].tolist(),)
    delta_beta = delta_beta.dropna(axis=0)
    update_nested_toml("preprocess.dbeta", f"delta_beta_avg_feature_num_remove_NaN{log_postfix}", delta_beta.shape[0])
    return delta_beta

def find_max_dBeta_grouped(group):
    idx_max = group["dbeta"].abs().idxmax()
    return group.loc[idx_max]

def get_dbeta_info(delta_beta: pd.DataFrame, dmp: pd.DataFrame, log_postfix:str="") -> pd.DataFrame:
    """
    Get the dmp information.
    
    Args:
        delta_beta (pd.DataFrame): Input DataFrame.
        dmp (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: dmp information.
    """
    dmp = dmp[["Unnamed: 0", "gene", "feature"]]
    update_nested_toml("preprocess.dbeta", f"dmp_before_dropna_shape_feature{log_postfix}", dmp.shape[0])

    dmp = dmp.dropna(axis=0)

    update_nested_toml("preprocess.dbeta", f"dmp_after_dropna_shape_feature{log_postfix}", dmp.shape[0])

    result = pd.merge(delta_beta, dmp, on="Unnamed: 0", how="left")
    update_nested_toml(
        "preprocess.dbeta", f"delta_beta_avg_feature_num_remove_NaN_join_dmp{log_postfix}", result.shape[0]
    )
    dbeta_info = result.groupby("gene", as_index=False).apply(
        find_max_dBeta_grouped, include_groups=False
    )

    dbeta_info.columns = ["gene", "ID", "dbeta", "feature"]
    dbeta_info = dbeta_info[["ID", "gene", "dbeta", "feature"]]
    
    return dbeta_info

def detect_threshold(
    dbeta_info: pd.DataFrame,
    config: dict,
    log_postfix:str=""
) -> tuple[pd.DataFrame, float]:
    """
    Detect the threshold of the given DataFrame.
    
    Args:
        dbeta_info (pd.DataFrame): Input DataFrame.
        config (dict): Configuration settings.

    Returns:
        pd.DataFrame: Threshold of the given DataFrame.
    """
    threshold = 1
    dbeta_TSS_threshold = dbeta_info[abs(dbeta_info["dbeta"]) > threshold]
    while True:
        dbeta_TSS_threshold = dbeta_info[abs(dbeta_info["dbeta"]) > threshold]
        count = dbeta_TSS_threshold.shape[0]
        if (
            config["preprocess"]["filtering"]["hyper"]["avg_dbeta_lower_bound"]
            <= count
            <= config["preprocess"]["filtering"]["hyper"]["avg_dbeta_upper_bound"]
        ):
            break
        threshold -= 0.01
    threshold = round(threshold, 2)
    update_nested_toml("preprocess.filtering", f"threshold{log_postfix}", threshold)
    logger.info(f"Threshold: {threshold}")
    return dbeta_TSS_threshold, threshold

def dbeta_graph(dbeta_TSS_threshold: pd.DataFrame, out_path: str):
    sns.kdeplot(dbeta_TSS_threshold["dbeta"])
    plt.xlabel("delta Beta value")
    plt.title("Density plot of delta Beta value")
    plt.savefig(out_path)
    plt.close()


import plotly.express as px
import pandas as pd
from sklearn.decomposition import PCA

def pca_graph(dbeta_info: pd.DataFrame, df: pd.DataFrame, out_path: str):
    
    beta_df = df.iloc[:-1, :]
    beta_df = beta_df[beta_df["Unnamed: 0"].isin(dbeta_info["ID"])]
    X = beta_df.iloc[:, 1:].dropna(axis=0).T
    y = df.iloc[-1, 1:].astype(int).to_list()

    logger.info(f"X shape: {X.shape}")
    logger.info(f"y shape: {len(y)}")

        
    pca = PCA(n_components=3)
    X_pca = pca.fit_transform(X)

    df = pd.DataFrame(
        {
            "Principal Component 1": X_pca[:, 0],
            "Principal Component 2": X_pca[:, 1],
            "Principal Component 3": X_pca[:, 2],
            "Class": y,
        }
    )
    fig = px.scatter_3d(
        df,
        x="Principal Component 1",
        y="Principal Component 2",
        z="Principal Component 3",
        color="Class",
        title="PCA of Dataset",
        color_continuous_scale="Viridis",
    )

    fig.update_layout(
        scene=dict(
            xaxis_title="Principal Component 1",
            yaxis_title="Principal Component 2",
            zaxis_title="Principal Component 3",
        )
    )

    fig.write_html(out_path)