import pandas as pd
import numpy as np
from logging.config import dictConfig
from api.logger import logging_config
from sklearn.model_selection import train_test_split
import logging

dictConfig(logging_config)

logger = logging.getLogger("process_norm")

def organize_dataset(
    df: pd.DataFrame, normal_count: int, tumor_count: int, sample_count: int
) -> pd.DataFrame:
    """
    add label column to the dataset
    remove duplicated columns if needed
    Parameters
    ----------
    df : pd.DataFrame
        dataset
    normal_count : int
        number of normal samples
    tumor_count : int
        number of tumor samples
    sample_count : int
        number of samples
    
    Returns
    -------
    pd.DataFrame
        organized dataset

    """
    ID = list(df.iloc[:, 0]) + ["label"]
    X = df.iloc[:, 1::sample_count].T
    if normal_count + tumor_count != X.shape[0]:
        raise ValueError(
            "normal_count + tumor_count does not match the number of samples"
        )
    y = [0] * normal_count + [1] * tumor_count
    X["label"] = y
    X.sort_values(by="label", inplace=True)
    X.index = range(X.shape[0])
    X = X.T
    X.insert(0, "Unnamed: 0", ID)
    return X


def split_dataset(
    df: pd.DataFrame, split_ratio: float, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    split the dataset into ratio_df and complement_df
    df should be organized dataset, the output of function organize_dataset

    Parameters
    ----------
    df : pd.DataFrame
        dataset
    split_ratio : float
        ratio of the dataset to split
    seed : int
        random seed

    """

    ID = list(df.iloc[:, 0])
    df_values = df.iloc[:, 1:].T
    complement_df, split_ratio_df = train_test_split(df_values, test_size=split_ratio, stratify=df_values["label"], random_state=seed)
    
    complement_df.sort_values(by="label", inplace=True)
    complement_df.index = range(len(complement_df))
    complement_df = complement_df.T
    complement_df.insert(0, "Unnamed: 0", ID)

    split_ratio_df.sort_values(by="label", inplace=True)
    split_ratio_df.index = range(len(split_ratio_df))
    split_ratio_df = split_ratio_df.T
    split_ratio_df.insert(0, "Unnamed: 0", ID)
    
    logger.info(f"complement_df feature: {complement_df.shape[0] - 1}")
    logger.info(f"complement_df sample (normal, tumor): {sum(complement_df.iloc[-1,1:].astype(int) == 0), sum(complement_df.iloc[-1,1:].astype(int) == 1)}")
    logger.info(f"ratio_df feature: {split_ratio_df.shape[0] - 1}")
    logger.info(f"ratio_df sample (normal, tumor): {sum(split_ratio_df.iloc[-1,1:].astype(int) == 0), sum(split_ratio_df.iloc[-1,1:].astype(int) == 1)}")
    return complement_df, split_ratio_df

def merge_datasets(
    df1: pd.DataFrame, 
    df2: pd.DataFrame,
) -> pd.DataFrame:
    """
    merge two datasets
    df1 and df2 should be organized datasets, the output of function organize_dataset
    
    Parameters
    ----------
    df1 : pd.DataFrame
        dataset 1
    df2 : pd.DataFrame
        dataset 2

    Returns
    -------
    pd.DataFrame
        merged dataset
    """

    df1.index = df1.iloc[:, 0].tolist()
    df2.index = df2.iloc[:, 0].tolist()
    df1 = df1.iloc[:, 1:]
    df2 = df2.iloc[:, 1:]
    concatdf = pd.concat([df1, df2], axis=1)
    concatdf = concatdf.dropna(how='any', axis=0)
    concatdf = concatdf.T
    concatdf = concatdf.sort_values(["label"])
    concatdf = concatdf.T
    concatdf.columns = range(concatdf.shape[1])
    concatdf.insert(0, "Unnamed: 0", concatdf.index)
    concatdf = concatdf.reset_index(drop=True)
    return concatdf


def inspect_nan(df: pd.DataFrame, mode: str = "column", remove: bool = False, verbose: bool = True) -> pd.DataFrame:

    """
    Find NaN values in a DataFrame.

    Parameters:
    ----------
    df (pd.DataFrame): DataFrame to search for NaNs.
    mode (str): Search mode. Choose 'column' to find NaNs in columns or 'row' to find NaNs in rows.

    Returns:
    ----------
        dict or list: Dictionary of columns with NaNs and their labels (column mode),
                      or a list of IDs with NaNs in their rows (row mode).
    """
    if mode == "column":
        labels = df.iloc[-1]  # Last row as labels
        nan_columns = df.isna().any()
        if verbose:
            report = {col: labels[col] for col in df.columns if nan_columns[col]}
            logger.info(f"Columns with NaNs: {report}")
        if remove:
            df = df.dropna(axis=1)
            df.columns = ["Unnamed: 0"] + list(range(df.shape[1] - 1))
        return df
    
    elif mode == "row":
        nan_rows = df.isna().any(axis=1)
        if verbose:
            report = df.loc[nan_rows, df.columns[0]].tolist()
            logger.info(f"Rows with NaNs: {report}")
        if remove:
            df = df.dropna(axis=0)
            df.index = range(df.shape[0])
        return df
    
    else:
        raise ValueError("Invalid mode. Choose 'column' or 'row'.")


def oversample(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    """
    oversample the dataset
    df should be organized dataset, the output of function organize_dataset

    Parameters
    ----------
    df : pd.DataFrame
        dataset
    seed : int
        random seed

    Returns
    -------
    pd.DataFrame
        oversampled dataset
    """

    # too tired to implement😿
    return pd.DataFrame()
