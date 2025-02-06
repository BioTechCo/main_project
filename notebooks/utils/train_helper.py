import pandas as pd
from sklearn.feature_selection import RFE, SequentialFeatureSelector
import random
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import (
    accuracy_score,
    recall_score,
    precision_score,
    f1_score,
    matthews_corrcoef,
    auc,
    roc_curve,
    confusion_matrix,
    fbeta_score,
    make_scorer,
)
import os

from .pipe import pipe_dec
from collections import defaultdict
from logging.config import dictConfig
from api.logger import logging_config
import logging
from .pipe import pipe_dec
import json
from typing import Literal, Union

dictConfig(logging_config)

logger = logging.getLogger("train_helper")


def set_parameters(model: object, param_grid: dict) -> object:
    grid_search = GridSearchCV(
        estimator=model, param_grid=param_grid, cv=5, n_jobs=-1, verbose=2
    )
    return grid_search


class GeneRank:
    def __init__(self, gene, rank):
        self.gene: str = gene
        self.rank: int = rank

    def __str__(self) -> str:
        return f"<Gene: {self.gene}, Rank: {self.rank}>"

    def __repr__(self) -> str:
        return f"<Gene: {self.gene}, Rank: {self.rank}>"

    def getGene(self) -> str:
        return self.gene

    def getRank(self) -> int:
        return self.rank


def read_selected_features(file_name: str) -> defaultdict:
    """
    Read the selected features from a txt file.
    The selected features are saved in a txt file for another selection process.
    The result is a dictionary with the gene name as the key and the selected features as the value.

    Parameters:
        file_name (str): The path to the txt file containing the selected features.

    Returns:
        dict: The dictionary of selected features for each gene.
    """
    gene_dict = defaultdict(list)
    with open(file_name, "r") as f:
        for line in f:
            selected_feature_names = line.split(",")[1:]
            selected_feature_names[-1] = selected_feature_names[-1].strip()
            gene_dict[line.split(",")[0]].append(selected_feature_names)
    return gene_dict


def read_selected_features_json(file_name: str) -> dict:
    """
    Read the selected features from a json file.
    The selected features are saved in a json file for simple model training.
    The result is a dictionary with the model name as the key and the selected features as the value, thereby removing the cluster information.

    Parameters:
        file_name (str): The path to the json file containing the selected features.

    Returns:
        dict: The dictionary of selected features for each model.
    """
    with open(file_name, "r") as f:
        gene_dict: dict[str, dict[str, list[str]]] = json.load(f)
    result_dict: dict[str, list[str]] = defaultdict(list)
    for model, cluster_dict in gene_dict.items():
        for _, genes in cluster_dict.items():
            result_dict[model] += genes
    return result_dict


class TrainHelper:
    def __init__(self, dbeta_info: pd.DataFrame) -> None:
        self.dbeta_info = dbeta_info
        self.selection_model = {}
        self.grid_estimators = {}  # only used by RFE
        self.Bagging_num = 5  # only used by RFE

    def setup_dbeta(self, gene_list: list[str]) -> None:
        """
        Setup the dbeta_info dataframe with the gene_list if there is a gene list selected by previous dataset

        Parameters:
            gene_list (list[str]): The list of genes to filter the dbeta_info dataframe.
        """
        self.dbeta_info = self.dbeta_info[self.dbeta_info["gene"].isin(gene_list)]

    def set_train_df(self, train_df: pd.DataFrame) -> None:
        self.train_df = train_df

    def set_train_validate_df(
        self, train_df: pd.DataFrame, validate_df: pd.DataFrame
    ) -> None:
        # used by rfe
        self.train_df = train_df
        self.validate_df = validate_df

    def generate_selected_features(
        self,
        gene_dict: dict[str, list[list[str]]],
        out_path: str,
        mode: Union[Literal["min", "max"], int] = "min",
        out_format: Literal["json", "txt"] = "json",
    ) -> None:
        """
        Generate the selected features of clusters based on the gene_dict.
        The selected features are saved in the out_path file.

        Parameters:
            gene_dict (dict[str, list[list[str]]): The dictionary of gene lists for each model.
            out_path (str): The path to save the selected features.
            mode (Literal["min", "max"]): The mode to select the features.
            If "min", the selected features are the genes with the lowest rank in each cluster.
            If "max", the selected features are the genes with the highest rank in each cluster.
            out_format (Literal["json", "txt"]): The output format of the selected features.
            If "json", the selected features are saved in a json file for simple model training.
            If "csv", the selected features are saved in a txt file for another selection process.
        """
        result_dict: dict[str, dict[str, list[str]]] = defaultdict(dict)
        result_set: set[str] = set()

        if mode == "min":
            category_dict: dict[str, dict[str, list[GeneRank]]] = defaultdict(dict)
            for model, gene_lists in gene_dict.items():
                record = []
                data = defaultdict(list)
                for i, gene_list in enumerate(gene_lists):
                    for gene in gene_list:
                        if gene not in record:
                            record.append(gene)
                            cluster = self.dbeta_info[self.dbeta_info["gene"] == gene][
                                "cluster"
                            ].values[0]
                            data[str(cluster)].append(GeneRank(gene, i))
                if len(data) == 0:
                    continue
                category_dict[model] = data
            for model, cluster_dict in category_dict.items():
                for cluster, genes in cluster_dict.items():
                    lowest_rank = min([gene.getRank() for gene in genes])
                    if out_format == "json":
                        result_dict[model][cluster] = [
                            gene.getGene()
                            for gene in genes
                            if gene.getRank() == lowest_rank
                        ]
                    elif out_format == "txt":
                        result_set.update(
                            [
                                gene.getGene()
                                for gene in genes
                                if gene.getRank() == lowest_rank
                            ]
                        )
        elif mode == "max":
            for model, gene_lists in gene_dict.items():
                data = defaultdict(list)
                for gene in gene_lists[-1]:
                    cluster = self.dbeta_info[self.dbeta_info["gene"] == gene][
                        "cluster"
                    ].values[0]
                    if out_format == "json":
                        data[str(cluster)].append(gene)
                    elif out_format == "txt":
                        result_set.add(gene)
                if len(data) == 0:
                    continue
                result_dict[model] = data
        elif isinstance(mode, int):
            data = defaultdict(list)
            gene_set_tmp = set()
            for _, gene_lists in gene_dict.items():
                for gene_list in gene_lists:
                    if len(gene_list) == mode:
                        for gene in gene_list:
                            gene_set_tmp.add(gene)
                    # filtered_gene_list.append(gene_list)
            for gene in gene_set_tmp:
                cluster = self.dbeta_info[self.dbeta_info["gene"] == gene][
                    "cluster"
                ].values[0]
                if out_format == "json":
                    data[str(cluster)].append(gene)
                elif out_format == "txt":
                    result_set.add(gene)
            result_dict["best"] = data
        else:
            raise ValueError("mode should be either 'min' or 'max'")
        if out_format == "json":
            with open(out_path, "w") as f:
                json.dump(result_dict, f)
        elif out_format == "txt":
            with open(out_path, "w") as f:
                f.write("gene\n")
                for gene in result_set:
                    f.write(f"{gene}\n")

    def set_train_validate(self) -> None:
        train_df_tt = self.train_df[
            self.train_df["Unnamed: 0"].isin(self.dbeta_info["ID"])
            | (self.train_df["Unnamed: 0"] == "label")
        ]
        validate_df_tt = self.validate_df[
            self.validate_df["Unnamed: 0"].isin(self.dbeta_info["ID"])
            | (self.validate_df["Unnamed: 0"] == "label")
        ]
        self.X_train = train_df_tt.iloc[:-1, 1:].T.values.tolist()
        self.X_validate = validate_df_tt.iloc[:-1, 1:].T.values.tolist()
        self.y_train = self.train_df.iloc[-1, 1:].astype(int)
        self.y_validate = self.validate_df.iloc[-1, 1:].astype(int)

        self.X_validate__ = []
        self.y_validate__ = []
        for seed in range(self.Bagging_num):
            self._validate_df_tt = validate_df_tt
            split_df = self._balance_dataframe_with_labels(validate_df_tt, seed=seed)
            self.split_df = split_df
            self.X_validate__.append(split_df.iloc[:-1, :].T.values.tolist())
            self.y_validate__.append(split_df.iloc[-1, :].astype(int))

    def set_selection_models(self, selection_models: dict) -> None:
        self.selection_models = selection_models

    def set_grid_estimators(self, grid_estimators: dict) -> None:
        # used by RFE
        self.grid_estimators = grid_estimators

    def select_feature_sfs(
        self,
        out_path: str,
        tol: float = 0.001,
        step: int = 1,
        n_features_to_select: Literal["auto", "cluster"] = "auto",
    ) -> None:
        """
        Select features using Sequential Feature Selector (SFS) with the selection models.
        The selected features are saved in the out_path/sfs/selected_feature_names_sfs.txt file.

        Parameters:
            out_path (str): The path to save the selected features.
            tol (float): The tolerance for the stopping criteriq.
            step (int): The step size for the number of features to select.
            n_features_to_select (int): The number of features to select. If "auto", the number of features to select is the number of unique clusters.
        """
        for selection_model_name, selection_model in self.selection_models.items():
            logger.info(f"Training {selection_model_name} with SFS")
            num_unique_clusters = len(set(self.dbeta_info["cluster"]))
            if n_features_to_select == "auto":
                step_ = n_features_to_select
            elif n_features_to_select == "cluster":
                step_ = num_unique_clusters
            else:
                raise ValueError(
                    "n_features_to_select should be either 'auto' or 'cluster'"
                )
            while 1:
                sfs = SequentialFeatureSelector(
                    estimator=selection_model,
                    direction="forward",
                    scoring=make_scorer(f1_score),
                    tol=tol,
                    n_jobs=-1,
                    n_features_to_select=step_,
                )
                sfs.fit_transform(self.X_train, self.y_train)
                self._save_selected_features(
                    sfs, selection_model_name
                ) >> self._append_to_file(f"{out_path}", is_first_header=False)
                if step_ == "auto":
                    logger.info(
                        f"Training finished with {self.selected_clusters} clusters selected"
                    )
                    break
                if self.selected_clusters == num_unique_clusters:
                    logger.info(
                        f"Training finished with {self.selected_clusters} clusters selected"
                    )
                    break
                step_ += step
                logger.info(
                    f"Training {selection_model_name} with {self.selected_clusters} clusters selected"
                )

    def select_feature_rfe(
        self,
        train_out_path: str,
        validate_out_path: str,
        selected_feature_path: str,
        feature_range: tuple,
    ) -> None:
        for selection_model_name, selection_model in self.selection_models.items():
            for feature_count in range(*feature_range):
                rfe = RFE(estimator=selection_model, n_features_to_select=feature_count)

                X_train_rfe = rfe.fit_transform(self.X_train, self.y_train)
                X_validate_rfe = rfe.transform(self.X_validate)

                X_validate_rfe__ = []
                for X_validate__item in self.X_validate__:
                    X_validate_rfe__.append(rfe.transform(X_validate__item))

                self._save_selected_features(
                    rfe, selection_model_name
                ) >> self._append_to_file(selected_feature_path, is_first_header=False)

                for estimator_name, grid_estimator in self.grid_estimators.items():
                    grid_estimator.fit(X_train_rfe, self.y_train)
                    self._fit_predict_append(
                        grid_estimator.best_estimator_,
                        X_train_rfe,
                        self.y_train,
                        selection_model_name,
                        estimator_name,
                        feature_count,
                        "rfe.csv",
                        "roc_curve.csv",
                        train_out_path,
                    )

                    self._fit_predict_append(
                        grid_estimator.best_estimator_,
                        X_validate_rfe,
                        self.y_validate,
                        selection_model_name,
                        estimator_name,
                        feature_count,
                        "rfe.csv",
                        "roc_curve.csv",
                        validate_out_path,
                    )
                    metric_avg = {}

                    # not planing to plot roc curve for each bagging anyway
                    for i, (X_validate_rfe__i, y_validate__i) in enumerate(
                        zip(X_validate_rfe__, self.y_validate__)
                    ):
                        y_pred_on_X = grid_estimator.best_estimator_.predict(
                            X_validate_rfe__i
                        )
                        y_proba_on_X = grid_estimator.best_estimator_.predict_proba(
                            X_validate_rfe__i
                        )[:, 1]
                        result = self._predict(
                            y=y_validate__i,
                            y_pred_on_X=y_pred_on_X,
                            y_proba_on_X=y_proba_on_X,
                        )
                        for key, value in result["metrics"].items():
                            if key not in metric_avg:
                                metric_avg[key] = value
                            else:
                                metric_avg[key] += value
                    for key, value in metric_avg.items():
                        metric_avg[key] = value / self.Bagging_num
                    pd.DataFrame(
                        [
                            {
                                "selection_model": selection_model_name,
                                "train_model": estimator_name,
                                "features": feature_count,
                                "accuracy": metric_avg["accuracy"],
                                "recall": metric_avg["recall"],
                                "specificity": metric_avg["specificity"],
                                "precision": metric_avg["precision"],
                                "f1_score": metric_avg["f1_score"],
                                "AUC": metric_avg["AUC"],
                                "MCC": metric_avg["MCC"],
                                "fbeta2_score": metric_avg["fbeta2_score"],
                            }
                        ]
                    ) >> self._append_to_file(f"{validate_out_path}/rfe_avg.csv")

    @pipe_dec
    def _append_to_file(
        self, data: pd.DataFrame, file_name: str, is_first_header: bool = True
    ) -> None:
        if not os.path.isfile(file_name):
            data.to_csv(file_name, index=False, header=is_first_header)
        else:
            data.to_csv(file_name, index=False, mode="a", header=False)

    def _fit_predict_append(
        self,
        estimator,
        X,
        y,
        selection_model_name,
        estimator_name,
        feature_count,
        method_path,
        roc_curve_path,
        out_path,
    ):
        y_pred_on_X = estimator.predict(X)
        y_proba_on_X = estimator.predict_proba(X)[:, 1]
        result = self._predict(
            y=y,
            y_pred_on_X=y_pred_on_X,
            y_proba_on_X=y_proba_on_X,
        )
        pd.DataFrame(
            [
                {
                    "selection_model": selection_model_name,
                    "train_model": estimator_name,
                    "features": feature_count,
                    **result["metrics"],
                }
            ]
        ) >> self._append_to_file(f"{out_path}/{method_path}")
        pd.DataFrame(
            [
                {
                    "selection_model": selection_model_name,
                    "train_model": estimator_name,
                    "features": feature_count,
                    **result["roc_curve"],
                }
            ]
        ) >> self._append_to_file(f"{out_path}/{roc_curve_path}")

    def _predict(
        self,
        y: pd.DataFrame,
        y_pred_on_X: pd.DataFrame,
        y_proba_on_X: pd.DataFrame,
    ) -> dict:
        tn, fp, _, _ = confusion_matrix(y, y_pred_on_X).ravel()
        accuracy = accuracy_score(y, y_pred_on_X)
        recall = recall_score(y, y_pred_on_X, zero_division=0)
        if (tn + fp) == 0:
            print("Warning: No true negative found. Setting specificity to 0.")
            specificity = -1
        else:
            specificity = tn / (tn + fp)
        precision = precision_score(y, y_pred_on_X, zero_division=0)

        fpr, tpr, threshold = roc_curve(y, y_proba_on_X)
        roc_auc = auc(fpr, tpr)
        fpr = fpr.tolist()
        tpr = tpr.tolist()

        return {
            "metrics": {
                "accuracy": accuracy,
                "recall": recall,
                "specificity": specificity,
                "precision": precision,
                "f1_score": f1_score(y, y_pred_on_X),
                "AUC": roc_auc,
                "MCC": matthews_corrcoef(y, y_pred_on_X),
                "fbeta2_score": fbeta_score(y, y_pred_on_X, beta=2),
            },
            "roc_curve": {
                "fpr": fpr,
                "tpr": tpr,
                "threshold": threshold,
            },
        }

    def _save_selected_features(
        self,
        selector,
        selection_model_name: str,
    ) -> None:
        index = self.dbeta_info.columns.get_loc("gene")
        cluster_index = self.dbeta_info.columns.get_loc("cluster")
        selected_features = self.dbeta_info.iloc[selector.support_, index]
        self.selected_clusters = len(
            set(self.dbeta_info.iloc[selector.support_, cluster_index])
        )

        selected_feature_df = pd.DataFrame(
            {
                "selection_model_name": [selection_model_name],
                **{
                    f"selected_feature_{i}": [selected_feature]
                    for i, selected_feature in enumerate(selected_features)
                },
            }
        )
        return selected_feature_df

    def _balance_dataframe_with_labels(self, df: pd.DataFrame, seed=42) -> pd.DataFrame:
        """
        Balances the dataframe to ensure equal normal and tumor samples for each ID
        by randomly selecting columns based on the labels in the last row.

        Parameters:
            df (pd.DataFrame): The input dataframe with:
                - 'ID' column as the first column (string).
                - Columns for samples, where the last row contains labels (0 for normal, 1 for tumor).
            seed (int): The random seed for reproducibility.

        Returns:
            pd.DataFrame: A balanced dataframe with equal normal and tumor samples per ID,
                        retaining the label row.
        """
        random.seed(seed)

        labels = df.iloc[-1, 1:]
        data = df.iloc[:-1, 1:]
        self._labels = labels
        self._data = data
        normal_columns = [col for col in data.columns if int(round(labels[col])) == 0]
        self.normal_columns = normal_columns
        tumor_columns = [col for col in data.columns if int(round(labels[col])) == 1]
        self.tumor_columns = tumor_columns
        min_columns = min(len(normal_columns), len(tumor_columns))
        self._min_columns = min_columns
        if len(normal_columns) > len(tumor_columns):
            selected_normal_columns = random.sample(normal_columns, min_columns)
            balanced_columns = selected_normal_columns + tumor_columns
        else:
            selected_tumor_columns = random.sample(tumor_columns, min_columns)
            balanced_columns = normal_columns + selected_tumor_columns
        self.balanced_columns = balanced_columns
        balanced_data = data[balanced_columns]
        balanced_labels = labels[balanced_columns]

        balanced_df = pd.concat([balanced_data, balanced_labels.to_frame().T])

        return balanced_df
