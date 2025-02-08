import os
import pandas as pd
import random
from .pipe import pipe_dec
from sklearn.metrics import (
    accuracy_score,
    recall_score,
    precision_score,
    f1_score,
    matthews_corrcoef,
    fbeta_score,
    roc_curve,
    auc,
    confusion_matrix,
)
from itertools import product
from logging.config import dictConfig
from api.logger import logging_config
import logging
from typing import Literal

dictConfig(logging_config)

logger = logging.getLogger("simple_model")


class SimpleModel:
    def __init__(
        self,
        train_df: pd.DataFrame,
        test_df: pd.DataFrame,
        gene_list: list,
        dbeta_info: pd.DataFrame,
    ):
        self.train_df = train_df
        self.test_df = test_df
        self.gene_list = gene_list
        self.dbeta_info = dbeta_info
        self.bagging_iterations = 10

    def setup_dbeta(self):
        self.dbeta_info = self.dbeta_info[self.dbeta_info["gene"].isin(self.gene_list)]

    def setup_train_test(self):
        """
        Prepare the training and testing dataframes by filtering out the genes
        that are not in the gene list.
        note that Xs are not yet transposed and converted to list
        this is done in the _prepare_data method where a gene list can be passed
        """
        train_df_tt = self.train_df[
            self.train_df["Unnamed: 0"].isin(self.dbeta_info["ID"])
            | (self.train_df["Unnamed: 0"] == "label")
        ]
        test_df_tt = self.test_df[
            self.test_df["Unnamed: 0"].isin(self.dbeta_info["ID"])
            | (self.test_df["Unnamed: 0"] == "label")
        ]
        self.pre_X_train = train_df_tt.iloc[:-1, :]
        self.pre_X_test = test_df_tt.iloc[:-1, :]
        self.y_train = self.train_df.iloc[-1, 1:].astype(int)
        self.y_test = self.test_df.iloc[-1, 1:].astype(int)

        self.pre_X_test__ = []
        self.y_test__ = []
        for seed in range(self.bagging_iterations):
            self._test_df_tt = test_df_tt
            split_df = self._balance_dataframe_with_labels(test_df_tt, seed=seed)
            self.split_df = split_df
            self.pre_X_test__.append(split_df.iloc[:-1, :])
            self.y_test__.append(split_df.iloc[-1, 1:].astype(int))

    def train(
        self,
        estimator_name: str = None,
        estimator: object = None,
        train_out_path: str = None,
        test_out_path: str = None,
        is_grid_search: bool = False,
    ):
        try:
            for combination in self.combs:
                logger.info(
                    f"Training for combination: {combination} with estimator: {estimator_name}"
                )
                X_train = self._prepare_data(self.pre_X_train, combination)
                estimator.fit(X_train, self.y_train)
                if is_grid_search:
                    self._estimator = estimator.best_estimator_
                else:
                    self._estimator = estimator
                self._fit_predict_append(
                    X_train,
                    self.y_train,
                    combination,
                    estimator_name,
                    self._estimator,
                    f"{estimator_name}_metrics.csv",
                    f"{estimator_name}_roc_curve.csv",
                    train_out_path,
                )

                X_test = self._prepare_data(self.pre_X_test, combination)
                self._fit_predict_append(
                    X_test,
                    self.y_test,
                    combination,
                    estimator_name,
                    self._estimator,
                    f"{estimator_name}_metrics.csv",
                    f"{estimator_name}_roc_curve.csv",
                    test_out_path,
                )

                # Bagging
                metric_avg = {}
                for i in range(self.bagging_iterations):
                    X_test__ = self._prepare_data(self.pre_X_test__[i], combination)
                    y_pred_on_X = self._estimator.predict(X_test__)
                    y_proba_on_X = self._estimator.predict_proba(X_test__)[:, 1]
                    result = self._predict(
                        y=self.y_test__[i],
                        y_pred_on_X=y_pred_on_X,
                        y_proba_on_X=y_proba_on_X,
                    )
                    # not drawing the roc curve for bagging
                    # pd.DataFrame(
                    #     [
                    #         {
                    #             "estimator_name": estimator_name,
                    #             **{f"gene_{i}": gene for i, gene in enumerate(combination)},
                    #             **result["roc_curve"]
                    #         }
                    #     ]
                    # ) >> self._append_to_file(f"{out_path}/{estimator_name}_roc_curve_avg_{i}.csv")

                    for key, value in result["metrics"].items():
                        if key not in metric_avg:
                            metric_avg[key] = value
                        else:
                            metric_avg[key] += value
                for key, value in metric_avg.items():
                    metric_avg[key] = value / self.bagging_iterations
                pd.DataFrame(
                    [
                        {
                            "estimator_name": estimator_name,
                            **{f"gene_{i}": gene for i, gene in enumerate(combination)},
                            **metric_avg,
                        }
                    ]
                ) >> self._append_to_file(
                    f"{test_out_path}/{estimator_name}_metrics_avg.csv"
                )

        except Exception as e:
            logger.error(f"Error in testing: {e}")
            raise e

    def setup_combinations(self):
        """
        Set up the combinations of genes based on the groups in the dbeta_info dataframe.
        """
        groups = self.dbeta_info["cluster"].unique()
        genes_by_group = {
            group: self.dbeta_info[self.dbeta_info["cluster"] == group]["gene"].tolist()
            for group in groups
        }

        self.combs = [
            combination
            for combination in product(*genes_by_group.values())
            if len(set(combination)) >= 1  # Ensure there is at least one gene
        ]

    def _prepare_data(self, df: pd.DataFrame, genes: tuple = None) -> list:
        if genes is None:
            return df.iloc[:, 1:].T.values.tolist()
        else:
            return (
                df[df["Unnamed: 0"].isin(list(self._get_ID(genes)))]
                .iloc[:, 1:]
                .T.values.tolist()
            )

    def _get_ID(self, genes: tuple) -> tuple:
        # logger.debug(f"Getting ID for genes: {tuple(self.dbeta_info[self.dbeta_info['gene'].isin(genes)]['ID'])}")
        return tuple(self.dbeta_info[self.dbeta_info["gene"].isin(genes)]["ID"])

    def _fit_predict_append(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        combination: tuple,
        estimator_name: str,
        estimator: object,
        evaluation_metric_path: str,
        roc_curve_path: str,
        out_path: str,
    ):
        """
        A helper function to fit the estimator, predict on the data, and append the results to a file.
        """
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
                    "estimator_name": estimator_name,
                    **{f"gene_{i}": gene for i, gene in enumerate(combination)},
                    **result["metrics"],
                }
            ]
        ) >> self._append_to_file(f"{out_path}/{evaluation_metric_path}")
        pd.DataFrame(
            [
                {
                    "estimator_name": estimator_name,
                    **{f"gene_{i}": gene for i, gene in enumerate(combination)},
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

    def _balance_dataframe_with_labels(self, df: pd.DataFrame, seed=42) -> pd.DataFrame:
        """
        Balances the dataframe to ensure equal normal and tumor samples for each ID
        by randomly selecting columns based on the labels in the last row.

        Parameters:
            df (pd.DataFrame): The input dataframe with:
                - 'ID' column as the first column (string), where the last row being 'label'.
                - Columns for samples, where the last row contains labels (0 for normal, 1 for tumor).
            seed (int): The random seed for reproducibility.

        Returns:
            pd.DataFrame: A balanced dataframe with equal normal and tumor samples per ID,
                        retaining the label row.
        """
        random.seed(seed)
        id = df.iloc[:, 0]
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

        balanced_df = pd.concat([id, balanced_df], axis=1)
        return balanced_df

    @pipe_dec
    def _append_to_file(
        self, data: pd.DataFrame, file_name: str, is_first_header: bool = True
    ) -> None:
        # check if folder exists
        if not os.path.exists(os.path.dirname(file_name)):
            os.makedirs(os.path.dirname(file_name))
        if not os.path.isfile(file_name):
            data.to_csv(file_name, index=False, header=is_first_header)
        else:
            data.to_csv(file_name, index=False, mode="a", header=False)
