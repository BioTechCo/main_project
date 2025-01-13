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

from logging.config import dictConfig
from api.logger import logging_config
import logging

dictConfig(logging_config)

logger = logging.getLogger('train_helper')


def set_parameters(model: object, param_grid: dict) -> object:
    grid_search = GridSearchCV(
        estimator=model, param_grid=param_grid, cv=5, n_jobs=-1, verbose=2
    )
    return grid_search


def read_selected_features(file_name: str) -> dict:
    gene_dict = {}
    with open(file_name, "r") as f:
        for line in f:
            selected_feature_names = line.split(",")[1:]
            selected_feature_names[-1] = selected_feature_names[-1].strip()
            gene_dict[line.split(",")[0]] = selected_feature_names
    return gene_dict


class TrainHelper():
    def __init__(self, train_df: pd.DataFrame, test_df: pd.DataFrame, TSS_threshold: pd.DataFrame) -> None:
        self.train_df = train_df
        self.test_df = test_df
        self.TSS_threshold = TSS_threshold
        self.selection_model = {}
        self.grid_estimators = {}
        self.Bagging_num = 5

    def set_train_test(self) -> None:
        train_df_tt = self.train_df[self.train_df["Unnamed: 0"].isin(
            self.TSS_threshold["ID"]) | (self.train_df["Unnamed: 0"] == "label")]
        test_df_tt = self.test_df[self.test_df["Unnamed: 0"].isin(
            self.TSS_threshold["ID"]) | (self.test_df["Unnamed: 0"] == "label")]
        self.X_train = train_df_tt.iloc[:-1, 1:].T.values.tolist()
        self.X_test = test_df_tt.iloc[:-1, 1:].T.values.tolist()
        self.y_train = self.train_df.iloc[-1, 1:].astype(int)
        self.y_test = self.test_df.iloc[-1, 1:].astype(int)

        self.X_test__ = []
        self.y_test__ = []
        for seed in range(self.Bagging_num):
            self._test_df_tt = test_df_tt
            split_df = self._balance_dataframe_with_labels(
                test_df_tt, seed=seed)
            self.split_df = split_df
            self.X_test__.append(split_df.iloc[:-1, :].T.values.tolist())
            self.y_test__.append(split_df.iloc[-1, :].astype(int))

    def set_selection_models(self, selection_models: dict) -> None:
        self.selection_models = selection_models

    def set_grid_estimators(self, grid_estimators: dict) -> None:
        self.grid_estimators = grid_estimators

    def select_feature_sfs(self, TrainOutPath: str, tol: float = 0.001, n_features_to_select="auto") -> None:
        for selection_model_name, selection_model in self.selection_models.items():
            logger.debug(f"Training {selection_model_name} with SFS")
            target = len(set(self.TSS_threshold['cluster']))
            counter = target
            while 1:
                sfs = SequentialFeatureSelector(
                    estimator=selection_model,
                    direction="forward",
                    scoring=make_scorer(f1_score),
                    tol=tol,
                    n_jobs=-1,
                    n_features_to_select=counter,
                )
                sfs.fit_transform(self.X_train, self.y_train)
                sfs.transform(self.X_test)
                self._save_selected_features(
                    sfs, selection_model_name
                ) >> self._append_to_file(f"{TrainOutPath}/sfs/selected_feature_names_sfs.txt", is_first_header=False)
                if self.selected_clusters == target:
                    logger.info(
                        f"Training finished with {self.selected_clusters} clusters selected")
                    break
                counter += target
                logger.debug(
                    f"Training {selection_model_name} with {self.selected_clusters} clusters selected")

    def train_rfe(self, TrainOutPath: str, TestOutPath: str, feature_range: tuple) -> None:
        for selection_model_name, selection_model in self.selection_models.items():
            for feature_count in range(*feature_range):
                rfe = RFE(estimator=selection_model,
                          n_features_to_select=feature_count)

                X_train_rfe = rfe.fit_transform(self.X_train, self.y_train)
                X_test_rfe = rfe.transform(self.X_test)

                X_test_rfe__ = []
                for X_test__item in self.X_test__:
                    X_test_rfe__.append(rfe.transform(X_test__item))

                self._save_selected_features(
                    rfe, selection_model_name
                ) >> self._append_to_file(f"{TrainOutPath}/selected_feature_names.txt", is_first_header=False)

                for estimator_name, grid_estimator in self.grid_estimators.items():
                    grid_estimator.fit(X_train_rfe, self.y_train)
                    self._fit_predict_append(grid_estimator.best_estimator_, X_train_rfe, self.y_train, selection_model_name,
                                             estimator_name, feature_count, "rfe", "roc_curve", TrainOutPath)

                    self._fit_predict_append(grid_estimator.best_estimator_, X_test_rfe, self.y_test, selection_model_name,
                                             estimator_name, feature_count, "rfe", "roc_curve", TestOutPath)

                    for i, (X_test_rfe__i, y_test__i) in enumerate(zip(X_test_rfe__, self.y_test__)):
                        self._fit_predict_append(grid_estimator.best_estimator_, X_test_rfe__i, y_test__i, selection_model_name,
                                                 estimator_name, feature_count, f"rfe_{i}", f"roc_curve_{i}", TestOutPath)

    @pipe_dec
    def _append_to_file(self, data: pd.DataFrame, file_name: str, is_first_header: bool = True) -> None:
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
        out_path
    ):
        y_pred_on_X = estimator.predict(X)
        y_proba_on_X = estimator.predict_proba(X)[:, 1]
        metrics = self._predict(
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
                    "accuracy": metrics["accuracy"],
                    "recall": metrics["recall"],
                    "specificity": metrics["specificity"],
                    "precision": metrics["precision"],
                    "f1_score": metrics["f1_score"],
                    "AUC": metrics["AUC"],
                    "MCC": metrics["MCC"],
                }
            ]
        ) >> self._append_to_file(f"{out_path}/{method_path}.csv")
        pd.DataFrame(
            [
                {
                    "selection_model": selection_model_name,
                    "train_model": estimator_name,
                    "features": feature_count,
                    "fpr": metrics["fpr"],
                    "tpr": metrics["tpr"],
                    "threshold": metrics["threshold"],
                }
            ]
        ) >> self._append_to_file(f"{out_path}/{roc_curve_path}.csv")

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
            "accuracy": accuracy,
            "recall": recall,
            "specificity": specificity,
            "precision": precision,
            "f1_score": f1_score(y, y_pred_on_X),
            "AUC": roc_auc,
            "MCC": matthews_corrcoef(y, y_pred_on_X),
            "fbeta2_score": fbeta_score(y, y_pred_on_X, beta=2),
            "fpr": fpr,
            "tpr": tpr,
            "threshold": threshold,
        }

    def _save_selected_features(
        self,
        selector,
        selection_model_name: str,
    ) -> None:
        index = self.TSS_threshold.columns.get_loc("gene")
        cluster_index = self.TSS_threshold.columns.get_loc("cluster")
        selected_features = self.TSS_threshold.iloc[selector.support_, index]
        self.selected_clusters = len(
            set(self.TSS_threshold.iloc[selector.support_, cluster_index]))

        selected_feature_df = pd.DataFrame({
            "selection_model_name": [selection_model_name],
            **{f"selected_feature_{i}": [selected_feature] for i, selected_feature in enumerate(selected_features)}
        })
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
        normal_columns = [col for col in data.columns if int(
            round(labels[col])) == 0]
        self.normal_columns = normal_columns
        tumor_columns = [col for col in data.columns if int(
            round(labels[col])) == 1]
        self.tumor_columns = tumor_columns
        min_columns = min(len(normal_columns), len(tumor_columns))
        self._min_columns = min_columns
        if len(normal_columns) > len(tumor_columns):
            selected_normal_columns = random.sample(
                normal_columns, min_columns)
            balanced_columns = selected_normal_columns + tumor_columns
        else:
            selected_tumor_columns = random.sample(tumor_columns, min_columns)
            balanced_columns = normal_columns + selected_tumor_columns
        self.balanced_columns = balanced_columns
        balanced_data = data[balanced_columns]
        balanced_labels = labels[balanced_columns]

        balanced_df = pd.concat([balanced_data, balanced_labels.to_frame().T])

        return balanced_df
