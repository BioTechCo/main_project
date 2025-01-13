
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


class SimpleModel():
    def __init__(
        self,
        train_df: pd.DataFrame,
        test_df: pd.DataFrame,
        gene_dict: dict,
        dbeta_info: pd.DataFrame
    ):
        self.train_df = train_df
        self.test_df = test_df
        self.gene_dict = gene_dict
        self.dbeta_info = dbeta_info
        self.bagging_iterations = 10

    def setup_dbeta(self, model_name: str):
        self.dbeta_info = self.dbeta_info[self.dbeta_info["gene"].isin(
            self.gene_dict[model_name])]
        print(self.dbeta_info['gene'])

    def setup_train_test(self):
        train_df_tt = self.train_df[self.train_df["Unnamed: 0"].isin(
            self.dbeta_info["ID"]) | (self.train_df["Unnamed: 0"] == "label")]
        test_df_tt = self.test_df[self.test_df["Unnamed: 0"].isin(
            self.dbeta_info["ID"]) | (self.test_df["Unnamed: 0"] == "label")]
        self.X_train = train_df_tt.iloc[:-1, 1:].T.values.tolist()
        self.X_test = test_df_tt.iloc[:-1, 1:].T.values.tolist()
        self.y_train = self.train_df.iloc[-1, 1:].astype(int)
        self.y_test = self.test_df.iloc[-1, 1:].astype(int)

        self.X_test__ = []
        self.y_test__ = []
        for seed in range(self.bagging_iterations):
            self._test_df_tt = test_df_tt
            split_df = self._balance_dataframe_with_labels(
                test_df_tt, seed=seed)
            self.split_df = split_df
            self.X_test__.append(split_df.iloc[:-1, :].T.values.tolist())
            self.y_test__.append(split_df.iloc[-1, :].astype(int))

    def train(
        self,
        estimator_name: str = None,
        estimator: object = None,
        train_out_path: str = None,
        test_out_path: str = None,
        is_grid_search: bool = False,
    ):
        estimator.fit(self.X_train, self.y_train)
        if is_grid_search:
            best_estimator = estimator.best_estimator_
        else:
            best_estimator = estimator
        self._fit_predict_append(
            self.X_train,
            self.y_train,
            best_estimator,
            f"{estimator_name}_metrics",
            f"{estimator_name}_roc_curve",
            train_out_path
        )
        self._fit_predict_append(
            self.X_test,
            self.y_test,
            best_estimator,
            f"{estimator_name}_metrics",
            f"{estimator_name}_roc_curve",
            test_out_path
        )
        for i in range(self.bagging_iterations):
            self._fit_predict_append(
                self.X_test__[i],
                self.y_test__[i],
                best_estimator,
                f"{estimator_name}_metrics_{i}",
                f"{estimator_name}_roc_curve_{i}",
                test_out_path
            )

    def _fit_predict_append(
        self,
        X,
        y,
        estimator,
        evaluation_metric_path,
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
                    "accuracy": metrics["accuracy"],
                    "recall": metrics["recall"],
                    "specificity": metrics["specificity"],
                    "precision": metrics["precision"],
                    "f1_score": metrics["f1_score"],
                    "AUC": metrics["AUC"],
                    "MCC": metrics["MCC"],
                }
            ]
        ) >> self._append_to_file(f"{out_path}/{evaluation_metric_path}.csv")
        pd.DataFrame(
            [
                {
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

    def _balance_dataframe_with_labels(
            self,
            df: pd.DataFrame,
            seed=42
    ) -> pd.DataFrame:
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

    @pipe_dec
    def _append_to_file(self, data: pd.DataFrame, file_name: str, is_first_header: bool = True) -> None:
        if not os.path.isfile(file_name):
            data.to_csv(file_name, index=False, header=is_first_header)
        else:
            data.to_csv(file_name, index=False, mode="a", header=False)
