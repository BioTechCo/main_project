import os
import pandas as pd

from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier

from sklearn.ensemble import RandomForestClassifier

from sklearn.svm import SVC

from sklearn.tree import DecisionTreeClassifier

from sklearn.metrics import (
    confusion_matrix,
    precision_score,
    accuracy_score,
    matthews_corrcoef,
    f1_score,
    recall_score,
    roc_curve,
    auc,

)
from itertools import product
from typing import Literal


class CombinationHelper:
    def __init__(
        self,
        train_df: pd.DataFrame,
        test_df: pd.DataFrame,
        dbeta: pd.DataFrame,
        gene_list: list[str],
    ) -> None:
        self.train_df = train_df
        self.test_df = test_df
        self.dbeta = dbeta
        self.gene_list = gene_list
        self.models = {}
        self.grid_estimators = {}

    def setup_dbeta(self):
        self.dbeta = self.dbeta[self.dbeta["gene"].isin(
            self.gene_list["gene"])]

    def setup_train_test(self):
        self.X_train = self.train_df[self.train_df['Unnamed: 0'].isin(
            self.dbeta['ID'])]
        self.X_test = self.test_df[self.test_df['Unnamed: 0'].isin(
            self.dbeta['ID'])]
        self.y_train = self.train_df.iloc[-1, 1:].astype(int)
        self.y_test = self.test_df.iloc[-1, 1:].astype(int)

    def setup_combinations(self):
        groups = self.dbeta['cluster'].unique()
        if len(groups) < 3:
            raise ValueError(
                "Not all four groups are present in the DataFrame.")

        genes_by_group = {
            group: self.dbeta[self.dbeta['cluster'] == group]['gene'].tolist() for group in groups}

        self.combs = [
            combination for combination in product(*genes_by_group.values())
            if len(set(combination)) == 4  # Ensure genes are distinct
        ]

    def setup_estimators(self):
        self.models = {
            "XGBoost": XGBClassifier(),
            "Random_Forest": RandomForestClassifier(),
            "SVM": SVC(),
            "Decision_Tree": DecisionTreeClassifier(),
        }

    def setup_grids(self):
        self.param_grids = {
            "XGBoost": {
                "learning_rate": [0.001, 0.01, 0.1],
                "n_estimators": [50, 100, 200],
                "max_depth": [2, 3, 4],
                "subsample": [0.6, 0.7, 0.8],
                "colsample_bytree": [0.6, 0.7, 0.8],
                "min_child_weight": [40, 50, 60],
            },
            "Random_Forest": {
                "n_estimators": [50, 100, 150, 200],
                "min_samples_split": [80, 90],
                "min_samples_leaf": [80, 90],
                "max_features": [5, 10, 20],
            },
            "SVM": {"kernel": ["rbf", "poly", "sigmoid"], "C": [50, 100, 200]},
            "Decision_Tree": {
                "min_samples_split": [70, 80, 90],
                "min_samples_leaf": [70, 80, 90],
                "max_features": [5, 10, 20],
            },
        }

    def setup_grid_estimator(self):
        for model_name, model in self.models.items():
            self.grid_estimators[model_name] = GridSearchCV(
                estimator=model, param_grid=self.param_grids[model_name], cv=5, n_jobs=-1, verbose=2
            )

    def get_ID(self, genes: tuple) -> tuple:
        return tuple(self.dbeta[self.dbeta['gene'].isin(genes)]['ID'])

    def prepare_data(self, type: Literal["train", "test"] = "train", genes: tuple = None) -> list:
        if type == "train":
            return self.X_train[self.X_train['Unnamed: 0'].isin(list(self.get_ID(genes)))].iloc[:, 1:].T.values.tolist()
        elif type == "test":
            return self.X_test[self.X_test['Unnamed: 0'].isin(list(self.get_ID(genes)))].iloc[:, 1:].T.values.tolist()
        else:
            raise ValueError("Invalid type. Must be 'train' or 'test'.")

    def append_to_file(self, file_name: str, data: pd.DataFrame) -> None:
        if not os.path.isfile(file_name):
            data.to_csv(file_name, index=False)
        else:
            data.to_csv(file_name, index=False, mode="a", header=False)

    def record(
        self,
        y_pred_on_X: pd.Series,
        y: pd.Series,
        model_name: str,
        best_score: float,
        best_params: dict,
        filepath: str,
        filename: str,
        combination: tuple
    ) -> None:

        tn, fp, _, _ = confusion_matrix(y, y_pred_on_X).ravel()
        accuracy = accuracy_score(y, y_pred_on_X)
        recall = recall_score(y, y_pred_on_X)
        if (tn + fp) == 0:
            print("Warning: No true negative found. Setting specificity to 0.")
            specificity = 0
        else:
            specificity = tn / (tn + fp)
        precision = precision_score(y, y_pred_on_X)
        fpr, tpr, _ = roc_curve(y, y_pred_on_X)
        roc_auc = auc(fpr, tpr)
        fpr = fpr.tolist()
        tpr = tpr.tolist()

        self.append_to_file(
            f"{filepath}/{filename}.csv",
            pd.DataFrame(
                [
                    {
                        "model": model_name,
                        "accuracy": accuracy,
                        "recall": recall,
                        "specificity": specificity,
                        "precision": precision,
                        "f1_score": f1_score(y, y_pred_on_X),
                        "AUC": roc_auc,
                        "MCC": matthews_corrcoef(y, y_pred_on_X),
                        "Mean cross-validated score": best_score,
                        "Best parameters": best_params,
                        **{f"ID{i+1}": gene for i, gene in enumerate(combination) if gene is not None}
                    }
                ]
            )
        )
        self.append_to_file(
            f"{filepath}/{filename}_fpr_tpr.csv",
            pd.DataFrame(
                [
                    {
                        "model": model_name,
                        "fpr": fpr,
                        "tpr": tpr,
                        "AUC": roc_auc,
                    }
                ]
            )
        )

    def train(self, train_folder: str, test_folder: str, filename: str, discard_overfitting: bool = False):
        for model_name, grid_estimator in self.grid_estimators.items():
            for combination in self.combs:
                prepared_X_train = self.prepare_data("train", combination)
                prepared_X_test = self.prepare_data("test", combination)
                grid_estimator.fit(prepared_X_train, self.y_train)

                best_params = grid_estimator.best_params_

                best_score = grid_estimator.best_score_
                y_pred_on_X_train = grid_estimator.best_estimator_.predict(
                    prepared_X_train)
                y_pred_on_X_test = grid_estimator.best_estimator_.predict(
                    prepared_X_test)

                if discard_overfitting:
                    train_acc = accuracy_score(self.y_train, y_pred_on_X_train)
                    test_acc = accuracy_score(self.y_test, y_pred_on_X_test)
                    if abs(train_acc - test_acc) > 0.1:
                        print("Train accuracy: ", round(train_acc, 6))
                        print("Test accuracy: ", round(test_acc, 6))
                        print(f"========={model_name} overfitting =========")
                        continue
                self.record(
                    y_pred_on_X_train,
                    self.y_train,
                    model_name,
                    best_score,
                    best_params,
                    train_folder,
                    filename,
                    combination
                )
                self.record(
                    y_pred_on_X_test,
                    self.y_test,
                    model_name,
                    best_score,
                    best_params,
                    test_folder,
                    filename,
                    combination
                )
