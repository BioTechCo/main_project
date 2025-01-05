import os
import pandas as pd

from itertools import combinations
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
from typing import List
from itertools import product
from typing import Literal


class CombinationHelper:
    def __init__(self,
                 y_train: pd.Series,
                 X_train: pd.DataFrame,
                 y_test: pd.Series,
                 X_test: pd.DataFrame,
                 dbeta: pd.DataFrame,
                 gene_list: List[str],
                 param_grids: dict) -> None:
        self.y_train = y_train
        self.X_train = X_train
        self.y_test = y_test
        self.X_test = X_test
        self.dbeta = dbeta
        self.gene_list = gene_list
        self.param_grids = param_grids

        self.dbeta = self.dbeta[self.dbeta["gene"].isin(
            self.gene_list["gene"])]
        self.X_train = self.X_train[self.X_train['Unnamed: 0'].isin(
            self.dbeta['ID'])]
        self.X_test = self.X_test[self.X_test['Unnamed: 0'].isin(
            self.dbeta['ID'])]

        groups = dbeta['cluster'].unique()
        if len(groups) < 3:
            raise ValueError(
                "Not all four groups are present in the DataFrame.")

        genes_by_group = {
            group: dbeta[dbeta['cluster'] == group]['gene'].tolist() for group in groups}

        self.combs = [
            combination for combination in product(*genes_by_group.values())
            if len(set(combination)) == 4  # Ensure genes are distinct
        ]

        self.models = {
            "XGBoost": XGBClassifier(),
            "Random_Forest": RandomForestClassifier(),
            "SVM": SVC(),
            "Decision_Tree": DecisionTreeClassifier(),
        }

    def prepare_data(self, type: Literal["train", "test"] = "train", genes: tuple = None) -> list:
        if type == "train":
            return self.X_train[self.X_train['Unnamed: 0'].isin(genes)].iloc[:, 1:].T.values.tolist()
        elif type == "test":
            return self.X_test[self.X_test['Unnamed: 0'].isin(genes)].iloc[:, 1:].T.values.tolist()
        else:
            raise ValueError("Invalid type. Must be 'train' or 'test'.")

    def set_parameters(self, model: object, model_name: str) -> object:
        return GridSearchCV(
            estimator=model, param_grid=self.param_grids[model_name], cv=5, n_jobs=-1, verbose=2
        )

    def append_to_file(file_name: str, data: pd.DataFrame) -> None:
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
            ),
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
            ),
        )

    def train(self, folder: str, filename: str, discard_overfitting: bool = False):
        for model_name, model in self.models.items():
            for combination in self.combs:
                y_train = self.y_train
                y_test = self.y_test
                X_train = self.prepare_data("train", combination)
                X_test = self.prepare_data("test", combination)
                model = self.set_parameters(model, model_name)
                model.fit(X_train, y_train)

                best_params = model.best_params_

                best_score = model.best_score_
                y_pred_on_X_train = model.best_estimator.predict(X_train)
                y_pred_on_X_test = model.best_estimator.predict(X_test)

                if discard_overfitting:
                    train_acc = accuracy_score(y_train, y_pred_on_X_train)
                    test_acc = accuracy_score(y_test, y_pred_on_X_test)
                    if abs(train_acc - test_acc) > 0.1:
                        print("Train accuracy: ", round(train_acc, 6))
                        print("Test accuracy: ", round(test_acc, 6))
                        print(f"========={model_name} overfitting =========")
                        continue
                self.record(
                    y_pred_on_X_train,
                    y_train,
                    model_name,
                    best_score,
                    best_params,
                    folder,
                    filename
                )
                self.record(
                    y_pred_on_X_test,
                    y_test,
                    model_name,
                    best_score,
                    best_params,
                    folder,
                    filename
                )
