from sklearn.model_selection import GridSearchCV
from sklearn.metrics import (
    accuracy_score,
    recall_score,
    precision_score,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
    auc,
    roc_curve,
    confusion_matrix,
)
import pandas as pd
import os

def set_parameters(model: object, param_grid: dict) -> object:
    grid_search = GridSearchCV(
        estimator=model, param_grid=param_grid, cv=5, n_jobs=-1, verbose=2
    )
    return grid_search

def append_to_file(file_name: str, data: pd.DataFrame) -> None:
        if not os.path.isfile(file_name):
            data.to_csv(file_name, index=False)
        else:
            data.to_csv(file_name, index=False, mode="a", header=False)

def record(
    best_estimator: object,
    X: pd.DataFrame,
    y: pd.Series,
    selection_model_name: str,
    train_model_name: str,
    feature_count: int,
    best_score: float,
    best_params: dict,
    filepath: str,
) -> None:
    
    y_pred_on_X = best_estimator.predict(X)

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

    append_to_file(
        f"{filepath}/rfe.csv",
        pd.DataFrame(
            [
                {
                    "selection_model": selection_model_name,
                    "train_model": train_model_name,
                    "features": feature_count,
                    "accuracy": accuracy,
                    "recall": recall,
                    "speficity": specificity,
                    "precision": precision,
                    "f1_score": f1_score(y, y_pred_on_X),
                    "J-index": recall + specificity - 1,
                    "AUC": roc_auc,
                    "MCC": matthews_corrcoef(y, y_pred_on_X),
                    "Mean cross-validated score": best_score,
                    "Best parameters": best_params,
                }
            ]
        ),
    )
    append_to_file(
        f"{filepath}/fpr_tpr.csv",
        pd.DataFrame(
            [
                {
                    "selection_model": selection_model_name,
                    "train_model": train_model_name,
                    "features": feature_count,
                    "fpr": fpr,
                    "tpr": tpr,
                    "AUC": roc_auc,
                }
            ]
        ),
    )
