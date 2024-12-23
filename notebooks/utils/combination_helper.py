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
)

def is_unique(lst):
    return len(lst) == len(set(lst))


def prepare_data(data, combination, label, cluster_num):
    X = []
    cluster_values = []
    for i in range(cluster_num):

        gene_values = data.loc[data[data.columns[0]].isin([combination[i]])]
        gene_values = gene_values.iloc[:, 1::].values.flatten().tolist()
        cluster_values.append(gene_values)

    for i in range(len(label)):
        sample = [cluster_values[j][i] for j in range(cluster_num)]
        X.append(sample)
    return X


def predict_comb(
    train_data, train_label, test_data, test_label, cluster_num, param_grids, final_gene
):
    combinations_list = list(combinations(final_gene["ID"], cluster_num))
    result = []
    models = {
        "XGBoost": XGBClassifier(),
        "RandomForest": RandomForestClassifier(),
        "SVM": SVC(),
        "Decision Tree": DecisionTreeClassifier(),
    }

    for model_name, model in models.items():
        # 取出所有可能的組合
        for combination in combinations_list:
            lst = []
            for i in range(cluster_num):
                lst.append(
                    int(
                        final_gene.loc[final_gene["ID"] == combination[i]][
                            "cluster"
                        ].iloc[0]
                    )
                )
            # 檢查是否在不同群
            if is_unique(lst):
                # 準備訓練 / 測試資料
                X_test = prepare_data(test_data, combination, test_label, cluster_num)
                X_train = prepare_data(
                    train_data, combination, train_label, cluster_num
                )

                # 模型參數調整
                param_grid = param_grids[model_name]
                grid_search = GridSearchCV(
                    estimator=model, param_grid=param_grid, cv=5, n_jobs=-1
                )
                grid_search.fit(X_train, train_label)
                print("Best Parameters:", grid_search.best_params_)
                print("Best Score:", grid_search.best_score_)
                model = grid_search.best_estimator_

                # 預測
                y_pred_train = model.predict(X_train)
                accuracy_train = accuracy_score(train_label, y_pred_train)
                # print("Train accuracy: ", round(accuracy_train, 2))
                y_pred = model.predict(X_test)
                accuracy = accuracy_score(test_label, y_pred)
                # print("Test accuracy: ", round(accuracy, 2))

                # 檢查是否過擬合
                if abs(accuracy_train - accuracy) > 0.1:
                    print("Train accuracy: ", round(accuracy_train, 2))
                    print("Test accuracy: ", round(accuracy, 2))
                    print(f"========={model_name} overfitting =========\n")
                else:
                    # 輸出結果
                    tn, fp, fn, tp = confusion_matrix(test_label, y_pred).ravel()
                    sensitivity = tp / (tp + fn)
                    specificity = tn / (tn + fp)
                    precision = precision_score(test_label, y_pred)
                    f1 = f1_score(test_label, y_pred)
                    mcc = matthews_corrcoef(test_label, y_pred)

                    sorted_combination = [x for _, x in sorted(zip(lst, combination))]
                    result.append(
                        [model_name]
                        + sorted_combination
                        + [
                            round(accuracy, 2),
                            round(sensitivity, 2),
                            round(specificity, 2),
                            round(precision, 2),
                            round(f1, 2),
                            round(mcc, 2),
                        ]
                    )
    id_columns = [f"ID{i+1}" for i in range(cluster_num)]
    result = pd.DataFrame(
        result,
        columns=["Model"]
        + id_columns
        + ["accuracy", "sensitivity", "specificity", "precision", "f1_score", "mcc"],
    )

    return result