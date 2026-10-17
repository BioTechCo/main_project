## init

train_df_file
訓練資料集
用於 section 1, 2, 3

dmps_files
用於 section 1

majority_out_path
用於 section 1, 2, 3, 4


## feature_selection

validate_df_file
rfe 驗證資料集

train_out_path 
訓練輸出

validate_out_path 
驗證輸出

dbeta_info_file
dbeta 資訊(有cluster)

training_param_file
訓練參數

### GO_distance_matrix
被 GO_distance_matrix.R 抓的 path

input_file
dbeta 資訊(無cluster)

base_out_dir
輸出資料夾

## clustering_visual

result_prefix
三種 clustering 結果的 prefix

dbeta_file
dbeta 資訊(無cluster)

bp_file
bp distance matrix

mf_file
mf distance matrix

cc_file
cc distance matrix

terms_count_file
每個分支的 term 數量

result_out_path
輸出資料夾

## simple model

dbeta_file
dbeta 資訊(有cluster)

selected_feature_file
選擇的 feature json 檔

example
{
    "best": {
        "4": [
            "SNORD115-10",
            "TMEM196",
            "MIR377",
            "PCK1",
            "PATE2",
            "RALYL"
        ],
        "2": [
            "ATP5G2",
            "C1orf150",
            "C1orf114",
            "PCDHB15",
            "ATG16L1",
            "DLC1",
            "DOC2A",
            "AIM2",
            "DCD"
        ],
        "3": [
            "SCARF1"
        ],
        "1": [
            "BHLHE23"
        ]
    }
}

train_out_path
訓練輸出

validate_out_path
驗證輸出

df_file
資料集

training_param_file
訓練參數

## others

preprocess.filtering.threshold
過濾閥值

dmp_before_dropna_shape_feature
dmp dropna 前的 feature 數量

dmp_after_dropna_shape_feature
dmp dropna 後的 feature 數量

delta_beta_avg_feature_num
delta beta 平均 feature 數量

delta_beta_avg_feature_num_remove_NaN
delta beta 平均 feature 數量(移除 NaN)

delta_beta_avg_feature_num_remove_NaN_join_dmp
delta beta 平均 feature 數量(移除 NaN, join dmp)