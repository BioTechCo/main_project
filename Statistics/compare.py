import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns



liquid = []

liquid.append(pd.read_csv("../data/liquid_group_1.csv"))
liquid.append(pd.read_csv("../data/liquid_group_2.csv"))
liquid.append(pd.read_csv("../data/liquid_group_3.csv"))
liquid.append(pd.read_csv("../data/liquid_group_4.csv"))

milk = []

milk.append(pd.read_csv("../data/milk_group_1.csv"))
milk.append(pd.read_csv("../data/milk_group_2.csv"))
milk.append(pd.read_csv("../data/milk_group_3.csv"))
milk.append(pd.read_csv("../data/milk_group_4.csv"))


merged_df = []
for i in range(0,4):
    merged_df.append(pd.merge(liquid[i],milk[i], on='gene', how='left',suffixes=('_liquid', '_milk')))
    print(merged_df)
    merged_df[i].to_csv("../data/liquid_milk_group"+str(i+1)+".csv",index=False)
    
    # 資料
# 情況 1：dbeta_tissue > 0 & dbeta_liquid > 0
# 情況 2：dbeta_tissue < 0 & dbeta_liquid > 0
# 情況 3：dbeta_tissue < 0 & dbeta_liquid < 0
# 情況 4：dbeta_tissue > 0 & dbeta_liquid < 0

title = ["left up","left down","right down","right up"]

for i in range(0,4):
    plt.figure(figsize=(8, 6))
    genes = merged_df[i]['gene']
    dbeta_liquid = merged_df[i]['dbeta_liquid_liquid']
    dbeta_milk = merged_df[i]['dbeta_liquid_milk']

    # 繪製散點圖
    plt.scatter(dbeta_liquid, dbeta_milk, color='blue')

    # 加標籤
    for j in range(len(genes)):
        plt.text(dbeta_liquid[j], dbeta_milk[j], genes[j], fontsize=9)

    plt.title(title[i])
    plt.xlabel('dbeta_liquid')
    plt.ylabel('dbeta_milk')
    plt.grid(True)



# plt.show()
