import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns



# d1_values = normalized_train.loc[normalized_train['Unnamed: 0'] == 'cg07792478']  MIR124-2
# d2_values = normalized_train.loc[normalized_train['Unnamed: 0'] == 'cg14231297']  ZSCAN18
# d3_values = normalized_train.loc[normalized_train['Unnamed: 0'] == 'cg04574090']  KCNE3


liquid = pd.read_csv("../data/breast_liquid_dbeta.csv")
# liquid = pd.read_csv("../data/breast_liquid_dbeta_all.csv")
tissue = pd.read_csv("../data/breast_train_tissue_dbeta.csv")

gene_list_name = 'gene'

# liquid['dbeta'] = liquid['dbeta'].abs()
# tissue['dbeta'] = tissue['dbeta'].abs()

merged_df = pd.merge(tissue, liquid, on='gene', how='inner',suffixes=('_tissue', '_liquid'))

# merged_df = pd.merge(tissue, liquid, on='ID', how='inner',suffixes=('_tissue', '_liquid'))

print(merged_df)
print()

plt.figure(figsize=(10, 6))
plt.subplot(2, 2, 1)
sns.histplot(merged_df['dbeta_tissue'], kde=True, color='blue')
plt.title('dbeta_tissue Histogram')
plt.xlabel('dbeta_tissue')
plt.ylabel('Frequency')

plt.subplot(2, 2, 2)
sns.histplot(merged_df['dbeta_liquid'], kde=True, color='green')
plt.title('dbeta_liquid Histogram')
plt.xlabel('dbeta_liquid')
plt.ylabel('Frequency')


scaler = MinMaxScaler()
merged_df[['dbeta_tissue_normalized', 'dbeta_liquid_normalized']] = scaler.fit_transform(merged_df[['dbeta_tissue', 'dbeta_liquid']])

print(merged_df.corr(method='pearson',numeric_only=True))
print(merged_df.corr(method='spearman',numeric_only=True))


plt.subplot(2, 2, 3)
sns.histplot(merged_df['dbeta_tissue_normalized'], kde=True, color='blue')
plt.title('dbeta_tissue Histogram')
plt.xlabel('dbeta_tissue')
plt.ylabel('Frequency')

plt.subplot(2, 2, 4)
sns.histplot(merged_df['dbeta_liquid_normalized'], kde=True, color='green')
plt.title('dbeta_liquid Histogram')
plt.xlabel('dbeta_liquid')
plt.ylabel('Frequency')



plt.figure(figsize=(8, 6))

genes = merged_df[gene_list_name]
dbeta_tissue = merged_df['dbeta_tissue']
dbeta_liquid = merged_df['dbeta_liquid']

# 繪製散點圖
plt.scatter(dbeta_tissue, dbeta_liquid, color='blue')

# 加標籤
for i in range(len(genes)):
    plt.text(dbeta_tissue[i], dbeta_liquid[i], genes[i], fontsize=9)

plt.title('Gene Differential Methylation')
plt.xlabel('dbeta_tissue')
plt.ylabel('dbeta_liquid')
plt.grid(True)





# 資料
plt.figure(figsize=(8, 6))
genes = merged_df[gene_list_name]
dbeta_tissue = merged_df['dbeta_tissue_normalized']
dbeta_liquid = merged_df['dbeta_liquid_normalized']

# 繪製散點圖
plt.scatter(dbeta_tissue, dbeta_liquid, color='blue')

# 加標籤
for i in range(len(genes)):
    plt.text(dbeta_tissue[i], dbeta_liquid[i], genes[i], fontsize=9)

plt.title('Gene Differential Methylation')
plt.xlabel('dbeta_tissue')
plt.ylabel('dbeta_liquid')
plt.grid(True)


plt.axhline(y=0.68, color='red', linestyle='--', linewidth=1)
plt.text(0.5, 0.7, 'y = 0.68', color='red', fontsize=10, rotation=0)



# 情況 1：dbeta_tissue > 0 & dbeta_liquid > 0
merged_df[(merged_df['dbeta_tissue'] > 0) & (merged_df['dbeta_liquid'] > 0)].to_csv("../data/liquid_group_1.csv",index=False)

# 情況 2：dbeta_tissue > 0 & dbeta_liquid < 0
merged_df[(merged_df['dbeta_tissue'] > 0) & (merged_df['dbeta_liquid'] < 0)] .to_csv("../data/liquid_group_2.csv",index=False)

# 情況 3：dbeta_tissue < 0 & dbeta_liquid < 0
merged_df[(merged_df['dbeta_tissue'] < 0) & (merged_df['dbeta_liquid'] < 0)].to_csv("../data/liquid_group_3.csv",index=False)

# 情況 4：dbeta_tissue < 0 & dbeta_liquid > 0
merged_df[(merged_df['dbeta_tissue'] < 0) & (merged_df['dbeta_liquid'] > 0)].to_csv("../data/liquid_group_4.csv",index=False)


plt.show()
