import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns


merge = pd.read_csv("../data/liquid_milk.csv")


print(merge)