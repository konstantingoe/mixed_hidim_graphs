"""plotting results"""
import json

import pandas as pd
import matplotlib.pyplot as plt

sim_path = "mixed_hidim_graphs/simulation_results/"


with open(sim_path + "binary_50.json", 'r+') as f:
    binary_50 = json.load(f)
    binary_50 = json.loads(binary_50[0])


my_df = pd.DataFrame(binary_50)


exploded = my_df.transpose().explode(["auc","frobenius", "tpr", "fpr"])
exploded["routine"] = exploded.index

exploded.boxplot(column=["auc", "frobenius"], by="routine", figsize=(12, 8))

plt.show()