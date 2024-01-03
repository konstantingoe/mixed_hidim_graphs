"""plotting results"""
import json

import pandas as pd
import matplotlib.pyplot as plt

sim_path = "mixed_hidim_graphs/simulation_results/"


with open(sim_path + "binary_50.json", 'r+') as f:
    binary_50 = json.load(f)
    binary_50 = json.loads(binary_50[0])


my_df = pd.DataFrame(binary_50)

my_df.explode(["oracle","fan", "poly", "mle"])



d = [50, 250, 750]
# latent Gaussian
oracle = {
    "d": [50, 250, 750],
    "Frobenius": [2.860, 3.201, 11.207],
    "sd": [0.091, 0.109, 0.140],
}

oracle_df = pd.DataFrame(oracle)


fig, ax = plt.subplots()
ax.errorbar(oracle_df["d"], oracle_df["Frobenius"], yerr=oracle_df["sd"], linestyle="")
plt.show()

