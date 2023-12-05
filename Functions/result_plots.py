"""plotting results"""
import pandas as pd
import matplotlib.pyplot as plt

result_dict = {
    "Frobenius": [2.860, 2.845, 2.892, 2.890],
    "FPR": [0.015, 0.019, 0.042, 0.044],
    "TPR": [0.335, 0.343, 0.356, 0.356],
    "AUC": [0.881, 0.864, 0.816, 0.811],
    "Routine": [
        "Oracle",
        "Oracle nonparanormal",
        "Polyserial ML",
        "Polyserial nonparanormal",
    ],
}

df = pd.DataFrame(result_dict)

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
