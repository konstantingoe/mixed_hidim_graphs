"""plotting results"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sim_path = "/Users/kgoebler/hume_revisions/mixed_hidim_graphs/simulation_results/"
sim_runs = 100
result_dict = {}

binary_import_string = ["binary_50", "binary_250", "binary_750"]
for result in binary_import_string:
    with open(sim_path + result + ".json", "r+") as f:
        interim = json.load(f)
        result_dict[result] = json.loads(interim[0])

    with open(sim_path + result + "_cubic.json", "r+") as f:
        interim = json.load(f)
        result_dict[result + "_cubic"] = json.loads(interim[0])

cubic_string = [x + "_cubic" for x in binary_import_string]


def get_metric_long(result_dict: dict, key_name_list: list, which_metric: str = "auc"):
    """get auc from result dict"""

    metric_list = []
    for string in key_name_list:
        temp = pd.DataFrame(result_dict[string]).loc[which_metric, :]
        reshaper = pd.DataFrame(
            {temp.index[i]: temp.iloc[i] for i in range(len(temp.index))}
        )
        metric_list.append(reshaper)

    metric_df = pd.concat(metric_list, axis=0)

    metric_df.insert(
        loc=0,
        column="dimension",
        value=np.concatenate(
            [([i] * reshaper.shape[0]) for i in [50, 250, 750]], axis=0
        ),
    )
    long_metric_df = metric_df.melt(id_vars=["dimension"], var_name="method")

    return long_metric_df


long_binary_auc_df = get_metric_long(result_dict, binary_import_string)
long_cubic_auc_df = get_metric_long(result_dict, cubic_string)

long_binary_frobenius_df = get_metric_long(
    result_dict, binary_import_string, which_metric="frobenius"
)
long_cubic_frobenius_df = get_metric_long(
    result_dict, cubic_string, which_metric="frobenius"
)

# Plotting
f, axs = plt.subplots(ncols=2, nrows=2, figsize=(12, 8), sharey="row", layout="tight")
g = sns.pointplot(
    x="dimension",
    y="value",
    hue="method",
    errorbar="sd",
    data=long_binary_auc_df,
    ax=axs[0, 0],
)
g.set(ylabel="AUC")
h = sns.pointplot(
    x="dimension",
    y="value",
    hue="method",
    errorbar="sd",
    data=long_cubic_auc_df,
    ax=axs[0, 1],
)
i = sns.pointplot(
    x="dimension",
    y="value",
    hue="method",
    errorbar="sd",
    data=long_binary_frobenius_df,
    ax=axs[1, 0],
)
i.set(ylabel="Frobenius Norm")
j = sns.pointplot(
    x="dimension",
    y="value",
    hue="method",
    errorbar="sd",
    data=long_cubic_frobenius_df,
    ax=axs[1, 1],
)
axs[0, 0].set(ylim=(0.5, 1))
axs[1, 0].set(ylim=(2, 6))

sns.move_legend(g, "upper right", bbox_to_anchor=(1.2, 0))
axs[0, 1].get_legend().set_visible(False)
axs[1, 1].get_legend().set_visible(False)
axs[1, 0].get_legend().set_visible(False)
sns.color_palette("colorblind")

plt.show()
