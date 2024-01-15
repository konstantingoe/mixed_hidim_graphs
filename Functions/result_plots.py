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
general_import_string = ["general_50", "general_250", "general_750"]

import_string = binary_import_string + general_import_string

for result in import_string:
    with open(sim_path + result + ".json", "r+") as f:
        interim = json.load(f)
        result_dict[result] = json.loads(interim[0])

    with open(sim_path + result + "_cubic.json", "r+") as f:
        interim = json.load(f)
        result_dict[result + "_cubic"] = json.loads(interim[0])

binary_cubic_string = [x + "_cubic" for x in binary_import_string]
general_cubic_string = [x + "_cubic" for x in general_import_string]


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
long_cubic_auc_df = get_metric_long(result_dict, binary_cubic_string)

long_binary_frobenius_df = get_metric_long(
    result_dict, binary_import_string, which_metric="frobenius"
)
long_cubic_frobenius_df = get_metric_long(
    result_dict, binary_cubic_string, which_metric="frobenius"
)

long_general_auc_df = get_metric_long(result_dict, general_import_string)
long_general_cubic_auc_df = get_metric_long(result_dict, general_cubic_string)
long_general_frobenius_df = get_metric_long(
    result_dict, general_import_string, which_metric="frobenius"
)
long_general_cubic_frobenius_df = get_metric_long(
    result_dict, general_cubic_string, which_metric="frobenius"
)

binary_df_list = [
    long_binary_auc_df,
    long_cubic_auc_df,
    long_binary_frobenius_df,
    long_cubic_frobenius_df,
]

general_df_list = [
    long_general_auc_df,
    long_general_cubic_auc_df,
    long_general_frobenius_df,
    long_general_cubic_frobenius_df,
]


# Plotting
def plot_simulation_results(
    result_data: list[pd.DataFrame], frobenius_ylim: tuple[float, float] = (2, 12)
):
    """plot simulation results"""
    _, axs = plt.subplots(
        ncols=2, nrows=2, figsize=(12, 8), sharey="row", layout="tight"
    )
    g = sns.pointplot(
        x="dimension",
        y="value",
        hue="method",
        errorbar="sd",
        data=result_data[0],
        ax=axs[0, 0],
        dodge=True,
    )
    g.set(ylabel="AUC")
    h = sns.pointplot(
        x="dimension",
        y="value",
        hue="method",
        errorbar="sd",
        data=result_data[1],
        ax=axs[0, 1],
        dodge=True,
    )
    i = sns.pointplot(
        x="dimension",
        y="value",
        hue="method",
        errorbar="sd",
        data=result_data[2],
        ax=axs[1, 0],
        dodge=True,
    )
    i.set(ylabel="Frobenius Norm")
    j = sns.pointplot(
        x="dimension",
        y="value",
        hue="method",
        errorbar="sd",
        data=result_data[3],
        ax=axs[1, 1],
        dodge=True,
    )
    axs[0, 0].set(ylim=(0.5, 1))
    axs[1, 0].set(ylim=frobenius_ylim)

    for ax_row in axs:
        for ax in ax_row:
            plt.setp(ax.collections, alpha=0.6)
            plt.setp(ax.lines, alpha=0.6)

    sns.move_legend(g, "upper right", bbox_to_anchor=(1.2, 0))
    axs[0, 1].get_legend().set_visible(False)
    axs[1, 1].get_legend().set_visible(False)
    axs[1, 0].get_legend().set_visible(False)
    sns.color_palette("colorblind")

    plt.show()


plot_simulation_results(result_data=binary_df_list)

plot_simulation_results(result_data=general_df_list)


# f, axs = plt.subplots(ncols=2, nrows=2, figsize=(12, 8), sharey="row", layout="tight")
# g = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="method",
#     errorbar="sd",
#     data=long_binary_auc_df,
#     ax=axs[0, 0],
#     dodge=True,
# )
# g.set(ylabel="AUC")
# h = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="method",
#     errorbar="sd",
#     data=long_cubic_auc_df,
#     ax=axs[0, 1],
#     dodge=True,
# )
# i = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="method",
#     errorbar="sd",
#     data=long_binary_frobenius_df,
#     ax=axs[1, 0],
#     dodge=True,
# )
# i.set(ylabel="Frobenius Norm")
# j = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="method",
#     errorbar="sd",
#     data=long_cubic_frobenius_df,
#     ax=axs[1, 1],
#     dodge=True,
# )
# axs[0, 0].set(ylim=(0.5, 1))
# axs[1, 0].set(ylim=(2, 12))

# for ax_row in axs:
#     for ax in ax_row:
#         plt.setp(ax.collections, alpha=0.6)
#         plt.setp(ax.lines, alpha=0.6)

# sns.move_legend(g, "upper right", bbox_to_anchor=(1.2, 0))
# axs[0, 1].get_legend().set_visible(False)
# axs[1, 1].get_legend().set_visible(False)
# axs[1, 0].get_legend().set_visible(False)
# sns.color_palette("colorblind")

# plt.show()
