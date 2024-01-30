"""plotting results"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from seaborn.palettes import color_palette
from matplotlib import gridspec

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
    long_metric_df = metric_df.melt(id_vars=["dimension"], var_name="Method")

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


def frobenius_sub_boxplot(data: pd.DataFrame, subfig_axs=None, subtitle: str = None):
    """boxplot function for subplots"""
    if subfig_axs:
        axs = subfig_axs.subplots(ncols=3, nrows=1)
    else:
        _, axs = plt.subplots(ncols=3, nrows=1, figsize=(12, 8), layout="tight")
    palette = color_palette("colorblind", len(data["Method"].unique()))
    # Set y-axis label for the leftmost plot
    axs[0].set_ylabel("Frobenius norm")

    dims = [50, 250, 750]
    counter = -1
    for ax in axs:
        counter += 1
        plt.setp(ax.collections, alpha=0.6)
        plt.setp(ax.lines, alpha=0.6)
        sns.boxplot(
            data=data[data["dimension"] == dims[counter]],
            x="Method",
            y="value",
            ax=ax,
            palette=palette,
        )
        ax.set_ylabel("")
        ax.set_xlabel(f"dimension = {dims[counter]}")
        ax.set(xticklabels=[])
        ax.tick_params(bottom=False)
    axs[0].set_ylabel("Frobenius norm")
    if subtitle:
        axs.set_title(subtitle)
    # plt.show()


def auc_sub_pointplot(
    data: pd.DataFrame, show_legend: bool = False, subfig_axs=None, subtitle: str = None
):
    """pointplot function for subplots"""
    if subfig_axs:
        axs = subfig_axs.subplots()
    else:
        _, axs = plt.subplots()
    palette = color_palette("colorblind", len(data["Method"].unique()))
    g = sns.pointplot(
        x="dimension",
        y="value",
        hue="Method",
        errorbar="sd",
        data=data,
        ax=axs,
        dodge=True,
        linewidth=2,
        palette=palette,
    )
    g.set(ylabel="AUC")
    axs.set(ylim=(0.5, 1))
    plt.setp(axs.collections, alpha=0.6)
    g.get_legend().set_visible(False)
    if show_legend:
        g.get_legend().set_visible(True)
    if subtitle:
        axs.set_title(subtitle)
    # plt.show()


# frobenius_sub_boxplot(long_general_cubic_frobenius_df)
# auc_sub_pointplot(long_general_cubic_auc_df, show_legend=True)

location = (
    "/Users/kgoebler/hume_revisions/mixed_hidim_graphs/paper/High-Dimensional Mixed Graphs EJS/Figures/",
)
plotname = ("simulation_results.pdf",)

fig = plt.figure(layout="constrained", figsize=(10, 4))
subfigs = fig.subfigures(nrows=2, ncols=2, wspace=0.07)

axs_left_top = auc_sub_pointplot(
    data=long_general_auc_df, subfig_axs=subfigs[0, 0], subtitle=r"$f(x) = x$"
)
axs_right_top = auc_sub_pointplot(
    data=long_general_cubic_auc_df,
    subfig_axs=subfigs[0, 1],
    subtitle=r"$f(x) = x^{1/3}$",
)
axs_left_bottom = frobenius_sub_boxplot(
    data=long_general_frobenius_df, subfig_axs=subfigs[1, 0]
)
axs_right_bottom = frobenius_sub_boxplot(
    data=long_general_cubic_frobenius_df, subfig_axs=subfigs[1, 1]
)

handles, labels = fig.get_axes()[0].get_legend_handles_labels()

alt_handles = color_palette(
    "colorblind", len(long_general_frobenius_df["Method"].unique())
)

fig.legend(handles, labels, loc=(0.48, 0.46))

plt.show()
fig.savefig(location + plotname)

sys.exit()


# Plotting
def plot_simulation_results(
    result_data: list[pd.DataFrame],
    frobenius_ylim: tuple[float, float] = (2, 12),
    save_fig: bool = False,
    location: str = "/Users/kgoebler/hume_revisions/mixed_hidim_graphs/paper/High-Dimensional Mixed Graphs EJS/Figures/",
    plotname: str = "simulation_results.pdf",
):
    """plot simulation results"""
    fig, axs = plt.subplots(
        ncols=2, nrows=2, figsize=(12, 8), sharey="row", layout="tight"
    )
    g = sns.pointplot(
        x="dimension",
        y="value",
        hue="Method",
        errorbar="sd",
        data=result_data[0],
        ax=axs[0, 0],
        dodge=True,
    )
    g.set(ylabel="AUC")
    h = sns.pointplot(
        x="dimension",
        y="value",
        hue="Method",
        errorbar="sd",
        data=result_data[1],
        ax=axs[0, 1],
        dodge=True,
    )
    i = sns.boxplot(
        x="dimension",
        y="value",
        hue="Method",
        # errorbar="sd",
        data=result_data[2],
        ax=axs[1, 0],
        dodge=True,
    )

    i.set(ylabel="Frobenius Norm")
    j = sns.pointplot(
        x="dimension",
        y="value",
        hue="Method",
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

    # Add individual titles to each subplot
    axs[0, 0].set_title(r"$f(x) = x$")
    axs[0, 1].set_title(r"$f(x) = x^{1/3}$")
    axs[1, 0].set_title(r"$f(x) = x$")
    axs[1, 1].set_title(r"$f(x) = x^{1/3}$")

    plt.show()
    if save_fig:
        fig.savefig(location + plotname)


plot_simulation_results(
    result_data=binary_df_list, save_fig=True, plotname="binary.pdf"
)

plot_simulation_results(
    result_data=general_df_list, save_fig=True, plotname="general.pdf"
)


# f, axs = plt.subplots(ncols=2, nrows=2, figsize=(12, 8), sharey="row", layout="tight")
# g = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="Method",
#     errorbar="sd",
#     data=long_binary_auc_df,
#     ax=axs[0, 0],
#     dodge=True,
# )
# g.set(ylabel="AUC")
# h = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="Method",
#     errorbar="sd",
#     data=long_cubic_auc_df,
#     ax=axs[0, 1],
#     dodge=True,
# )
# i = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="Method",
#     errorbar="sd",
#     data=long_binary_frobenius_df,
#     ax=axs[1, 0],
#     dodge=True,
# )
# i.set(ylabel="Frobenius Norm")
# j = sns.pointplot(
#     x="dimension",
#     y="value",
#     hue="Method",
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
