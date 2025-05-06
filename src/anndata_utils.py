import pandas as pd
import matplotlib.pyplot as plt
import anndata
import seaborn as sns

def plot_histograms_by_group(adata, var_names, obs_names, bins=20, legend_loc="upper right"):
    """
    Plot histograms for specified features (var_names) stratified by sample groups (obs_names).

    Parameters:
        adata (AnnData): Input AnnData object.
        var_names (list of str): List of feature names to plot.
        obs_names (str): Column name in adata.obs to stratify the histogram by.
        bins (int): Number of bins for the histograms. Defaults to 20.
        legend_loc (str): Location for the legend. Defaults to "upper right".

    Returns:
        None
    """
    if adata.is_view:
        adata = adata.copy()

    if obs_names not in adata.obs.columns:
        raise ValueError(f"'{obs_names}' not found in adata.obs.")

    # Ensure stratification column is categorical
    adata.obs[obs_names] = adata.obs[obs_names].fillna("Unknown").astype("category")

    for var_name in var_names:
        if var_name not in adata.var_names:
            print(f"Warning: '{var_name}' not found in adata.var_names. Skipping.")
            continue

        # Extract data
        data = pd.DataFrame({
            obs_names: adata.obs[obs_names],
            var_name: adata[:, var_name].X.flatten()
        })

        # Plot histogram stratified by group
        plt.figure(figsize=(8, 6))
        histplot = sns.histplot(
            data=data,
            x=var_name,
            hue=obs_names,
            bins=bins,
            kde=True,
            stat="count",
            multiple="stack",
            alpha=0.3
        )

        # Explicitly fetch and set legend
        handles, labels = histplot.get_legend_handles_labels()
        if handles and labels:
            plt.legend(handles, labels, title=obs_names, loc=legend_loc)
        else:
            print(f"Warning: No legend found for '{var_name}' stratified by '{obs_names}'.")

        # Add labels and title
        plt.title(f"Histogram of '{var_name}' Stratified by '{obs_names}'")
        plt.xlabel(var_name)
        plt.ylabel("Count")
        plt.tight_layout()
        plt.show()

def plot_boxplots_by_group(adata, var_names, obs_names):
    """
    Plot boxplots for specified features (var_names) stratified by sample groups (obs_names).

    Parameters:
        adata (AnnData): Input AnnData object.
        var_names (list of str): List of feature names to plot.
        obs_names (str): Column name in adata.obs to stratify the boxplot by.

    Returns:
        None
    """
    if adata.is_view:
        adata = adata.copy()

    if obs_names not in adata.obs.columns:
        raise ValueError(f"'{obs_names}' not found in adata.obs.")

    # Ensure stratification column is categorical
    adata.obs[obs_names] = adata.obs[obs_names].fillna("Unknown").astype("category")

    for var_name in var_names:
        if var_name not in adata.var_names:
            print(f"Warning: '{var_name}' not found in adata.var_names. Skipping.")
            continue

        # Extract data
        data = pd.DataFrame({
            obs_names: adata.obs[obs_names],
            var_name: adata[:, var_name].X.flatten()
        })

        # Plot boxplot
        plt.figure(figsize=(10, 6))
        boxplot = sns.boxplot(
            data=data,
            x=obs_names,
            y=var_name,
            showfliers=False,
            palette="Set3"
        )
        sns.stripplot(
            data=data,
            x=obs_names,
            y=var_name,
            color="black",
            alpha=0.4,
            jitter=True,
            size=3
        )

        # Add labels and title
        plt.title(f"Boxplot of '{var_name}' Stratified by '{obs_names}'")
        plt.xlabel(obs_names)
        plt.ylabel(var_name)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()
