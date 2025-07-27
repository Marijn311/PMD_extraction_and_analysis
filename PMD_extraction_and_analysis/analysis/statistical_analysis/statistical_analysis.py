import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import itertools
import seaborn as sns
from scipy.stats import ranksums


"""
This script compares the distributions of PMDs between two or more patient groups.
For each PMD, the difference in distribution between groups is assessed using the Wilcoxon rank-sum test and Cohen's d effect size.
To account for multiple comparisons across 113 cerebral regions, Bonferroni correction is applied.
PMDs are whose adjusted p-values are above a specified threshold (P_THRES) or whose Cohen's d values are below a specified threshold (D_THRES) are excluded dropped.
The remaining PMDs which have statistically significant difference between groups are included in the final results table.
Additionally, the script provides the option to visualize PMD distributions using boxplots.

NOTE: Only the first two datasets in the DATASETS dictionary are used for to calculate the P and D value. 
Hence, the results table which shows the relevant PMDs with P-values and D-values is a comparison between the first two datasets.
There is an option to add more datasets to the DATASETS dictionary, such that you can visualize the PMD distributions of more than two datasets in the boxplot.
If you want the relevant PMDs for other datasets comparisons, you can change the order of the datasets in the DATASETS dictionary.
"""

LABELS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "dummy_dataset", "atlas", "HO_labels.xlsx"))
DATASET_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "dummy_dataset", "dataset"))
DATASETS = {
    'group_001': 'group_001/pmd_dataset.nc',
    'group_002': 'group_002/pmd_dataset.nc',
}

P_THRES = 0.05
D_THRES = 0.3 

PATIENTS_TO_EXCLUDE = []

# Analyse PMD distribution
TYPE = 'distribution'
BONF_CORR_NUM = 113  
ALL_METRICS = ['mean', 'median', 'std', 'iqr', 'skew', 'kurt', 'hart'] # You might not want to include every metric or img_type
ALL_IMG_TYPES = ['CBF', 'CBV_LC', 'MTT', 'Tmax'] 
REGIONS_TO_EXCLUDE = [] 

# # Analyse PMD assymmetry 
# TYPE = 'asymmetry'
# BONF_CORR_NUM = 56  
# ALL_METRICS = ['mean_diff', 'median_diff', 'std_diff', 'iqr_diff', 'skew_diff', 'kurt_diff', 'hart_diff'] 
# ALL_IMG_TYPES = ['CBF', 'CBV_LC', 'MTT', 'Tmax'] 
# REGIONS_TO_EXCLUDE = [104] # Exclude lateral ventricles and the brainstem. the brainstem has no counterpart in the other hemisphere.

SHOW_PLOT = False
if SHOW_PLOT == True:
    ALL_METRICS = ['median'] 
    ALL_IMG_TYPES = ['CBF'] 


def wilcoxon_rank_sum_test(group1, group2, characteristic):
    group1_data = group1[characteristic]
    group2_data = group2[characteristic]
    _, p = ranksums(group1_data, group2_data)
    return p


def load_dataset(dataset_name):
   
    if dataset_name not in DATASETS:
        raise ValueError(f"Dataset {dataset_name} not found. Available datasets: {list(DATASETS.keys())}")

    # Load the dataset    
    data_path = os.path.join(DATASET_PATH, DATASETS[dataset_name])
    data = xr.open_dataarray(data_path)

    # Keep only the selected image type and metric
    data = data.sel(metric=METRIC).drop_vars('metric')
    data = data.sel(img_type=IMG_TYPE).drop_vars('img_type')

    # Set the attribute called name to dataset_name
    data.attrs['name'] = dataset_name

    # Exclude patients in PATIENTS_TO_EXCLUDE
    data = data.where(~data.patient.isin(PATIENTS_TO_EXCLUDE), drop=True)
    return data

def add_line_breaks(label, max_words_per_line= 2):
    """Add line breaks after every max_words_per_line words in the label. This is to make the plotting nicer."""
    words = label.split()
    return '\n'.join([' '.join(words[i:i + max_words_per_line]) for i in range(0, len(words), max_words_per_line)])


def plot_comparison(datasets, parameter):
    """Create box plots for the specified parameter between two groups, ordered by p-values."""
    
    # Load the labels (the names that correspond to the brain region numbers) 
    labels = pd.read_excel(LABELS_PATH)
    labels = labels[['id_remapped', 'name']]

    dfs = []
    for dataset in datasets:
        name = dataset.attrs['name']
        dataset = dataset.to_pandas()
        dataset.attrs['name'] = name # The xarray does not copy the name attribute to the pandas dataframe, so this is a little workaround
        df = pd.DataFrame(dataset).melt()
        df['Group'] = dataset.attrs['name'] # To distinguish between the patient groups in the boxplot
        dfs.append(df)
        # Reshaping the data frames using melt:
        # This transforms the data into a long format with only two columns: 'region' and 'value'.
        # Each row represents a single data point, with 'region' indicating the brain region and 'value' being the actual measurement.
        # Note that the patient ID is not retained in this transformation.
        # This format is suitable for creating boxplots, where 'region' serves as the categorical x-axis and 'value' as the numerical y-axis.

    df = pd.concat([df for df in dfs])
    
    # Calculate p-values and cohen's d-value for each region
    df['p'] = np.nan
    df['d'] = np.nan

    for region in df['region'].unique():
        values1 = datasets[0].sel(region=region).values
        values2 = datasets[1].sel(region=region).values
        # Remove the nan values from the arrays before calculating the statistics
        values1 = values1[~np.isnan(values1)]
        values2 = values2[~np.isnan(values2)]

        if len(values1) != 0 and len(values2) != 0:
            
            # Calculate P-value using Wilcoxon rank-sum test
            p_value = wilcoxon_rank_sum_test(
                pd.DataFrame({parameter: values1}),
                pd.DataFrame({parameter: values2}),
                parameter
            )
            
            # Bonferroni correction for multiple comparisons
            p_value = p_value * BONF_CORR_NUM  

            # Calculate Cohen's d
            """Cohen's d is designed for comparing two groups. 
            It takes the difference between two means and expresses it in standard deviation . 
            It tells you how many standard deviations lie between the two means.
            Cohen's d is a measure of effect size. AKA if the difference is significant, how big is the difference.
            This tells us more about how useful the difference is in practice to differentiate between groups.
            """
            mean1 = np.mean(values1)
            mean2 = np.mean(values2)
            std1 = np.std(values1)
            std2 = np.std(values2)
            n1 = len(values1)
            n2 = len(values2)
            # Use the pooled standard deviation
            pooled_std = np.sqrt(((n1 - 1) * std1 ** 2 + (n2 - 1) * std2 ** 2) / (n1 + n2 - 2))
            d_value = (mean1 - mean2) / pooled_std

            # Set the statistics in the df
            df.loc[df['region'] == region, 'p'] = p_value
            df.loc[df['region'] == region, 'd'] = d_value
        else:
            df.loc[df['region'] == region, 'p'] = np.nan
            df.loc[df['region'] == region, 'd'] = np.nan
    
    # Remove the regions that have nan p-values. 
    # If there is a nan p-value, then there is not a single datapoint for this region in one of the groups. 
    regions_with_nan_pvalues = df[df['p'].isnull()]['region'].unique()
    df = df[~df['region'].isin(regions_with_nan_pvalues)]
    df = df[~df['region'].isin(REGIONS_TO_EXCLUDE)] # This is also a good time to remove the regions that were set by the user

    # Get unique regions ordered by their p-values, we will use this to order the boxplots based on the p-values and to set the x-axis labels
    ordered_regions = df.groupby('region')['p'].first().sort_values().index

    # Set the max regions to the number of regions that have a p-value smaller than p threshold and a d-value larger than d threshold
    max_regions = len(df[(df['p'] < P_THRES) & (df['d'].abs() > D_THRES)]['region'].unique())
    print(f"For {METRIC} {IMG_TYPE} there are {max_regions} regions with statistically significant differences between the groups.")
  
    # Filter DataFrame to only include regions up to max_regions
    ordered_regions = ordered_regions[:max_regions]
    df = df[df['region'].isin(ordered_regions)].reset_index(drop=True)

    if SHOW_PLOT == True:
        # Set style and create figure
        plt.style.use('default')
        sns.set_style("ticks")
        _, ax = plt.subplots(figsize=(30 / 2.54, 15 / 2.54)) 
        colors = ['#ff7f0e', '#1f77b4','#2ca02c']  # orange and blue and green
        font_size = 14

        # Create box plot using region for x-axis positioning
        sns.boxplot(x='region', y='value', hue='Group', data=df, 
                    palette=colors,
                    showfliers=False,  # this shows the outliers with an extra dark circle, but all values are always included
                    boxprops={'facecolor': 'none'},
                    order=ordered_regions)  # Specify the order explicitly

        # Add individual points
        sns.stripplot(x='region', y='value', hue='Group', data=df,
                    dodge=True, size=4, alpha=0.4, jitter=0.2, 
                    palette=colors,
                    order=ordered_regions)  # Specify the order explicitly

        # Set custom x-axis with region numbers, region names, p-values, and d-values
        new_labels = []
        for r in ordered_regions:
            label_row = labels.loc[labels["id_remapped"] == r, "name"]
            if not label_row.empty:
                region_name = label_row.values[0]
            else:
                region_name = "Unknown"
            region_name = add_line_breaks(region_name)
            p_value = df[df["region"] == r]["p"].iloc[0]
            d_value = df[df["region"] == r]["d"].iloc[0]
            new_labels.append(f'{region_name}\n(p={p_value:.3f}\nd={d_value:.2f})')

        ax.set_xticks(range(len(ordered_regions)))
        ax.set_xticklabels(new_labels, rotation=0, ha='center')

        # Set the x-axis limits to show only the plotted regions
        ax.set_xlim(-0.5, max_regions - 0.5)
        # Set the y-axis limits to ymax plus 1 to make room for the legend
        y_max = df['value'].max() + 1.5
        y_min = df['value'].min() 
        ax.set_ylim(y_min, y_max)

        # Customize appearance
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel(f'{METRIC.capitalize()} r{IMG_TYPE.upper()}', fontsize=font_size)
        ax.tick_params(axis='both', which='major', labelsize=font_size)
        ax.set_xlabel("")  # Remove the general x-axis label

        boxplot_handles, boxplot_labels = ax.get_legend_handles_labels()
        boxplot_handles = boxplot_handles[:len(datasets)]  # Only take the first set of handles (boxplot)
        boxplot_labels = boxplot_labels[:len(datasets)]  # Only take the first set of labels (boxplot)
        ax.legend(boxplot_handles, boxplot_labels, loc='upper left', fontsize=font_size, frameon=False)
        plt.tight_layout()
        plt.savefig(f'PMD_extraction_and_analysis/analysis/statistical_analysis/{METRIC}_{parameter}_boxplot.png', dpi=600, bbox_inches='tight')
        plt.show()

    # Get the regions, p-values, and d-values for the regions that are plotted
    plotted_regions = []
    plotted_ps = []
    plotted_ds = []
    
    # Get the p-values and d-values associated with the plotted regions
    for r in df['region'].unique():
        label_row = labels.loc[labels["id_remapped"] == r, "name"]
        region_name = label_row.values[0]
        # Get the p-value and d-value for the region
        p_value = df[df["region"] == r]["p"].iloc[0]
        d_value = df[df["region"] == r]["d"].iloc[0]
        p_value = f"{p_value:.3f}"
        d_value = f"{d_value:.3f}"
        plotted_regions.append(region_name)
        plotted_ps.append(p_value)
        plotted_ds.append(d_value)

    return plotted_regions, plotted_ps, plotted_ds

  
def analyze_datasets(datasets_to_analyse: list, parameter: str):
    """Analyze and compare multiple datasets."""
    datasets = []
    for dataset_name in datasets_to_analyse:
        if dataset_name is None:
            continue
        dataset = load_dataset(dataset_name)
        datasets.append(dataset)
    regions,p,d = plot_comparison(datasets, parameter)
    return regions,p,d




if __name__ == "__main__":

    assert len(DATASETS) >= 2, "There should be at least 2 keys in DATASETS"
    # Define the datasets to compare, set undefined datasets to None
    datasets_to_analyse = list(DATASETS.keys())

    # Create a DataFrame to store results with metrics as rows and image types as columns
    results = pd.DataFrame(index=ALL_METRICS, columns=ALL_IMG_TYPES)  

    # Iterate over all combinations of image type and metrics
    for img_type, metric in itertools.product(ALL_IMG_TYPES, ALL_METRICS):
        global IMG_TYPE, METRIC
        IMG_TYPE = img_type
        METRIC = metric

        regions,p,d = analyze_datasets(datasets_to_analyse, IMG_TYPE)
        
        # Format the results with line breaks between regions and their corresponding p and d values
        formatted_results = "\n".join([
            f"{region} (p:{p_val}, d:{d_val})\n"
            for region, p_val, d_val in zip(regions, p, d)
        ])
        results.loc[metric, img_type] = formatted_results

    output_filename = f"PMD_extraction_and_analysis/analysis/statistical_analysis/PMD_statistical_analysis_results_{TYPE}.xlsx"
    results.to_excel(output_filename, index=True, header=True)

    # Save the results to an Excel file
    print("Done!")