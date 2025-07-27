from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
from config import *
import shap


def load_single_pmd_dataset(data, label):
    """
    Given a single xarray dataset, load all features (PMDs) and the corresponding labels for each patient.
    Also return the names of the PMDs.
    """

    # Initiliaze storage for the pmds, labels and pmd names
    pmds = []
    labels = []
    pmd_names = []

    # Convert the xarray dataset to a xarray dataarray 
    # It is somehow very important that this line comes before the drop_sel lines. 
    # If not we get weird indexing errors.
    data = data.to_array(dim="new_dim")
        
    # Drop the unwanted patients
    for patient in data.coords['patient'].values:
        if patient in PATIENTS_TO_EXCLUDE:
            data = data.drop_sel(patient=patient)
    
    # Loop over all patients in the dataset
    patients = data.coords['patient'].values
    for patient in patients:

        # Select the data of this patient
        patient_data = data.sel(patient=patient)
        
        # Drop the unwanted PMDs
        for img_type in patient_data.coords['img_type'].values:
            if img_type not in IMG_TYPES_TO_INCLUDE:
                patient_data = patient_data.drop_sel(img_type=img_type)
        for metric in patient_data.coords['metric'].values:
            if metric not in METRICS_TO_INCLUDE:
                patient_data = patient_data.drop_sel(metric=metric)
        for region in patient_data.coords['region'].values:
            if region not in REGIONS_TO_INCLUDE:
                patient_data = patient_data.drop_sel(region=region)
        
        # Generate PMD names corresponding to the flattened array
        img_types = patient_data.coords['img_type'].values
        metrics = patient_data.coords['metric'].values
        regions = patient_data.coords['region'].values
        pmd_names_patient = []
        for metric in metrics:
            for img_type in img_types:
                for region in regions:
                    pmd_names_patient.append(f"{patient}||{img_type}_{metric}_{region}")
        
        
        # Flatten the array to shape [num_pmds]
        patient_data = patient_data.values.flatten()
        pmds.append(patient_data)
        labels.append(label)
        pmd_names.append(pmd_names_patient)
    
    return pmds, labels, pmd_names



def load_and_preprocess_pmd_dataset(datasets):
    """
    Given two xarray datasets, load all features (PMDs) and the corresponding labels for each dataset.
    Combine the PMDs and labels from both datasets into a x and y that can be used for ML.
    Then apply some simple preprocessing steps such as z-score normalization and imputing NaN values.

    Args:
        datasets (list of xarray.Dataset): list containing multiple datasets. Each entry is a xarray dataset for one patient group. e.g. stroke or epilepsy.	
  
    
    Returns:
        x (numpy.ndarray): A 2D numpy array of shape (num_samples, num_pmds) containing the PMDs for all patients.
        y (numpy.ndarray): A 1D numpy array of shape (num_samples,) containing the labels for all patients. e,g. (0: stroke, 1: epilepsy)
        pmd_names (numpy.ndarray): A 2D numpy array of shape (num_samples, num_pmds) containing the names of the PMDs.
    """

    pmds_list = []
    labels_list = []
    pmd_names_list = []
    for count, dataset in enumerate(datasets):
        pmds, labels, pmd_names = load_single_pmd_dataset(dataset, count)
        pmds_list.append(pmds)
        labels_list.append(labels)
        pmd_names_list.append(pmd_names)

        
    # # For manually checking the correspondence using the excel files
    # name = pmd_names_list[0][-1][23]
    # value = pmds_list[0][-1][23]
    # print(name)
    # print(value)

    # Combine the pmds and labels from both datasets into a single array
    x = np.array([item for sublist in pmds_list for item in sublist])
    y = np.array([item for sublist in labels_list for item in sublist])
    pmd_names = np.array([item for sublist in pmd_names_list for item in sublist])

    """
    Calculate the mean and standard deviation for every PMD.
    We ignore the NaN values when calculating the mean and std.
    Some NaN values exist because:
    -Kurtosis and skewedness are NaN when std is 0.
    -A difference (assymmetry) PMD is NaN when at least one of the two regions that form a pair is NaN. Or when their is no pair e.g. for the brainstem.
    -This might happen when Kurtosis or skewedness is nan. Or in some patients a small region might not exist at all due to bad registration. 
    """

    mean = np.nanmean(x, axis=0)
    std = np.nanstd(x, axis=0)
    
    # For some PMDs the value in the entire dataset is NaN. So the result of np.nanstd is also NaN.
    # This is the case for the difference PMDs in region 104 (brainstem) 
    # These PMDs will be removed from the dataset completely.
    indices_0_std = np.argwhere(std == 0)
    if DEBUG_PRINTS == True:
        print("\nThe following PMDs have a NaN value across the entire dataset. These PMDs will be removed:")
        for index in indices_0_std:
            index = index[0] #get the scalar value from the array
            print(pmd_names[0][index].split("||")[1])
    x = np.delete(x, indices_0_std, axis=1)
    pmd_names = np.delete(pmd_names, indices_0_std, axis=1)
    std = np.delete(std, indices_0_std)
    mean = np.delete(mean, indices_0_std)

    # Apply z-score normalization to normalize the PMDs
    x = (x - mean) / std

    if DEBUG_PRINTS == True:
        plt.hist(x.flatten(), bins=100)
        plt.title("Histogram of the PMDs after normalization (all PMDs combined in same histogram)")
        plt.show()
    
    # Count NaN values per class
    for class_label in np.unique(y):
        class_nan_count = np.sum(np.isnan(x[y == class_label]))
        total_class_count = np.sum(y == class_label) * x.shape[1]
        if DEBUG_PRINTS == True:
            print(f"Number of NaN values across all pmds and all patients in class {class_label}: {class_nan_count}")
            print(f"Total number of pmd times number of patients in class {class_label}: {total_class_count}")
            print(f"Percentage of NaN values in class {class_label}: {class_nan_count / total_class_count * 100:.2f}%\n")
        if class_nan_count / total_class_count * 100 > 3:
            raise ValueError(f"Warning: More than 3% of the data in class {class_label} is NaN. You may need to check the data quality.")
      
    # Use median imputation to fill the NaN values in the dataset
    indices = np.argwhere(np.isnan(x))
    for ind in indices:
        x[ind[0], ind[1]] = np.nanmedian(x[:, ind[1]])
         
        if DEBUG_PRINTS == True:
            print(f"NaN values in the data:")
            print(f"{pmd_names[ind[0]][ind[1]]}\n")

    return x, y, pmd_names


def pca(x_train, x_test, n_components):
    """
    Perform PCA on the training and test datasets to reduce the number of features (PMDs)
    by computing the principal components and projecting the data onto them.
    Args:
        x_train (numpy.ndarray): The training data
        x_test (numpy.ndarray): The test data
        n_components (int): Number of principal components to keep.
    
    Returns:
        tuple: Transformed training and test pmd matrices with reduced dimensions.
    """
    pca = PCA(n_components=n_components)
    x_train_pca = pca.fit_transform(x_train)
    x_test_pca = pca.transform(x_test)
    if DEBUG_PRINTS == True:
        print(f"\nPCA applied with {n_components} components.")
        print(f"Original shape: {x_train.shape}")
        print(f"Transformed shape: {x_train_pca.shape}")
        print(f"Explained variance ratio of the components: {pca.explained_variance_ratio_}")
        print(f"Total explained variance: {np.sum(pca.explained_variance_ratio_)}")
    return x_train_pca, x_test_pca


def shap_analysis(explainer, all_shap_values, pmd_names):	 
    """
    Perform SHAP analysis on the model predictions and plot the results.
    Args:
        explainer (shap.Explainer): The SHAP explainer object for the model.
        all_shap_values (numpy.ndarray): The SHAP values for all samples and PMDs.
        pmd_names (list): The names of the PMDs.

    (SHAP analysis is not really usefull in a model with 6300+ features.)
    """

    if NR_PC > 0:
        raise ValueError("SHAP analysis with PCA features is uninterpretable.")

    if NR_SHAP > 50:
        raise ValueError("SHAP analysis with more than 50 features results in unreadable plots. Please set NR_SHAP to a lower value.")


    # Keep only the PMD names for the first sample, since we dont care about the patient name.
    pmd_names = pmd_names[0,:] # This is a list of lists, so we need to flatten it.

    # In case of CV, the shap values are passed as nr_of_folds x nr_of_sample_per_fold x nr_of_pmds. So reshaping is needed.
    if all_shap_values.ndim == 3:
        all_shap_values = all_shap_values.reshape(all_shap_values.shape[0] * all_shap_values.shape[1], all_shap_values.shape[2])

    # Split every string in pmd_names on the || and keep only the second element, since we don't care about the patient name.
    pmd_names = [x.split("||")[1] for x in pmd_names]   

    shap.summary_plot(all_shap_values, max_display=NR_SHAP, feature_names=pmd_names)
    shap.summary_plot(all_shap_values, max_display=NR_SHAP, plot_type='bar', feature_names=pmd_names) 
    # decision plot don't work for 6300+ features, and with PCA features the info is not that useful.
    # shap.decision_plot(base_value=explainer.expected_value, shap_values=all_shap_values, feature_names=pmd_names, feature_display_range=slice(-1,-NR_SHAP,-1), ignore_warnings=True)

    # Calculate which percentage of the total absolute shap values comes from the top n pmds.
    # Get the indices of the top n pmds
    top_n_indices = np.argsort(np.abs(all_shap_values).mean(axis=0))[::-1][:NR_SHAP]
    total_sum = np.abs(all_shap_values).sum()
    top_n_sum = np.abs(all_shap_values[:, top_n_indices]).sum()
    top_n_percentage = top_n_sum / total_sum * 100
    print(f"\nTop {NR_SHAP} pmds contribute {top_n_percentage:.2f}% to the total absolute SHAP values.")

    # Get the pmd names of the top n pmds
    top_n_pmd_names = [pmd_names[i] for i in top_n_indices]
    print(f"\nTop {NR_SHAP} pmds: {top_n_pmd_names}")


    print("\nWe have assigned each PMD a rank based on SHAP value they have.")
    print("The least important PMD is ranked 1, the most important PMD gets the highest rank.")

    pmds_ranked = np.argsort(np.abs(all_shap_values).mean(axis=0))[::-1]
    ranked_pmd_names = [pmd_names[i] for i in pmds_ranked]
    
    # Do some formatting reasons.
    ranked_pmd_names = [x.replace("_LC_", "_") for x in ranked_pmd_names]
    ranked_pmd_names = [x.replace("_diff", "diff") for x in ranked_pmd_names]
    
    
    # Create a 2D array with PMD information and rankings
    pmd_ranking = []
    for i, pmd_idx in enumerate(pmds_ranked):
        pmd_name = ranked_pmd_names[i]
        pmd_parts = pmd_name.split("_")
        img_type = pmd_parts[0]
        metric = pmd_parts[1]
        region = pmd_parts[2]
        rank = len(pmds_ranked) - i  # Highest rank for most important PMD
        pmd_ranking.append([pmd_name, img_type, metric, region, rank])
    
    pmd_ranking = np.array(pmd_ranking)

    # print the avrage rank of all pmds per img_type and metric, and region
    unique_img_types = np.unique(pmd_ranking[:, 1])
    unique_metrics = np.unique(pmd_ranking[:, 2])   
    unique_regions = np.unique(pmd_ranking[:, 3])
    # Calculate average ranks for image types and sort by rank
    img_type_ranks = []
    for img_type in unique_img_types:
        indices = np.where(pmd_ranking[:, 1] == img_type)[0]
        if len(indices) > 0:
            avg_rank = np.mean(pmd_ranking[indices, 4].astype(float))
            img_type_ranks.append((img_type, avg_rank))
    
    # Sort image types by average rank (highest rank first)
    img_type_ranks.sort(key=lambda x: x[1], reverse=True)
    
    # Normalize image type ranks
    if img_type_ranks:
        max_img_rank = max(rank for _, rank in img_type_ranks)
        min_img_rank = min(rank for _, rank in img_type_ranks)
        rank_range = max_img_rank - min_img_rank
        
        for img_type, avg_rank in img_type_ranks:
            if rank_range > 0:
                normalized_rank = (avg_rank - min_img_rank) / rank_range
            else:
                normalized_rank = 0.5  # All ranks are the same
            print(f"The average rank of all PMDs with Image Type: {img_type} -> {avg_rank:.2f} (normalized: {normalized_rank:.3f})")
    
    # Calculate average ranks for metrics and sort by rank
    metric_ranks = []
    for metric in unique_metrics:
        indices = np.where(pmd_ranking[:, 2] == metric)[0]
        if len(indices) > 0:
            avg_rank = np.mean(pmd_ranking[indices, 4].astype(float))
            metric_ranks.append((metric, avg_rank))
    
    # Sort metrics by average rank (highest rank first)
    metric_ranks.sort(key=lambda x: x[1], reverse=True)
    
    # Normalize metric ranks
    if metric_ranks:
        max_metric_rank = max(rank for _, rank in metric_ranks)
        min_metric_rank = min(rank for _, rank in metric_ranks)
        rank_range = max_metric_rank - min_metric_rank
        
        for metric, avg_rank in metric_ranks:
            if rank_range > 0:
                normalized_rank = (avg_rank - min_metric_rank) / rank_range
            else:
                normalized_rank = 0.5  # All ranks are the same
            print(f"The average rank of all PMDs with Metric: {metric} -> {avg_rank:.2f} (normalized: {normalized_rank:.3f})")
    
    # Calculate average ranks for regions and sort by rank
    region_ranks = []
    for region in unique_regions:
        indices = np.where(pmd_ranking[:, 3] == region)[0]
        if len(indices) > 0:
            avg_rank = np.mean(pmd_ranking[indices, 4].astype(float))
            region_ranks.append((region, avg_rank))
    
    # Sort regions by average rank (highest rank first)
    region_ranks.sort(key=lambda x: x[1], reverse=True)
    
    # Normalize region ranks
    if region_ranks:
        max_region_rank = max(rank for _, rank in region_ranks)
        min_region_rank = min(rank for _, rank in region_ranks)
        rank_range = max_region_rank - min_region_rank
        
        for region, avg_rank in region_ranks:
            if rank_range > 0:
                normalized_rank = (avg_rank - min_region_rank) / rank_range
            else:
                normalized_rank = 0.5  # All ranks are the same
            print(f"The average rank of all PMDs with Region: {region} -> {avg_rank:.2f} (normalized: {normalized_rank:.3f})")




def load_single_metadataset(metadata, label):
    """
    Given a single pandas dataframe, load all features (metadata) and the corresponding labels for each patient.
    Also return the names of the metadata features.
    """

    # metadata should a np.ndarray of shape (num_samples, num_features)
    # labels should be a np.ndarray of shape (num_samples,)
    # metadata_names should be an np.ndarray of shape (nr_samples, num_features)
    
    # Drop the unwanted patients
    metadata = metadata[~metadata["patient_id"].isin(PATIENTS_TO_EXCLUDE)]
            

    metadata = metadata.drop(columns=["patient_id"])
    metadata_names = metadata.columns.values.tolist()
    #convert metadata_names to array of shape (num_samples, num_features) (the feature names are the same for each sample, but we expect this format)
    metadata_names = np.array([metadata_names] * metadata.shape[0])
    metadata = metadata.to_numpy()
    labels= np.full(metadata.shape[0], label)

    return metadata, labels, metadata_names



def load_and_preprocess_metadata_dataset(metadatasets):
    """ 
    This file should take the pandas dataframes with the metadata per group and load them as a x and y. which can be used for ML
    such as logistic regression.
    """

    metadata_list = []
    labels_list = []
    metadata_names_list = []
    for count, dataset in enumerate(metadatasets):
        metadata, labels, metadata_names = load_single_metadataset(dataset, count)
        metadata_list.append(metadata)
        labels_list.append(labels)
        metadata_names_list.append(metadata_names)
    
    # Combine the metadata and labels from both datasets into a single array
    x = np.array([item for sublist in metadata_list for item in sublist])
    y = np.array([item for sublist in labels_list for item in sublist])
    metadata_names = np.array([item for sublist in metadata_names_list for item in sublist])

    """
    Calculate the mean and standard deviation for every metadata.
    We ignore the NaN values when calculating the mean and std.
    """

    mean = np.nanmean(x, axis=0)
    std = np.nanstd(x, axis=0)
    
    # Apply z-score normalization to normalize the metadata
    x = (x - mean) / std

    if DEBUG_PRINTS == True:
        plt.hist(x.flatten(), bins=100)
        plt.title("Histogram of the metadata after normalization (all metadata combined in same histogram)")
        plt.show()
    
    # Count NaN values per class
    for class_label in np.unique(y):
        class_nan_count = np.sum(np.isnan(x[y == class_label]))
        total_class_count = np.sum(y == class_label) * x.shape[1]
        if DEBUG_PRINTS == True:
            print(f"Number of NaN values across all metadata and all patients in class {class_label}: {class_nan_count}")
            print(f"Total number of metadata times number of patients in class {class_label}: {total_class_count}")
            print(f"Percentage of NaN values in class {class_label}: {class_nan_count / total_class_count * 100:.2f}%\n")
        if class_nan_count / total_class_count * 100 > 3:
            raise ValueError(f"Warning: More than 3% of the data in class {class_label} is NaN. You may need to check the data quality.")
      
    # Use median imputation to fill the NaN values in the dataset
    indices = np.argwhere(np.isnan(x))
    for ind in indices:
        x[ind[0], ind[1]] = np.nanmedian(x[:, ind[1]])
         
        if DEBUG_PRINTS == True:
            print(f"NaN values in the data:")
            print(f"{metadata_names[ind[0]][ind[1]]}\n")


    a=4
    return x, y, metadata_names