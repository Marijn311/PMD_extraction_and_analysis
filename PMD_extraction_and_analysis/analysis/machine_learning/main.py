import xarray as xr
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold	
import random
from metrics import *
from utils import *
from config import *	
from models import *
from augment import *


# Print some info of the configuration
print(" ")
print("Model: ", MODEL)
print("Validation method: ", VAL_METHOD)
print("class names: ", DATASET_PATHS.keys())
print("Included image types: ", IMG_TYPES_TO_INCLUDE)
print("Included metrics: ", METRICS_TO_INCLUDE)
print("Included regions: ", REGIONS_TO_INCLUDE)
print(" ")

# Load the data
if USE_PMDS == True:
    datasets = []
    for key, value in DATASET_PATHS.items():
        dataset = xr.open_dataset(value)
        datasets.append(dataset)
    x, y, feature_names = load_and_preprocess_pmd_dataset(datasets)

if USE_METADATA == True:
    metadatasets = []
    for key, value in DATASET_PATHS.items():
        value = value.replace("pmds.nc", "metadata.xlsx") 
        metadataset = pd.read_excel(value)
        metadatasets.append(metadataset)
    x_meta, y_meta, feature_names_meta = load_and_preprocess_metadata_dataset(metadatasets)
    

if USE_PMDS == True and USE_METADATA == True:
    # Check if the number of samples is the same
    if x.shape[0] != x_meta.shape[0]:
        raise ValueError("The number of samples in the PMD dataset and metadata dataset do not match.")
    # Concatenate the PMD and metadata 
    x = np.concatenate((x, x_meta), axis=1)
    feature_names = np.concatenate((feature_names, feature_names_meta), axis=1)

if USE_PMDS == False and USE_METADATA == True:
    x = x_meta
    y = y_meta
    feature_names = feature_names_meta

# Print the distribution of the classes
unique, counts = np.unique(y, return_counts=True)
for label, count in zip(unique, counts):
    print(f"Class {label}: {count} samples ({count/len(y)*100:.2f}%)")

# Initialize metrics storage
test_scores = {'test_true': [], 'test_prob': []}

# Initialize SHAP values storage
all_shap_values = []

# Define a model
model = select_model(y)

# Determine the number of folds 
# For tt, the number of folds is set to 1
# For loo, the number of folds is set to the number of train samples
# For cv, the number of folds is set to NR_FOLDS in config.py
if VAL_METHOD == "tt": 
    NR_FOLDS = 1
    kf = [(None, None)]  # Placeholder for manual split
if VAL_METHOD == "cv" or VAL_METHOD == "loo":
    if VAL_METHOD == "loo":
        NR_FOLDS = len(x)
    kf = KFold(n_splits=NR_FOLDS, shuffle=True, random_state=RANDOM_SEED).split(x)

for fold, (train_indices, test_indices) in enumerate(kf):
    print(f"\nFold: {fold + 1}")

    if VAL_METHOD == "tt":
        # Take NR_TEST_SAMPLES samples from each class for the test set
        x_test = np.zeros((NR_TEST_SAMPLES * NUM_CLASSES, x.shape[1]))
        y_test = np.zeros(NR_TEST_SAMPLES * NUM_CLASSES)
        for i in range(NUM_CLASSES):  # Loop over the classes
            indices = np.where(y == i)[0]  # Get the indices of the samples of this class
            indices = random.sample(list(indices), NR_TEST_SAMPLES)  # Randomly take NR_TEST_SAMPLES
            x_test[i * NR_TEST_SAMPLES:(i + 1) * NR_TEST_SAMPLES] = x[indices]  # Add to test set
            y_test[i * NR_TEST_SAMPLES:(i + 1) * NR_TEST_SAMPLES] = y[indices]  # Add to test labels
            x_train = np.delete(x, indices, axis=0)  # Remove from training set
            y_train = np.delete(y, indices)
    else:
        # Split the data into train and test sets for this fold
        x_train, x_test = x[train_indices], x[test_indices]
        y_train, y_test = y[train_indices], y[test_indices]

    # Do PCA if configured
    if NR_PC > 0:
        x_train, x_test = pca(x_train, x_test, NR_PC)
    
    # Augment the train data if configured
    if AUGMENT_TYPE == 'noise':
        x_train, y_train = aug_noise(x_train, y_train)
    if AUGMENT_TYPE == 'smote':
        x_train, y_train = aug_smote(x_train, y_train)
    if AUGMENT_TYPE == 'tabpfn':
        x_train, y_train = aug_tabpfn(x_train, y_train)

    # Train the model
    model.fit(x_train, y_train)

    # Make predictions
    test_pred = model.predict(x_test)
    test_prob = model.predict_proba(x_test)

    # Add to results to dictionary
    test_scores['test_prob'].append(test_prob)
    test_scores['test_true'].append(y_test)

    if DEBUG_PRINTS == True:
        print(f"True labels: {y_test}")
        print(f"Predicted labels: {test_pred}")
        print(f"Predicted probabilities: {test_prob}")

    # Generate and store the SHAP values if configured
    if NR_SHAP > 0:
        explainer = shap.Explainer(model, x_train)
        shap_values = explainer(x_test)
        all_shap_values.append(shap_values.values)

# Generate SHAP summary plot if configured
if NR_SHAP > 0:
    all_shap_values = np.array([shap_values for shap_values in all_shap_values])  # Unpack the list
    all_shap_values = np.squeeze(all_shap_values)  # Remove singleton dimension
    shap_analysis(explainer, all_shap_values, feature_names)

# Aggregate metrics across all folds
test_true = np.concatenate(test_scores['test_true'])
test_prob = np.concatenate(test_scores['test_prob'])

if NUM_CLASSES == 2:
    # Generate metrics with confidence intervals
    metrics_with_ci(test_true, test_prob)
else:
    # Generate metrics for multiclass problems
    metrics_multiclass(test_true, test_prob)

print("\nDone!")
