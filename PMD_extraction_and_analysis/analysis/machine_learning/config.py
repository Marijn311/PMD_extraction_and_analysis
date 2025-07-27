import os

MODEL = 'lr' # See the available models in models.py
VAL_METHOD = "loo" # "tt" or "loo" or "cv"
NR_TEST_SAMPLES = 10 # Number of samples of each class to include in the test set in case VAL_METHOD == "tvt"
NR_FOLDS = 5 # For cross-validation (when chosen loo this is automatically set to the number of samples in the dataset)
NR_PC = 0 # Number of principal components. If 0, no PCA is applied

USE_PMDS = True # If True, the PMDs are used as features
USE_METADATA = False # If True, the metadata is used as features
# If both are True, the features are concatenated

# In my experiments none of the augmentations worked well.
AUGMENT_TYPE = 'none' # 'smote' or 'noise' or 'tabpfn' or 'oversample' or 'none'
# SHAP analysis is not really useful for a model with 6300+ features.
NR_SHAP = 0 # Number of top features to show in the SHAP summary plot. If 0, no SHAP summary plot is made. Plotting all features is not supported yet for 6300+ features

# You can chose to include only subsets of the PMDs
IMG_TYPES_TO_INCLUDE = ['CBF', 'CBV_LC', 'MTT', 'Tmax'] 
METRICS_TO_INCLUDE = ["mean", "median", "std", "iqr", "hart", "mean_diff", "median_diff", "std_diff", "iqr_diff", "hart_diff"]  # The model for my paper performs way better if we exlude skew and kurtosis. This is likely because these are the features with a lot of nans. (They divide by the std, which is zero in many tmax cases.)  
REGIONS_TO_INCLUDE = [i for i in range(1, 117+1)]

# You can also choose to exclude some patients due to bad quality data or other reasons
PATIENTS_TO_EXCLUDE = [ ]

# If set to True, extra debug prints and plots are shown
DEBUG_PRINTS = False

# Define data paths
DATASET_PATHS = {
    "group_001": os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "dummy_dataset", "dataset", "group_001", "pmd_dataset.nc")),
    "group_002": os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "dummy_dataset", "dataset", "group_002", "pmd_dataset.nc")),
}

NUM_CLASSES = len(DATASET_PATHS)
RANDOM_SEED = 42
