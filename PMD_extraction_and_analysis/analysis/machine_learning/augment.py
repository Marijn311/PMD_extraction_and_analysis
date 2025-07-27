# Augmentation did not prove useful in my dataset/experiments.

# import torch
# from tabpfn_extensions import TabPFNClassifier, TabPFNRegressor, unsupervised
# from imblearn.over_sampling import SMOTE
# import numpy as np
# from config import *


def aug_smote(x_train, y_train):
    """
    SMOTE (Synthetic Minority Over-sampling Technique) is an oversampling technique that generates synthetic samples 
    for the minority class by interpolating between existing samples.

    -SMOTE automatically finds the minority class by looking at the given lables.
    
    -SMOTE generate new samples for each feature independently.
    
    -A nearest neighbor in SMOTE is determined based on the overall feature space, not on individual features separately. 
    The algorithm typically uses Euclidean distance (or another distance metric) to find the k-nearest neighbors in multi-dimensional space.
    
    -SMOTE uses the corresponding feature from these fixed k-neighbours to interpolate between the samples of the minority class.
    
    -The formula for the interpolation is as follows:
    new_sample = x_i + (x_j - x_i) * random_number
    where:
    - x_i is the feature vector of the original sample
    - x_j is the feature vector of the nearest neighbor
    - random_number is a random number between 0 and 1
    This is done for each nearest neighbor X_j
    """
    
    # Implement smote
    smote = SMOTE(sampling_strategy='minority', k_neighbors=20, random_state=RANDOM_SEED)	#i think more neighbours is better 
    x, y = smote.fit_resample(x_train, y_train)
    # The output x,y of SMOTE are already a mix of the original and the synthetic samples.
    # The number of samples that are generated is equal to the difference between the number of samples in the majority and minority class.
    # So the output x,y of SMOTE are perfectly balanced. 

    return x, y



def aug_noise(x_train, y_train):
    """
    This function augments the data by adding noise to the minority class samples and saving these as new samples.
    The noise is generated from a normal distribution with mean 0 and std 0.1.
    The number of samples to generate is equal to the difference between the number of samples in the majority and minority class.
    The output x,y are perfectly balanced.
    """

    # Calculate the number of samples that need to be added to each the minority class
    unique, counts = np.unique(y_train, return_counts=True)
    class_counts = dict(zip(unique, counts))
    minority_class_id = min(class_counts, key=class_counts.get)
    majority_class_id = max(class_counts, key=class_counts.get)
    
    # You can either manually set the number of samples to generate, or let the function calculate the number of samples you need to balance the classes.
    num_samples_to_gen = class_counts[majority_class_id] - class_counts[minority_class_id]
    # num_samples_to_gen = 25

    # Get the indices of the samples of the minority class
    indices = np.where(y_train == minority_class_id)[0]

    # Randomly select num_samples_to_gen samples fromm the minority class.
    # We do this with replacement since we want to generate more samples than we have.
    indices_to_add_noise_to = np.random.choice(indices, num_samples_to_gen, replace=True)

    # Make an empty array for the augmented samples
    x_noise = np.zeros((indices_to_add_noise_to.shape[0], x_train.shape[1]))
    y_noise = np.zeros(indices_to_add_noise_to.shape[0])

    # Add noise to the samples.
    # The noise is generated from a normal distribution with mean 0 and std 0.1.
    for count, index in enumerate(indices_to_add_noise_to):
        x_noise[count] = x_train[index] + np.random.normal(0, 0.1, x_train[index].shape)
        y_noise[count] = y_train[index]

    # Add the augmented samples to the original dataset
    x = np.concatenate((x_train, x_noise), axis=0)
    y = np.concatenate((y_train, y_noise), axis=0)

    return x, y


def aug_tabpfn(x_train, y_train): 
    
    # Initialize unsupervised model
    model_unsupervised = unsupervised.TabPFNUnsupervisedModel(
        tabpfn_clf=TabPFNClassifier(n_estimators=1),
        tabpfn_reg=TabPFNRegressor(n_estimators=1),
    )

    # Convert data to torch tensors
    x_tensor = x_train.clone().detach().to(dtype=torch.float32) if isinstance(x_train, torch.Tensor) else torch.tensor(x_train, dtype=torch.float32)
    y_tensor = y_train.clone().detach().to(dtype=torch.float32) if isinstance(y_train, torch.Tensor) else torch.tensor(y_train, dtype=torch.float32)

    # Take only the train data for the class 1 for which we want to generate synthetic data
    x_tensor_class1 = x_tensor[y_tensor == 1]

    # This model generates new samples feature by feature.
    # Generating 6300 features would take impossibly long. 
    # I have to run this in combination with PCA to make it somewhat manageable.
    # But with PCA = 130 it takes about 30 min per fold. So 150 folds would take 37.5 hours.
    # This is not OK. So i will test it in normal mode first.
    # Or try this on the HPC cluster. On laptop it takes too long.
    assert x_tensor_class1.shape[1] < 200, "The number of features is too high. Please use PCA."

    # Fit the unsupervised model, to the class 1 data
    model_unsupervised.fit(x_tensor_class1)

    # Generate synthetic data
    # This function seems to get slower for every steps (feature). so features 1-5 are very quick, 5-10 are slow, and 10-25 are very slow.
    x_synthetic = model_unsupervised.generate_synthetic_data(n_samples=92)
    y_synthetic = torch.ones(x_synthetic.shape[0], dtype=torch.float32)

    print(f"Number of class 1 samples: {x_tensor_class1.shape[0]}")
    print(f"Number of synthetic class 1 samples generated: {x_synthetic.shape[0]}")

    # Concat the original data with the synthetic data
    x = torch.cat((x_tensor, x_synthetic), dim=0)
    y = torch.cat((y_tensor, y_synthetic), dim=0)

    return x, y

def aug_oversample(x_train, y_train):
    """Simply determine the difference between the number of samples in the majority and minority class.
    Then randomly duplicate samples from the minority class until the classes are balanced.
    """
    # Calculate the number of samples that need to be added to each the minority class
    unique, counts = np.unique(y_train, return_counts=True)
    class_counts = dict(zip(unique, counts))
    minority_class_id = min(class_counts, key=class_counts.get)
    majority_class_id = max(class_counts, key=class_counts.get)

    # You can either manually set the number of samples to generate, or let the function calculate the number of samples you need to balance the classes.
    num_samples_to_gen = class_counts[majority_class_id] - class_counts[minority_class_id]

    # Get the indices of the samples of the minority class
    indices = np.where(y_train == minority_class_id)[0]

    # Randomly select num_samples_to_gen samples fromm the minority class.
    # We do this with replacement since we want to generate more samples than we have.
    indices_to_copy = np.random.choice(indices, num_samples_to_gen, replace=True)

    # Make an empty array for the augmented samples
    x_copied = np.zeros((indices_to_copy.shape[0], x_train.shape[1]))
    y_copied = np.zeros(indices_to_copy.shape[0])

    # create copied samples
    for count, index in enumerate(indices_to_copy):
        x_copied[count] = x_train[index]
        y_copied[count] = y_train[index]

    # Add the augmented samples to the original dataset
    x = np.concatenate((x_train, x_copied), axis=0)
    y = np.concatenate((y_train, y_copied), axis=0)

    return x, y

