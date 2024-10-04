import numpy as np

def calculate_alpha(data):
    """
    Calculate alpha values for the Dirichlet-Multinomial distribution.
    
    :param data: 2D numpy array where each row is a sample and each column is a species
    :return: tuple (alpha_values, alpha_0)
    """
    mean_proportions = np.mean(data, axis=0)
    variance_proportions = np.var(data, axis=0)
    
    # Avoid division by zero
    mask = variance_proportions != 0
    alpha = np.zeros_like(mean_proportions)
    alpha[mask] = mean_proportions[mask] * (((mean_proportions[mask] * (1 - mean_proportions[mask])) / variance_proportions[mask]) - 1)
    
    # Handle cases where variance is zero
    alpha[~mask] = 1e6  # Set a large value for species with zero variance
    
    alpha_0 = np.sum(alpha)
    
    return alpha, alpha_0