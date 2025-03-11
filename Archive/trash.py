# Analyze the algorithm.

import sys
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
import symnmf

SEED = 1234
ERROR_MSG = "An Error Has Occurred"
SEPERATOR = ','

def load_data(filename):
    """
    Load data points from the given file.
    Args: filename: Path to the input file
    Returns: numpy.ndarray: Data points
    """
    with open(filename, 'r') as f:
        return np.array([list(map(float, line.strip().split(SEPERATOR))) for line in f])

def symnmf_clustering(X, k):
    """
    Perform SymNMF clustering.
    
    Args:
        X: Data points
        k: Number of clusters
    
    Returns:
        numpy.ndarray: Cluster assignments for each data point
    """
    # Get the normalized similarity matrix
    W = symnmf.norm(X)
    
    # Initialize H as specified in section 1.4.1
    np.random.seed(SEED)
    m = np.mean(W)
    H = np.random.uniform(0, 2 * np.sqrt(m/k), size = (len(X), k))
    
    # Apply SymNMF to get the association matrix H
    H = symnmf.symnmf(H, W)
    
    # Get cluster assignments based on maximum association score (see 1.5)
    labels = np.argmax(H, axis=1)
    return labels

def main():
    """    
    Expects two command line arguments:
    1. k: Number of clusters
    2. filename: Path to the input file
    """
    # Parse CMD args
    if len(sys.argv) != 3:
        print(ERROR_MSG)
        return
    
    try:
        k = int(sys.argv[1])
    except ValueError:
        print(ERROR_MSG)
        return
    filename = sys.argv[2]
    
    try:
        X = load_data(filename)
    except:
        print(ERROR_MSG)
        return
    
    # Check if k is valid
    if k < 1 or k >= len(X):
        print(ERROR_MSG)
        return
    
    # Perform K-means clustering and calculate silhouette score
    try:
        kmeans_labels = kmeans_general(X, k) #FIX yamit
        kmeans_score = silhouette_score(X, kmeans_labels)
    except:
        print(ERROR_MSG)
        return
    
    # Perform SymNMF clustering and calculate silhouette score
    try:
        symnmf_labels = symnmf_clustering(X, k)
        symnmf_score = silhouette_score(X, symnmf_labels)
    except:
        print(ERROR_MSG)
        return
    
    # Print the results
    print(f"nmf: {symnmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")

if __name__ == "__main__":
    main()