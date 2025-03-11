# UNDONE

import sys
import numpy as np
from sklearn.metrics import silhouette_score
import symnmf
import kmeans

# Constants
ERROR_MSG = "An Error Has Occurred"
SEPARATOR = ","
SEED = 1234

def load_data(filename):
    """    
    Returns: numpy.ndarray: Data points
    """
    try:
        return np.loadtxt(filename, delimiter=SEPARATOR)
    except:
        print(ERROR_MSG)
        sys.exit(1)

def kmeans_clustering(X, k, filename):
    """
    Returns: numpy.ndarray: Cluster assignments for each data point
    """
    # Convert numpy array to format expected by kmeans module
    centroids = kmeans.kmeans_general(k, kmeans.DEFAULT_ITER, filename)
    
    # Get cluster assignments by finding the nearest centroid for each point
    labels = []
    for point in X:
        min_dist = float('inf')
        nearest_cluster = 0
        
        # Convert centroids linked list to a list for easier processing
        centroid_vectors = []
        current = centroids.head
        while current:
            centroid_vectors.append(current.vector)
            current = current.next
        
        # Find the nearest centroid
        for i, centroid in enumerate(centroid_vectors):
            # Use the distance function from kmeans module
            dist = kmeans.dist_mju(list(point), centroid)
            if dist < min_dist:
                min_dist = dist
                nearest_cluster = i
        labels.append(nearest_cluster)
    
    return np.array(labels)

def symnmf_clustering(X, k):
    """
    Returns: numpy.ndarray: Cluster assignments for each data point
    """
    # Get the normalized similarity matrix
    W = symnmfmodule.norm(X.tolist())
    
    # Initialize H as specified in section 1.4.1
    np.random.seed(SEED)
    m = np.mean(W)
    H_init = np.random.uniform(0, 2 * np.sqrt(m/k), size=(len(X), k))
    
    # Apply SymNMF to get the association matrix H
    H = symnmfmodule.symnmf(W, H_init.tolist())
    
    # Get cluster assignments based on maximum association score (section 1.5)
    labels = np.argmax(H, axis = 1)
    return labels

def main():
    """    
    Expects 2 CMD args:
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
        kmeans_labels = kmeans_clustering(X, k, filename)
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
        sys.exit(1)
    
    # Print the results
    print(f"nmf: {symnmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")

if __name__ == "__main__":
    main()