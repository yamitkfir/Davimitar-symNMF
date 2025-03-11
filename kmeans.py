DEFAULT_ITER = 300
MAX_ITER = 1000
EPSILON = 0.0001
ERROR_MSG = "An Error Has Occurred"
ITER_MSG = "Invalid maximum iteration!"
CLUSTER_MSG = "Invalid number of clusters!"

class Node:
    """Node for linked list implementation"""
    def __init__(self, vector=None):
        self.vector = vector
        self.next = None


class LinkedList:
    """Linked list data structure for points and centroids"""
    def __init__(self):
        self.head = None
        self.size = 0
    
    def append(self, vector):
        new_node = Node(vector)
        if not self.head:
            self.head = new_node
        else:
            current = self.head
            while current.next:
                current = current.next
            current.next = new_node
        self.size += 1
    
    def get_at_index(self, index):
        if index >= self.size:
            return None
        current = self.head
        for _ in range(index):
            current = current.next
        return current.vector


def read_file(file_name):
    """Read points from file into linked list"""
    original_points = LinkedList()
    try:
        with open(file_name, 'r') as f:
            for line in f:
                if line.strip():  # Skip empty lines
                    vector = [float(x) for x in line.strip().split(',')]
                    original_points.append(vector)
        return original_points
    except:
        print(ERROR_MSG)
        sys.exit(1)


def dist_mju(vector1, vector2):
    """Calculate Euclidean distance between two vectors"""
    if len(vector1) != len(vector2):
        return float('inf')
    return sum((a - b) ** 2 for a, b in zip(vector1, vector2)) ** 0.5


def assign_cluster(original_points, current_centroids, K):
    """Assign points to nearest centroid"""
    clusters = [LinkedList() for _ in range(K)]  # initialize K empty LL to be clusters
    current = original_points.head
    
    while current:
        min_dist = float('inf')
        closest_cluster = 0
        
        # Find nearest centroid
        centroid_node = current_centroids.head
        for i in range(K):
            dist = dist_mju(current.vector, centroid_node.vector)
            if dist < min_dist:
                min_dist = dist
                closest_cluster = i
            centroid_node = centroid_node.next
        
        # Assign to nearest cluster
        clusters[closest_cluster].append(current.vector)
        current = current.next
    
    return clusters


def update_cluster(clusters, dimension, K):
    """Calculate new centroids based on clusters"""
    next_centroids = LinkedList()
    
    for cluster in clusters:
        if cluster.size == 0:
            continue
            
        # Calculate mean for each dimension
        new_centroid = [0.0] * dimension
        count = 0
        current = cluster.head
        
        while current:
            for i in range(dimension):
                new_centroid[i] += current.vector[i]
            count += 1
            current = current.next
        
        # Compute average
        if count > 0:
            new_centroid = [x/count for x in new_centroid]
            next_centroids.append(new_centroid)
    
    return next_centroids


def kmeans_general(K, max_iter, input_filename, epsilon = EPSILON):
    # Read data
    original_points = read_file(input_filename)
    if original_points.size == 0:
        print(ERROR_MSG)
        sys.exit(1)
    
    # Validate K
    if K <= 1 or K >= original_points.size:
        print(CLUSTER_MSG)
        sys.exit(1)
    
    # Initialize centroids with first K points
    current_centroids = LinkedList()
    for i in range(K):
        vector = original_points.get_at_index(i)
        if vector:
            current_centroids.append(vector.copy())
    
    # Main iteration loop
    iteration = 0
    while iteration < max_iter:
        # Assign points to clusters
        clusters = assign_cluster(original_points, current_centroids, K)
        
        # Update centroids
        dimension = len(original_points.head.vector)
        next_centroids = update_cluster(clusters, dimension, K)
        
        # Check convergence
        max_delta = 0.0
        curr_centroid = current_centroids.head
        next_centroid = next_centroids.head
        
        while curr_centroid and next_centroid:
            delta = dist_mju(curr_centroid.vector, next_centroid.vector)
            max_delta = max(max_delta, delta)
            curr_centroid = curr_centroid.next
            next_centroid = next_centroid.next
        
        if max_delta < epsilon:
            break
        
        # Update centroids for next iteration
        current_centroids = next_centroids
        iteration += 1
    
    return current_centroids


def print_results(centroids):
    """Print final centroids in desired format"""
    current = centroids.head
    while current:
        print(','.join(f'{x:.4f}' for x in current.vector))
        current = current.next


def main():
    # input must be 3 or 4 arguments
    if len(sys.argv) not in [3, 4]:
        print(ERROR_MSG)
        sys.exit(1)
        
    try:
        K = int(sys.argv[1])
    except ValueError:
        print(CLUSTER_MSG)
        sys.exit(1)
    try:
        # "if iter is not provided, use default value"
        iter = int(sys.argv[2]) if len(sys.argv) == 4 else DEFAULT_ITER
    except ValueError:
        print(ITER_MSG)
        sys.exit(1)
        
    try:
        filename = sys.argv[-1]
        
        # Validate iter
        if iter <= 1 or iter >= MAX_ITER:
            print(ITER_MSG)
            sys.exit(1)
            
        # Run algorithm
        final_centroids = kmeans_general(K, iter, filename)
        
        # Print results
        print_results(final_centroids)
        
    except ValueError:
        print(ERROR_MSG)
        sys.exit(1)


if __name__ == "__main__":
    import sys
    main()