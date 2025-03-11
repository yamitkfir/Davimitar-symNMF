# Python interface of your code
import numpy as np
import sys
import symnmfmodule  # Import our C module

RANDOM_SEED = 1234
ERROR_MSG = "An Error Has Occurred"
SEPERATOR = ','

def main():
    np.random.seed(RANDOM_SEED)
    # Read CMD args
    if len(sys.argv) != 4:
        print(ERROR_MSG)
        sys.exit(1)
    try:
        k = int(sys.argv[1])
    except ValueError:
        print(ERROR_MSG)
        sys.exit(1)
    goal = sys.argv[2]
    file_name = sys.argv[3]
    
    # Check if goal is valid
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print(ERROR_MSG)
        sys.exit(1)
    
    # Read data points from file
    try:
        data_points = np.loadtxt(file_name, delimiter = SEPERATOR)
        # np.loadtxt() is a NumPy function that
        # reads data from a text file and creates a NumPy array.
    except:
        print(ERROR_MSG)
        sys.exit(1)
    
    # Check if k is valid
    if k >= len(data_points) or k <= 0:
        print(ERROR_MSG)
        sys.exit(1)
    
    # Call fitting function according to goal
    try:
        if goal == "sym":
            result = symnmfmodule.sym(data_points.tolist())
        elif goal == "ddg":
            result = symnmfmodule.ddg(data_points.tolist())
        elif goal == "norm":
            result = symnmfmodule.norm(data_points.tolist())
        elif goal == "symnmf":
            # For symnmf, first of all get the normalized similarity matrix W
            W = symnmfmodule.norm(data_points.tolist())
            
            # Initialize H as described in 1.4.1
            n = len(data_points)
            m = np.mean(W)
            H_init = np.random.uniform(0, 2 * np.sqrt(m/k), size=(len(X), k))
            
            # Call symnmf function with initial H and W
            result = symnmfmodule.symnmf(W, H_init.tolist())
    except:
        print(ERROR_MSG)
        sys.exit(1)
    
    # Pring the resulting matrix
    for row in result:
        print(SEPERATOR.join(["%.4f" % val for val in row]))

if __name__ == "__main__":
    main()