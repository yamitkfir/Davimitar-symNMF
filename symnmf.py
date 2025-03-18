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
        print(ERROR_MSG + "14")
        print(sys.argv)
        sys.exit(1)
    try:
        k = int(sys.argv[1])
    except ValueError:
        print(ERROR_MSG + "19")
        sys.exit(1)
    goal = sys.argv[2]
    file_name = sys.argv[3]
    
    # Check if goal is valid
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print(ERROR_MSG + "26")
        sys.exit(1)
    
    # Read data points from file
    try:
        data_points = np.loadtxt(file_name, delimiter = SEPERATOR)
        # np.loadtxt() is a NumPy function that
        # reads data from a text file and creates a NumPy array.
    except:
        print(ERROR_MSG + "35")
        sys.exit(1)
    
    # Check if k is valid
    if k >= len(data_points) or k <= 0:
        print(ERROR_MSG + "40")
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
            H_init = np.random.rand(n,k)
            for row in H_init:
                die = (SEPERATOR.join(["%.4f" % val for val in row])) 
                # TEMP TODO for some reason only works with this, as if to take a timeout or smth...
            H_init *= 2 * np.sqrt(m / k)
            
            # TEMP Print initial H
            print("Initial H2 matrix:")
            for row in H_init:
                print(SEPERATOR.join(["%.4f" % val for val in row]))

            result = symnmfmodule.symnmf(W, H_init.tolist())
    except Exception as e:
        print(f"{ERROR_MSG}: {e}")
        sys.exit(1)
    
    # Pring the resulting matrix
    for row in result:
        print(SEPERATOR.join(["%.4f" % val for val in row]))

if __name__ == "__main__":
    main()