# Python interface of your code
import math
import numpy as np
import pandas as pd
import sys
import symnmfmodule  # Import our C module

RANDOM_SEED = 1234
ERROR_MSG = "An Error Has Occurred"
SEPERATOR = ','

def initH(n, k, W):
    m = np.mean(W)
    H_init = np.random.uniform(0, 2 * np.sqrt(m/k), size=(n, k))
    return H_init
    # m = np.average(W)
    # high = 2 * math.sqrt(m/k)
    # basel = np.random.uniform(0, np.nextafter(high, high+1), (n, k))
    # return basel

def check_validity(goal, k, data_points):
    '''
    Checks the validity of the goal and k inputs. If one is invalid, prints error and terminates program.
    '''
    # Check if goal is valid
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print(ERROR_MSG)
        sys.exit(1)
    
    # Check if k is valid, given it's an integer
    if  k >= len(data_points) or k <= 0:
        print(ERROR_MSG)
        sys.exit(1)

def main():
    np.random.seed(RANDOM_SEED)
    if len(sys.argv) != 4: # Read CMD args
        print(ERROR_MSG)
        print(sys.argv)
        sys.exit(1)
    try:
        k = int(sys.argv[1])
    except ValueError:
        print(ERROR_MSG)
        sys.exit(1)
    goal = sys.argv[2]
    file_name = sys.argv[3]
    try: # Read data points from file
        data_points = np.loadtxt(file_name, delimiter = SEPERATOR)
    except:
        print(ERROR_MSG)
        sys.exit(1)
    check_validity(goal, k, data_points)
    try: # Call fitting function according to goal
        if goal == "sym":
            result = symnmfmodule.sym(data_points.tolist())
        elif goal == "ddg":
            result = symnmfmodule.ddg(data_points.tolist())
        elif goal == "norm":
            result = symnmfmodule.norm(data_points.tolist())
        elif goal == "symnmf": # For symnmf, first of all get the normalized similarity matrix W
            W = symnmfmodule.norm(data_points.tolist())
            n = len(data_points)
            H_init = initH(n, k, W).tolist()
            result = symnmfmodule.symnmf(W, H_init)
    except Exception as e:
        print(f"{ERROR_MSG}: {e}")
        sys.exit(1)
    for row in result: # Print the resulting matrix
        print(SEPERATOR.join(["%.4f" % val for val in row]))

if __name__ == "__main__":
    main()