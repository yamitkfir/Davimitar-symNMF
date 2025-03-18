#!/bin/bash

# Define colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
RESET='\033[0m'

# Function to run a test
run_test() {
    local command=$1
    local expected_file=$2
    
    echo -e "\nRunning: ${command}"
    
    # Create a temporary file for the output
    temp_output=$(mktemp)
    
    # Run the command and save output to the temp file
    eval $command > $temp_output
    
    # Check if expected output file exists
    if [ -f "$expected_file" ]; then
        # Compare the output with the expected output
        if diff -q $temp_output $expected_file &>/dev/null; then
            echo -e "${GREEN}Test passed!${RESET}"
        else
            echo -e "${RED}Test failed!${RESET}"
            echo "Output difference:"
            diff $temp_output $expected_file
        fi
    else
        # If expected file doesn't exist, create it for future reference
        echo -e "${RED}Expected file $expected_file does not exist. Creating it with current output.${RESET}"
        cp $temp_output $expected_file
    fi
    
    # Clean up the temporary file
    rm $temp_output
}

# Check that required executables exist
if ! [ -f "symnmf.py" ] || ! [ -f "symnmf.c" ]; then
    echo -e "${RED}Error: Required source files (symnmf.py or symnmf.c) not found${RESET}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p output

# Compile C module first
echo "Compiling C module..."
python3 setup.py build_ext --inplace

# Make the symnmf C executable
echo "Compiling symnmf executable..."
make

# Run tests for the C executable
echo -e "\n----- Testing C executable -----"
for input_file in input_1.txt input_2.txt; do
    # Test sym goal
    run_test "./symnmf sym $input_file" "output/c_sym_${input_file%.txt}.out"
    
    # Test ddg goal
    run_test "./symnmf ddg $input_file" "output/c_ddg_${input_file%.txt}.out"
    
    # Test norm goal
    run_test "./symnmf norm $input_file" "output/c_norm_${input_file%.txt}.out"
done

# Run tests for the Python module
echo -e "\n----- Testing Python module -----"
for input_file in input_1.txt input_2.txt; do
    # Test sym goal
    run_test "python3 symnmf.py sym $input_file" "output/py_sym_${input_file%.txt}.out"
    
    # Test ddg goal
    run_test "python3 symnmf.py ddg $input_file" "output/py_ddg_${input_file%.txt}.out"
    
    # Test norm goal
    run_test "python3 symnmf.py norm $input_file" "output/py_norm_${input_file%.txt}.out"
    
    # For symnmf, we need to specify k (number of clusters)
    for k in 2 3; do
        run_test "python3 symnmf.py $k symnmf $input_file" "output/py_symnmf_k${k}_${input_file%.txt}.out"
    done
done

echo -e "\nAll tests completed!"
