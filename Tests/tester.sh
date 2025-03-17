#!/bin/bash

# Configuration constants
SRC_DIR="/home/developer/sp/Davimitar-symNMF" # Directory containing source code files
# INPUT_DIR="${SRC_DIR}/tests/HW1_tests" # Directory containing input files
INPUT_DIR="${SRC_DIR}/tests/Claude" # Directory containing input files

# Input files to test
INPUT_FILES=("input_1.txt" "input_2.txt")

# Cluster values to test for symnmf
K_VALUES=(2 3)

OUTPUT_DIR="${INPUT_DIR}/expected_output" # Directory for expected output files
ACTUAL_OUTPUT_DIR="${INPUT_DIR}/test_output" # Directory for actual test outputsGREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
RESET='\033[0m'

# Function to run a test
run_test() {
    local command=$1
    local expected_file=$2
    local test_name=$(basename "$expected_file" .out)
    local output_file="${ACTUAL_OUTPUT_DIR}/${test_name}.txt"
    
    echo -e "\nRunning: ${command}"
    
    # Create test output directories if they don't exist
    mkdir -p "${ACTUAL_OUTPUT_DIR}"
    mkdir -p "$(dirname "${expected_file}")"
    
    # Run the command and capture stdout to file while displaying it
    eval $command | tee "${output_file}" 2> "${output_file}.err"
    
    # Check if there was an error
    if [ -s "${output_file}.err" ]; then
        echo -e "${RED}Command produced errors:${RESET}"
        cat "${output_file}.err"
    fi
    
    # Verify that the output file was created and has content
    if [ ! -s "${output_file}" ]; then
        echo -e "${RED}Error: Output file was not created or is empty: ${output_file}${RESET}"
        return 1
    else
        echo -e "Output saved to: ${output_file}"
    fi
    
    # Check if expected output file exists     
    if [ -f "$expected_file" ]; then         
        # Compare the output with the expected output         
        if diff -q "${output_file}" "$expected_file" &>/dev/null; then             
            echo -e "${GREEN}Test passed!${RESET}"         
        else             
            echo -e "${RED}Test failed!${RESET}"             
            echo "Output difference:"             
            echo -e "${YELLOW}Expected output:${RESET}"             
            cat "$expected_file"             
            echo -e "${YELLOW}Actual output:${RESET}"             
            cat "${output_file}"         
        fi
    else
        # If expected file doesn't exist, create it for future reference
        echo -e "${YELLOW}First run - creating expected output file: $expected_file${RESET}"
        cp "${output_file}" "$expected_file"
        echo -e "${GREEN}Created reference file for future tests.${RESET}"
    fi
    
    # Clean up error file if it's empty
    if [ ! -s "${output_file}.err" ]; then
        rm "${output_file}.err"
    fi
}

# Check that required executables exist
if ! [ -f "${SRC_DIR}/symnmf.py" ] || ! [ -f "${SRC_DIR}/symnmf.c" ]; then
    echo -e "${RED}Error: Required source files (symnmf.py or symnmf.c) not found in ${SRC_DIR}${RESET}"
    exit 1
fi

# Create output directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${ACTUAL_OUTPUT_DIR}

# Clear existing test output directory to avoid confusion with old results
echo "Clearing previous test output files..."
rm -f ${ACTUAL_OUTPUT_DIR}/*.txt ${ACTUAL_OUTPUT_DIR}/*.err

# Compile C module first
echo "Compiling C module..."
cd ${SRC_DIR} && python3 setup.py build_ext --inplace

# Make the symnmf C executable
echo "Compiling symnmf executable..."
cd ${SRC_DIR} && make

# Debug - check if the input files exist
echo -e "\n----- Checking for input files -----"
for input_file in "${INPUT_FILES[@]}"; do
    if [ -f "${INPUT_DIR}/${input_file}" ]; then
        echo -e "${GREEN}Found input file: ${INPUT_DIR}/${input_file}${RESET}"
        # Show first few lines
        head -n 3 "${INPUT_DIR}/${input_file}"
    else
        echo -e "${RED}Input file not found: ${INPUT_DIR}/${input_file}${RESET}"
    fi
done

# Run tests for the C executable
echo -e "\n----- Testing C executable -----"
for input_file in "${INPUT_FILES[@]}"; do
    if [ ! -f "${INPUT_DIR}/${input_file}" ]; then
        echo -e "${RED}Skipping tests for missing input file: ${input_file}${RESET}"
        continue
    fi
    
    # Test sym goal
    # run_test "${SRC_DIR}/symnmf sym ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/c_sym_${input_file%.txt}.out"
    
    # Test ddg goal
    # run_test "${SRC_DIR}/symnmf ddg ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/c_ddg_${input_file%.txt}.out"
    
    # Test norm goal
    # run_test "${SRC_DIR}/symnmf norm ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/c_norm_${input_file%.txt}.out"
done

# Run tests for the Python module
echo -e "\n----- Testing Python module -----"
for input_file in "${INPUT_FILES[@]}"; do
    if [ ! -f "${INPUT_DIR}/${input_file}" ]; then
        echo -e "${RED}Skipping tests for missing input file: ${input_file}${RESET}"
        continue
    fi
    
    # Test sym goal
    # run_test "python3 ${SRC_DIR}/symnmf.py 1 sym ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/py_sym_${input_file%.txt}.out"
    
    # Test ddg goal
    # run_test "python3 ${SRC_DIR}/symnmf.py 1 ddg ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/py_ddg_${input_file%.txt}.out"
    
    # Test norm goal
    # run_test "python3 ${SRC_DIR}/symnmf.py 1 norm ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/py_norm_${input_file%.txt}.out"
    
    # For symnmf, we need to specify k (number of clusters)
    for k in "${K_VALUES[@]}"; do
        run_test "python3 ${SRC_DIR}/symnmf.py $k symnmf ${INPUT_DIR}/${input_file}" "${OUTPUT_DIR}/py_symnmf_k${k}_${input_file%.txt}.out"
    done
done

# Count the number of files created
num_files=$(find ${ACTUAL_OUTPUT_DIR} -name "*.txt" | wc -l)
echo -e "\nAll tests completed! Created ${num_files} output files in ${ACTUAL_OUTPUT_DIR}"

# # List the created files
# echo -e "\nCreated test output files:"
# find ${ACTUAL_OUTPUT_DIR} -name "*.txt" | sort | while read file; do
#     file_size=$(du -h "$file" | cut -f1)
#     echo "  - $(basename "$file") ($file_size)"
done