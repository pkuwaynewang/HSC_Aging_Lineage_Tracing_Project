import sys

def merge_files(file1, file2, output_file):
    # Open both files
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        # Read all lines from both files
        file1_lines = f1.readlines()
        file2_lines = f2.readlines()

        # Prepare to store the results
        results = []

        # Iterate over the lines from both files
        for elem1, line2 in zip(file1_lines, file2_lines):
            # Strip newline characters and split file2 elements
            elements2 = line2.strip().split(',')
            elem1 = elem1.strip()

            # Append paired elements to results
            for element in elements2:
                results.append([element, elem1])

    # Write the results to a new file
    with open(output_file, 'w') as output:
        for pair in results:
            output.write(f"{pair[0]}\t{pair[1]}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python merge.py <file1.txt> <file2.txt>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = "merged_output.txt"
    merge_files(file1, file2, output_file)
    print("Merge completed. Output file:", output_file)

if __name__ == "__main__":
    main()
