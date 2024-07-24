import random

def select_random_lines(file_path, n):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return random.sample(lines, n)

def save_lines_to_file(lines, output_file_path):
    with open(output_file_path, 'w') as file:
        for line in lines[:-1]:  # Process all lines except the last one
            file.write(line.strip() + '\n')
        if lines:  # Check if there are any lines
            file.write(lines[-1].strip())  # Write the last line without an additional newline
