import random


#This function selects n lines randomly of a text file.
#Input = path to the file from which you desire to extract n random lines, and number of lines.
#Output = List of random lines

def select_random_lines(file_path, n):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return random.sample(lines, n)

#This function writes each element of a list in a new line of a text file.
#Input = list, and path where it is going to save the output text file.
#Output = no output.
   
def save_lines_to_file(lines, output_file_path):
    with open(output_file_path, 'w') as file:
        for line in lines[:-1]:  # Process all lines except the last one
            file.write(line.strip() + '\n')
        if lines:  # Check if there are any lines
            file.write(lines[-1].strip())  # Write the last line without an additional newline
