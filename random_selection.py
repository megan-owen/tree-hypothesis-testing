import random # Import the random module for random sampling

def select_random_lines(file_path, n): #Define a function to select n random lines from a file
    with open(file_path, 'r') as file: # Open the specified file in read mode
        lines = file.readlines() # Read all lines from the file into a list
    return random.sample(lines, n) # Return a list of n randomly sampled lines from the file

def save_lines_to_file(lines, output_file_path): # Define a function to save lines to a specified file
    with open(output_file_path, 'w') as file: # Open the specified file in write mode
        for line in lines[:-1]:  # Process all lines except the last one
            file.write(line.strip() + '\n') # Write each line without leading/trailing whitespace, followed by a newline
        if lines:  # Check if there are any lines
            file.write(lines[-1].strip())  # Write the last line without an additional newline
