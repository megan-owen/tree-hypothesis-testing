import os


def trees_to_one(file_path_1, file_path_2, num_trees_per_group): # put an r in front of the string to ignore escape char
    full_Tree1 = open("allTreesGroupOne.txt", "w")  # create a file called allTrees1
    full_Tree2 = open("allTreesGroupTwo.txt", "w")  # create a file called allTrees2

    dir1 = file_path_1
    dir2 = file_path_2

    # print(dir1)
    # dir1 = r"simphy_test_program\simphy_test_program\testing2\001"
    # dir2 = r"simphy_test_program\simphy_test_program\testing2\002"

    for i in range(num_trees_per_group):
        files1 = (os.listdir(dir1))[i]
        files2 = (os.listdir(dir2))[i]
        # will create an array of files in the directory use the index to choose how many files are copied

        file_name1 = dir1 + r"\{file}"
        file1 = open(file_name1.format(file=files1))  # placing the actual file (.tree) into the string
        file_name2 = dir2 + r"\{file}"

        file2 = open(file_name2.format(file=files2))

        for line in file1:  # copy over the file to the full tree file1
            full_Tree1.write(line)

        for line in file2:  # copy over the file to the full tree file2
            full_Tree2.write(line)