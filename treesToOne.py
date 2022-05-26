import os


def trees_to_one(file_path_1, file_path_2, num_trees_per_group):
    full_Tree1 = open("allTreesGroupOne.txt", "w")  # create a file called allTrees1
    full_Tree2 = open("allTreesGroupTwo.txt", "w")  # create a file called allTrees2
    # Don't use any names that contain escape characters \t or numbers like \001

    dir1 = "{file_path}"
    dir1 = dir1.format(file_path=file_path_1)

    dir2 = "{file_path}"
    dir2 = dir2.format(file_path=file_path_2)

    # dir1 = "simphy_test_program\simphy_test_program\practice2\p001"
    # dir2 = "simphy_test_program\simphy_test_program\practice2\p002"

    for i in range(num_trees_per_group):
        files1 = (os.listdir(dir1))[i]
        files2 = (os.listdir(dir2))[i]
        # will create an array of files in the directory use the index to choose how many files are copied

        file_name1 = dir1 + "\{file}"
        file1 = open(file_name1.format(file=files1))
        file_name2 = dir2 + "\{file}"

        file2 = open(file_name2.format(file=files2))

        for line in file1:  # copy over the file to the full tree file1
            full_Tree1.write(line)

        for line in file2:  # copy over the file to the full tree file2
            full_Tree2.write(line)
