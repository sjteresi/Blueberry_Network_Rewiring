# This file is to analyze the process file. Both as a whole, and within
# individual modules.
import os
from combine_data import modules

# from module_classes import *
from module_classes import function_to_generate_all_modules

"""
Each data has:
    self.module - the module color it belongs to
    self.df - the dataframe it has
    self.file - the file name. Not that important.
    self.file_type   = This is the file type. This
        seems like it is a bit of a hack though.
"""
# Takes in a dictionary, sorts it by the value from least to greatest, then
# returns an array of lists.
def sortByValue(my_dict):
    """
    This allows me to sort a dictionary. It doesn't *really* sort a
    dictionary, but it basically does.
    """
    sortedByValue = sorted(my_dict.items(), key=lambda t: t[1])
    return sortedByValue


# Makes a file listing the different proteins found in a module.
def module_data(data, folder="Module_Data/"):
    if os.path.exists(folder) == False:  # We gotta create the folder if it
        os.mkdir(folder)  # doesn't exist.
    # Creates a file in the folder we want
    Unique_Terms = open(
        folder + "genes_" + data.module + "_" + data.file_type + ".txt", "w"
    )
    Gene_Counts = {}
    for genes in data.df["matching proteins in your network (IDs)"]:
        # Code below counts up all the genes in each row of the module.
        # Not the most efficient, but higher efficiency is un-needed.
        for gene in genes:
            if gene not in Gene_Counts:
                Gene_Counts[gene] = 1
            else:
                Gene_Counts[gene] += 1
    Sorted_Gene_Counts = sortByValue(Gene_Counts)
    string_of_genes = ""
    for Gene in reversed(Sorted_Gene_Counts):
        string_of_genes += Gene[0] + ": " + str(Gene[1]) + "\n"
    Unique_Terms.write(string_of_genes)
    Unique_Terms.close()


# Makes a file with the term description and number of times shown up.
def term_data(list_of_files, file_type="Process"):
    # First, we make the different containers for our data
    Unique_Terms = open("term_data.txt", "w")
    Term_Counts = {}
    # Second, we fill in all our data
    for data in list_of_files:
        if data.file_type != file_type:
            # This skips the current file if it isn't a process file.
            continue
        for desc in data.df["term description"]:
            if desc not in Term_Counts:
                Term_Counts[desc] = 1
            else:
                Term_Counts[desc] += 1
    # Third, we sort our data
    Sorted_Unique_Terms = sortByValue(Term_Counts)
    # Lastly, we write our data
    string_of_unique_terms = ""
    for desc in reversed(Sorted_Unique_Terms):
        string_of_unique_terms += desc[0] + ":\t" + str(desc[1]) + "\n"
    Unique_Terms.write(string_of_unique_terms)
    Unique_Terms.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="...")
    parser.add_argument("generation_path", type = str, help = """Path to the
                        folder used for generatng the classes.""")
    parser.add_argument("output_folder",type = str, help = """Path to the
                        folder that the data will be outputted to."""   )
    args = parser.parse_args()
    generation_path = args.generation_path
    output_folder = args.output_folder

    modules = function_to_generate_all_modules(generation_path):

    #term_data(modules)  # TODO edit?
    for data in modules:  # TODO edit?
        # We only want to work with the Process files.
        if data.file_type != "Process":
            continue
        # Loops through every file. One file at a time.
        # The code below is supposed to count how many times any given gene shows up.
        module_data(data, output_folder)
