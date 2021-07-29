# This file is to analyze the process file. Both as a whole, and within
# individual modules.
from combine_data import modules
"""
Each data has:
    self.module - the module color it belongs to
    self.df - the dataframe it has
    self.file - the file name. Not that important.
    self.file_type   = This is the file type. This
        seems like it is a bit of a hack though.
"""
#Takes in a dictionary, sorts it by the value from least to greatest, then
#returns an array of lists.
def sortByValue(my_dict):
    sortedByValue = sorted(my_dict.items(), key = lambda t:t[1])
    return sortedByValue

# Makes a file listing the different proteins found in a module.
def module_data(data):
    Unique_Terms = open("genes_" + data.module + "_" + data.file_type + ".txt",'w')
    Gene_Counts = {}
    for genes in data.df["matching proteins in your network (labels)"]:
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

#Unique_Terms = open("Unique_Terms.txt",'w')
#Gene_Counts = {}
for data in modules:
    #We only want to work with the Process files.
    if data.file_type != "Process":
        continue
    #Loops through every file. One file at a time.
    #The code below is supposed to count how many times any given gene shows up.
    module_data(data)


    """
    for genes in data.df["matching proteins in your network (labels)"]:
        for gene in genes:
            if gene not in Gene_Counts:
                Gene_Counts[gene] = 1
            else:
                Gene_Counts[gene] += 1
    """

"""
Sorted_Gene_Counts = sortByValue(Gene_Counts)
string_of_genes = ""
for Gene in reversed(Sorted_Gene_Counts):
    string_of_genes += Gene[0] + ": " + str(Gene[1]) + "\n"

Unique_Terms.write(string_of_genes)
Unique_Terms.close()
"""
