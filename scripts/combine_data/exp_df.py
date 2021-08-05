import pandas as pd
import os

#This path contains all of the experiments between drapier, liberty, and the
#control and treatment sets.
#MAGIC:
experimental_path = "/home/alder/Documents/research/Blueberry_Data/Diff_Ex/EdgeR_Output/All_Hap/Bonferroni"

#This path leads to the file that tells what module every blueberry gene
#belongs to.
#MAGIC:
blueberry_module_path = "/home/alder/Documents/research/Blueberry_Data/Blueberry_Data/Genes_and_ModuleColors.tsv"

def add_data(table,direction_file,key):
    #Adds three columns to the table for the direction file.
    name = direction_file[0:-14]
    xcol1 = name + " Total Diff Ex"
    xcol2 = name + " Upregulated"
    xcol3 = name + " Downregulated"
    table[xcol1] = 0
    table[xcol2] = 0
    table[xcol3] = 0

    exp_data = pd.read_csv(experimental_path + "/" +  direction_file, sep = "\t",index_col = "Gene_Name")
    exp_data = exp_data[exp_data["Direction_Differentially_Regulated"] != 0]
    for index, row in exp_data.iterrows():
        table[xcol1][key.loc[index]] += 1
        if row[0] == 1:
            table[xcol2][key.loc[index]] += 1
        if row[0] == -1:
            table[xcol3][key.loc[index]] += 1

    return table

def create_table():

    key_df = pd.read_csv(blueberry_module_path, sep = "\t",index_col = "Gene_Names")

    #This will be the basis for the table. From here the code will add columns
    #three at a time. This will happen for each experimental file.
    #NOTE: this is the number of blueberry genes, not arabidopsis genes.
    table = pd.DataFrame(key_df.value_counts())
    table.columns = ["# of genes"]
    i = 0
    for direction_file in os.listdir(experimental_path):
        if direction_file.endswith("Direction.tsv") != True:
            #This prevents it from reading the wrong files.
            continue
        print(i * progress_step , "%")
        i += 1
        table = add_data(table,direction_file,key_df)

    return table,key_df

progress_step = 100.0/28 #Unnecesarry and magic, but this shows the progress.
table,key_df = create_table()
#print(key_df.loc["augustus_masked-VaccDscaff1-processed-gene-1.5"]["moduleColor"])
table.to_csv("Diff_Ex_Table.csv",sep = "\t")
print(table)
