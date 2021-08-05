import os
import pandas as pd

# Takes in a folder where the files are located, and ouputs a dataframe.
def create_df(folder):
    genes = [] #Unique list of genes
    df_list = [] #List of dataframes (one for each file)

    for file_ in os.listdir(folder):
        if file_.endswith(".tsv"): #This is the file type we want.
            df_list.append(
                pd.value_counts(
                    pd.read_csv(
                        os.path.join(folder,file_),
                        names = [file_[0:-4]]
                    ).iloc[:,0]
                )
            )
            #That's a mouthful! First it checks if it is a tsv. Then it reads it as
            #a csv, joining the path to the folder and the file name. At this
            #point, it is basically a list of genes. From here, we run value counts
            #on it to determine the number of times each gene is found in the
            #module. Lastly, we append it to the list of dataframes.
            genes = set(list(genes) + list(df_list[-1].index)) #Removes Repeats

    #Creates a dataframe with all the possible genes.
    total_df = pd.DataFrame(index = genes)
    for df in df_list:
        total_df = total_df.join(df)
    total_df = total_df.fillna(0) #Replaces the genes that didn't exist in one df with 0s.
    total_df.columns = [x.split("_")[0] for x in total_df.columns]
    return total_df

def save_df(df):
    df.to_csv("module_genes",sep = "\t")

def open_df():
    df = pd.read_csv("module_genes", sep = "\t")
    return df


folder = "/home/alder/Documents/research/Blueberry_Data/WGCNA_Data/modulecolors_AT_w_duplicates"
dataframe = create_df(folder)
#save_df(dataframe)
#df = open_df()
#print(df)
#print(dataframe["bisque4"].sort_value())
print(dataframe.sum())
