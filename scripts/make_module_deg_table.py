import pandas as pd
import numpy as np
import os

genes_colors_path = "/home/alder/Documents/research/Blueberry_Data/"
data_folder = "/home/alder/Documents/research/Blueberry_Data/Diff_Ex/EdgeR_Output/"
hap = "All_Hap"
method = "Bonferroni"
combined_path = data_folder + hap + '/' + method + '/'

data  = pd.read_csv(genes_colors_path + "Genes_and_ModuleColors.tsv", delimiter="\t")

# Create the result dataframe
result = data.groupby("moduleColor", as_index=False).count()
result.rename({"moduleColor" : "Module ID", "Gene_Names" : "Number of Total Genes"},
    axis = 1, inplace = True)

# TODO: result2 and result3 not implemented yet

# I intend for this result to have the actual sum of dif ex.
# So that if the dif ex data had [-1,-1,1] it would return -1 instead of 3
result2 = result.copy()

# This one will have the the sum of the absolute values for result and result2.
# Then one more column that shows the percentage for that module.
result3 = result.copy()

for f in os.listdir(combined_path):
    # Only use tsv files
    if f[-4:] != ".tsv":
        continue
    
    # Reset the data for the new iteration
    genes_colors = data.copy()
    genes_colors2 = data.copy()

    # Read in the comparison data

    df = pd.read_csv(combined_path + f, delimiter = '\t')

    # Replace the negatives since we only care about if it is differentially
    # expressed or not.
    genes_colors["Diff Ex"] = df["Direction_Differentially_Regulated"].replace({-1:1})
    genes_colors2["Diff Ex"] = df["Direction_Differentially_Regulated"]

    result["Diff_Ex_" + f[:-4]] = genes_colors[["moduleColor","Diff Ex"]].groupby("moduleColor",
        as_index =False).sum()["Diff Ex"]
    result2["Diff_Ex_" + f[:-4]] = genes_colors2[["moduleColor","Diff Ex"]].groupby("moduleColor",
        as_index=False).sum()["Diff Ex"].abs()  # NOTE: The absolute value accounts
                                                # for something.

# Result is the sum of absolute values.
# Result2 is the absolute value of sums.




print(result)
result.to_csv("modules_diff_ex_table.tsv", sep="\t")

print(np.sum(np.abs([-1,1,-5,5])))
result["total"] = result.drop(columns = ["Module ID", "Number of Total Genes"]).apply(np.sum,axis = 1)
result2["total"] = result2.drop(columns = ["Module ID", "Number of Total Genes"]).apply(np.sum,axis = 1)
print(result["total"])
print(result2["total"])


result3["Diff_Ex_Absolute"] = result["total"]
result3["Diff_Ex_Cancellation"] = result2["total"]

result3["% Expressed Same Direction"] = 1 - (result["total"] - result2["total"]) / 2 / result["total"]
print(result3)

result3.to_csv("modules_diff_ex_summary.tsv", sep="\t")
result3.dropna(inplace = True)
result3.to_csv("modules_diff_ex_summary_shortened.tsv", sep="\t")
