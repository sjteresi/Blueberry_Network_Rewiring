from combine_data import modules
"""
Each data has:
    self.module - the module color it belongs to
    self.df - the dataframe it has
    self.file - the file name. Not that important.
    type(self).__class__   = This is the file type. This
        seems like it is a bit of a hack though.
"""

Gene_Counts = {}
total = 0.0
n = 0.0
too_high = 0
for data in modules:
    #Loops through every file. One file at a time.
    #print(data.module)

    #The code below is supposed to count how many times any given gene shows up.
    if data.
    for genes in data.df["matching proteins in your network (labels)"]:
        print(genes)
        for gene in genes:
            print(gene)
            if gene not in Gene_Counts:
                Gene_Counts[gene] = 1
            else:
                Gene_Counts[gene] += 1
"""
    #I couldn't think of anything else that we would be looking for.
    for fdr in data.df["false discovery rate"]:
        n += 1
        total += fdr
        if fdr >= .05:
            too_high += 1

#print(Gene_Counts)
print(total)
print(n)
print(total/n)
print(too_high)
"""

print(Gene_Counts)
