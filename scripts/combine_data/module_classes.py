import pandas as pd
import os


def fix_ID(str1):
    return str1.split(".")[1]


def function_to_generate_all_modules(PATH):
    """
    TODO get input args
    """
    # PATH = "../../../Blueberry_Data/Enrichment_Data/Gene_Analysis/"
    modules = []
    for directory in os.listdir(PATH):
        EnrichmentFactory.PATH = PATH + directory + "/"
        for file_ in os.listdir(PATH + directory):
            modules.append(EnrichmentFactory.build_data(file_))
    return modules


class EnrichmentFactory:
    PATH = ""

    @staticmethod
    def build_data(file_string):
        split_string = file_string.split(".")
        module_type = split_string[0]
        type_name = split_string[1]
        file_type = split_string[2]
        if file_type != "tsv":
            raise NotImplementedError("File type is %s instead of tsv" % type_name)
        return EnrichmentData(module_type, file_string, type_name)


class EnrichmentData:
    def __init__(self, module, file_string, type_name):
        self.module = module
        self.file_string = file_string
        self.csv_in()
        self.file_type = type_name

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentFactory.PATH + self.file_string, delimiter="\t")
        self.df = self.df.loc[
            self.df["false discovery rate"] < 0.05
        ]  # MAGIC: Column name and fdr cutoff.
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.
        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")
        self.df[column1] = self.df[column1].apply(lambda x: [fix_ID(i) for i in x])

    def csv_out(self):
        # not sure about this one
        return 0
