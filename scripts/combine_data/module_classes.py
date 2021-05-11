import pandas as pd
from abc import ABC, abstractmethod


class EnrichmentFactory:
    @staticmethod
    def build_data(file_string):
        split_string = files_string.split(".")
        module_type = split_string[0]
        type_name = split_string[1]
        file_type = split_string[2]
        if file_type != "tsv":
            raise NotImplementedError("File type is %s instead of tsv" % type_name)
        if type_name == "Component":
            return Component(module_type, file_string)
        if type_name == "Smart":
            return Smart(module_type, file_string)
        if type_name == "Keyword":
            return Keyword(module_type, file_string)
        if type_name == "Function":
            return Function(module_type, file_string)
        if type_name == "Pfam":
            return Pfam(module_type, file_string)
        if type_name == "Process":
            return Process(module_type, file_string)
        if type_name == "InterPro":
            return InterPro(module_type, file_string)
        if type_name == "RCTM":
            return RCTM(module_type, file_string)
        if type_name == "PMID":
            return PMID(module_type, file_string)
        if type_name == "NetworkNeighborAL":
            return NetworkNeighborAL(module_type, file_string)
        if type_name == "KEGG":
            return KEGG(module_type, file_string)
        # TODO Add More Types

        raise NotImplementedError("No class found for typename: %s" % type_name)


class EnrichmentData(ABC):
    PATH = ""

    @abstractmethod
    def csv_in():
        # reads this type of csv
        pass

    @abstractmethod
    def csv_out(
        self,
    ):
        # self.df.to_csv(self.file +)
        # pd.to_csv(filename,sep = "\t",header = true, index = false)
        pass


class Component(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.)
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.
        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        # not sure about this one
        return 0


class Smart(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class Keyword(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class Function(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class Pfam(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class Process(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class InterPro(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class RCTM(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class PMID(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class NetworkNeighborAL(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0


class KEGG(EnrichmentData):
    def __init__(self, module, file):
        self.module = module
        self.file = file
        self.csv_in()

    def csv_in(self):
        self.df = pd.read_csv(EnrichmentData.PATH + self.file, delimiter="\t")
        column1 = (
            "matching proteins in your network (IDs)"  # MAGIC: Column name from file.
        )
        column2 = "matching proteins in your network (labels)"  # MAGIC: Column name from file.

        self.df[column1] = self.df[column1].str.split(",")
        self.df[column2] = self.df[column2].str.split(",")

    def csv_out(self):
        return 0
