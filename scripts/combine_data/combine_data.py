#!/usr/bin/env python

"""
Combine module enrichment data and differential expression data
"""

__author__ = "Scott Teresi, Alder Fulton"

# IMPORTS go here
import os
from module_classes import * 
PATH = "/home/alder/Documents/research/Blueberry_Data/Enrichment_Data/Gene_Analysis/"

modules = []
for directory in os.listdir(PATH):
    EnrichmentData.PATH = PATH + directory + "/"
    for file_ in os.listdir(PATH+directory):
        modules.append(EnrichmentFactory.build_data(file_) )









