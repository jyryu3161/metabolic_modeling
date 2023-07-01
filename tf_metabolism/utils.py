import tempfile
import argparse
import copy
import os
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model

def argument_parser(version=None):
    parser = argparse.ArgumentParser()    
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-c1', '--input_omics_file1', required=True, help="Input omics file for condition 1")
    parser.add_argument('-c2', '--input_omics_file2', required=True, help="Input omics file for condition 2")
    parser.add_argument('-r', '--related_sample', default=False, action='store_true', help="Flag for related sample")

    return parser

def update_cobra_model(cobra_model):
    temp = tempfile.NamedTemporaryFile(prefix='temp_cobra_', suffix='.xml', dir='./', delete=True)
    temp_outfile = temp.name
    temp.close()

    write_sbml_model(cobra_model, temp_outfile, use_fbc_package=False)
    cobra_model = read_sbml_model(temp_outfile)
    os.remove(temp_outfile)
    return cobra_model

def check_input_file_format(omics_file1, omics_file2):
    
    
    return