import logging
import os
import glob
import time
import pkg_resources
import shutil
import copy
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
import pandas as pd

from tf_metabolism.metabolic_simulation import flux_prediction
from tf_metabolism.metabolic_simulation import flux_sum
from tf_metabolism.metabolic_simulation import Simulator
from tf_metabolism.omics_integration import tINIT
from tf_metabolism import utils

from tf_metabolism.statistical_analysis import statistical_comparison
from tf_metabolism.statistical_analysis import enrichment

from tf_metabolism.metabolic_model import model_editing
from tf_metabolism import __version__

def reconstruct_GEM(biomass_reaction, generic_model_file, universal_model_file, medium_file, output_dir, omics_file, present_metabolite_file, essential_reaction_file, metabolic_task_file):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    use_normalization = False
    integration_method = 'INIT'
    scoring_method = 'rank'
    
    generic_cobra_model = read_sbml_model(generic_model_file)
    universal_model = read_sbml_model(universal_model_file)
    
    generic_cobra_model = model_editing.make_medium_model(generic_cobra_model, medium_file)
    universal_model = model_editing.make_medium_model(universal_model, medium_file)
    
    tINIT.omics_preprocessing(generic_cobra_model, output_dir, omics_file, False)
    tINIT.reconstruct_GEMs(generic_cobra_model, universal_model, biomass_reaction, integration_method, scoring_method, present_metabolite_file, essential_reaction_file, metabolic_task_file, medium_file, output_dir)

def predict_metabolic_fluxes(output_dir, flux_output_dir, model_dir, generic_model_file):
    if not os.path.isdir(flux_output_dir):
        os.mkdir(flux_output_dir)
    
    output_folder_basename = os.path.basename(flux_output_dir)
    model_file_dir = glob.glob('%s/reconstructed_models/*'%(model_dir))
    flux_profiles = {}
    for each_dir in model_file_dir:
        each_model_file = glob.glob(each_dir+'/functional_*.xml')[0]
        basename = os.path.basename(each_model_file).split('.xml')[0].strip()
        sample_name = basename.split('functional_')[1].strip()
        
        omics_file = '%s/splited_omics_data/%s.csv'%(model_dir, sample_name)
        predicted_flux = flux_prediction.calculate_flux(flux_output_dir, omics_file, generic_model_file, each_model_file)   
        flux_profiles[sample_name] = predicted_flux
    
    df = pd.DataFrame.from_dict(flux_profiles)
    df.to_csv(output_dir+'/%s.csv'%(output_folder_basename))
    return

def predict_enriched_metabolic_pathways(output_dir, cobra_model, comparison_result_df):
    up_df = comparison_result_df[comparison_result_df['Status']=='UP']
    up_reactions = list(up_df.index)
    up_pathway_df = enrichment.pathway_enrichment_analysis(cobra_model, up_reactions)
    up_pathway_df = up_pathway_df[up_pathway_df['FDR corrected P']<0.05]
    up_pathway_df.to_csv(output_dir+'/Up_regulated_pathways.csv')
    
    down_df = comparison_result_df[comparison_result_df['Status']=='DOWN']
    down_reactions = list(down_df.index)
    down_pathway_df = enrichment.pathway_enrichment_analysis(cobra_model, down_reactions)
    down_pathway_df = down_pathway_df[down_pathway_df['FDR corrected P']<0.05]
    down_pathway_df.to_csv(output_dir+'/Down_regulated_pathways.csv')
    return
   
def predict_enriched_transcription_factors(transcript_id_info, trrust, cobra_model, output_dir):
    gene_transcript_info = {}
    with open(transcript_id_info, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            gene_id = sptlist[0].strip()
            transcript_id = sptlist[1].strip()
            gene_transcript_info[transcript_id] = gene_id
    
    up_tf_gene_info = {}
    with open(trrust, 'r') as fp:
        fp.readline()
        for line in fp:
            if 'Activation' in line:
                sptlist = line.strip().split('\t')
                tf = sptlist[0].strip()
                gene = sptlist[2].strip()
                if tf not in up_tf_gene_info:
                    up_tf_gene_info[tf] = [gene]
                else:
                    up_tf_gene_info[tf].append(gene)
                    up_tf_gene_info[tf] = list(set(up_tf_gene_info[tf]))
    
    down_tf_gene_info = {}
    with open(trrust, 'r') as fp:
        fp.readline()
        for line in fp:            
            if 'Repression' in line:
                sptlist = line.strip().split('\t')
                tf = sptlist[0].strip()
                gene = sptlist[2].strip()
                if tf not in down_tf_gene_info:
                    down_tf_gene_info[tf] = [gene]
                else:
                    down_tf_gene_info[tf].append(gene)
                    down_tf_gene_info[tf] = list(set(down_tf_gene_info[tf]))
                    
    changed_flux_info_file = output_dir+'/Differentially_changed_fluxes.csv'
    changed_flux_info_df = pd.read_csv(changed_flux_info_file, index_col=0)
    
    ################
    up_reactions = changed_flux_info_df[changed_flux_info_df['Status']=='UP'].index
    up_transcripts = []
    for each_reaction in cobra_model.reactions:
        if each_reaction.id in up_reactions:
            transcripts = [gene.id for gene in each_reaction.genes]
            up_transcripts+=transcripts
    
    up_transcripts = list(set(up_transcripts))
    up_genes = []
    
    for transcript in up_transcripts:
        if transcript in gene_transcript_info:
            up_genes.append(gene_transcript_info[transcript]) 
    up_genes = list(set(up_genes))
    
    tf_enrichment_df = enrichment.tf_enrichment_analysis(up_tf_gene_info, up_genes)
    tf_enrichment_df.to_csv(output_dir+'/TF_up_enrichment.csv')
    
    ################
    down_reactions = changed_flux_info_df[changed_flux_info_df['Status']=='DOWN'].index
    down_transcripts = []
    for each_reaction in cobra_model.reactions:
        if each_reaction.id in down_reactions:
            transcripts = [gene.id for gene in each_reaction.genes]
            down_transcripts+=transcripts
    
    down_transcripts = list(set(down_transcripts))
    down_genes = []
    
    for transcript in down_transcripts:
        if transcript in gene_transcript_info:
            down_genes.append(gene_transcript_info[transcript])
        
    down_genes = list(set(down_genes))
    
    tf_enrichment_df = enrichment.tf_enrichment_analysis(down_tf_gene_info, down_genes)
    tf_enrichment_df.to_csv(output_dir+'/TF_down_enrichment.csv')
    return
    
def main():    
    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    
    parser = utils.argument_parser(version=__version__)
    
    options = parser.parse_args()
    output_dir = options.output_dir
    
    omics_file1 = options.input_omics_file1
    omics_file2 = options.input_omics_file2
    related_sample_flag = False
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    output_dir_c1 = '%s/condition1'%(output_dir)
    output_dir_c2 = '%s/condition2'%(output_dir)
    
    ## Load data    
    metabolic_task_file = pkg_resources.resource_filename('tf_metabolism', 'data/MetabolicTasks.csv')
    medium_file = pkg_resources.resource_filename('tf_metabolism', 'data/RPMI1640_medium.txt')
    present_metabolite_file = pkg_resources.resource_filename('tf_metabolism', 'data/essential_metabolites.txt')
    essential_reaction_file = pkg_resources.resource_filename('tf_metabolism', 'data/essential_reactions.txt')
    
    # Check omics data format and choose proper metabolic model
    generic_model_file_name = utils.check_input_file_format(omics_file1, omics_file2)
    generic_model_file_name = pkg_resources.resource_filename('tf_metabolism', 'data/Recon2M.2_Entrez_Gene.xml')
    generic_model_file = generic_model_file_name
    universal_model_file = generic_model_file_name
    
    biomass_reaction = 'biomass_reaction'
    cobra_model = read_sbml_model(generic_model_file)
    
    # Find differentially expressed genes or transcripts
    omics1_df = pd.read_csv(omics_file1, index_col=0)
    omics2_df = pd.read_csv(omics_file2, index_col=0)
    
    generic_cobra_model = read_sbml_model(generic_model_file)
    taget_genes = []
    for each_gene in generic_cobra_model.genes:
        taget_genes.append(int(each_gene.id))
    
    taget_genes = list(set(taget_genes)&set(omics1_df.index))
    
    omics1_df = omics1_df.loc[taget_genes]
    omics2_df = omics2_df.loc[taget_genes]

    comparison_result_df = statistical_comparison.two_grouped_data_comparison(omics1_df, omics2_df, related_sample_flag, 0.05)
    comparison_result_df.to_csv(output_dir+'/Differentially_expressed_genes.csv')
    
# Reconstruct GEMs
    reconstruct_GEM(biomass_reaction, generic_model_file, universal_model_file, medium_file, output_dir_c1, omics_file1, present_metabolite_file, essential_reaction_file, metabolic_task_file)
    reconstruct_GEM(biomass_reaction, generic_model_file, universal_model_file, medium_file, output_dir_c2, omics_file2, present_metabolite_file, essential_reaction_file, metabolic_task_file)
    
    ## Predict metabolix fluxes using each condition specific GEM
    output_dir_f1 = '%s/flux1'%(output_dir)
    output_dir_f2 = '%s/flux2'%(output_dir)
    predict_metabolic_fluxes(output_dir, output_dir_f1, output_dir_c1, generic_model_file)
    predict_metabolic_fluxes(output_dir, output_dir_f2, output_dir_c2, generic_model_file)
    
    # ## Statistical analysis of metabolic flux
    flux_file1 = '%s/flux1.csv'%(output_dir)
    flux_file2 = '%s/flux2.csv'%(output_dir)
    flux1_df = pd.read_csv(flux_file1, index_col=0)
    flux2_df = pd.read_csv(flux_file2, index_col=0)
    comparison_result_df = statistical_comparison.two_grouped_data_comparison(flux1_df, flux2_df, related_sample_flag, 0.05)
    comparison_result_df.to_csv(output_dir+'/Differentially_changed_fluxes.csv')
    
    # Perform enrichment analysis for metabolic pathways
    predict_enriched_metabolic_pathways(output_dir, cobra_model, comparison_result_df)
    
    ## Calculate flux_sum
    result = flux_sum.calculate_flux_sum(cobra_model, flux_file1, flux_file2, related=related_sample_flag, p_value_cutoff=0.05)
    result.to_csv('%s/flux_sum_results.csv'%(output_dir))
    
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
    
if __name__ == '__main__':
    main()
