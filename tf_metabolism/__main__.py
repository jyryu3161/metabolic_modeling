import logging
import os
import glob
import time
import numpy as np
import pkg_resources
import joblib
import shutil
import copy
import umap
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
import pandas as pd
import tqdm

from tf_metabolism.metabolic_simulation import flux_prediction
from tf_metabolism.metabolic_simulation import flux_sum
from tf_metabolism.metabolic_simulation import Simulator
from tf_metabolism.omics_integration import tINIT
from tf_metabolism import utils

from tf_metabolism.statistical_analysis import statistical_comparison
from tf_metabolism.statistical_analysis import enrichment

from tf_metabolism.metabolic_model import model_editing
from tf_metabolism import __version__

import matplotlib.pyplot as plt
import seaborn as sns

from adjustText import adjust_text

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
    
def run_targeting_simulation(output_dir, targeting_result_dir):
    target_folders = glob.glob(output_dir+'/condition2/reconstructed_models/*')
    for each_folder in target_folders:
        basename = os.path.basename(each_folder)
        model_file = glob.glob(each_folder+'/functional_%s.xml'%(basename))[0]
        flux_file = glob.glob(output_dir+'/flux2/Flux_prediction_%s.csv'%(basename))[0]

        flux_dist = {}
        with open(flux_file, 'r') as fp:
            for line in fp:
                sptlist = line.strip().split(',')
                reaction = sptlist[0].strip()
                flux = float(sptlist[1].strip())
                flux_dist[reaction] = flux

        cobra_model = read_sbml_model(model_file)

        obj = Simulator.Simulator()
        obj.load_cobra_model(cobra_model)

        results = {}
        for each_gene in tqdm.tqdm(cobra_model.genes):
            results[each_gene] = {}
            flux_constraints = {}
            target_reactions = []
            for each_reaction in each_gene.reactions:
                target_reactions.append(each_reaction.id)
                flux_constraints[each_reaction.id] = [0.0, 0.0]        

            status, obj_value, perturbed_flux_dist = obj.run_MOMA(wild_flux=flux_dist, flux_constraints=flux_constraints)
            results[each_gene] = perturbed_flux_dist

        output_df = pd.DataFrame.from_dict(results)
        output_df.to_csv(targeting_result_dir+'/MOMA_target_results_%s.csv'%(basename))
    return

# 이상치를 찾기 위한 IQR 계산 함수
def find_outliers_iqr(df):
    outliers_indices = []
    for column in df.columns:
        Q1 = df[column].quantile(0.25)
        Q3 = df[column].quantile(0.75)
        IQR = Q3 - Q1
        outlier_step = 1.5 * IQR
        outliers_list_col = df[(df[column] < Q1 - outlier_step) | (df[column] > Q3 + outlier_step)].index
        outliers_indices.extend(outliers_list_col)
    outliers_indices = list(set(outliers_indices))
    return outliers_indices

def visualize_targeting_results(cobra_model, output_dir, output_viz_dir):
    model_candidate_reactions = []
    for each_reaction in cobra_model.reactions:
        compartments = []
        for each_met in each_reaction.reactants + each_reaction.products:
            compartments.append(each_met.compartment)
        compartments = list(set(compartments))
        if len(compartments) == 1:
            model_candidate_reactions.append(each_reaction.id)
        
    de_df = pd.read_csv(output_dir+'/Differentially_changed_fluxes.csv', index_col=0)
    de_df = de_df[de_df['P']<0.05]
    
    target_reactions = de_df.index
    target_reactions = list(set(model_candidate_reactions) & set(target_reactions))
    
    # Load the data
    file_path_1 = output_dir+'/flux1.csv'
    file_path_2 = output_dir+'/flux2.csv'

    df1 = pd.read_csv(file_path_1, index_col=0)
    df2 = pd.read_csv(file_path_2, index_col=0)
    
    # Merge the two dataframes
    df1 = df1.loc[target_reactions]
    df2 = df2.loc[target_reactions]
    
    merged_df = pd.merge(df1, df2, left_index=True, right_index=True)
    
    # Extract the 'C' and 'H' columns
    c_cols = [col for col in df1.columns]
    t_cols = [col for col in df2.columns]
    
    c_data = merged_df[c_cols].T#.values
    t_data = merged_df[t_cols].T#.values
    merged_df = merged_df.T
    
    # Combine the 'C' and 'H' data
    combined_data = np.concatenate((c_data, t_data), axis=0)
    labels = np.concatenate((np.zeros(c_data.shape[0]), np.ones(t_data.shape[0])))
    
    # Use UMAP to reduce dimensionality
    umap_model = umap.UMAP()
    umap_model.fit(merged_df)
    umap_result = umap_model.transform(merged_df)
    umap_result_df = pd.DataFrame(umap_result, index=merged_df.index)
    umap_result_df['label'] = labels
    umap_result_df.to_csv(output_viz_dir + '/UMAP_Visualization_raw_data.csv')
    
    label_0_data = umap_result_df[umap_result_df['label'] == 0]
    outliers_indices_label_0 = find_outliers_iqr(label_0_data.iloc[:, 1:3])
    outliers_data_label_0 = label_0_data.loc[outliers_indices_label_0]
    non_outliers_data_label_0 = label_0_data.drop(outliers_indices_label_0)
#     center_x = non_outliers_data_label_0.iloc[:, 1].mean()
#     center_y = non_outliers_data_label_0.iloc[:, 2].mean()
    
    # Plot the results
    plt.figure(figsize=(10, 7))
    plt.scatter(umap_result[labels == 0, 0], umap_result[labels == 0, 1], color='skyblue', label='Control', alpha=0.5)
    plt.scatter(umap_result[labels == 1, 0], umap_result[labels == 1, 1], color='lightcoral', label='Test', alpha=0.5)
    plt.title('UMAP Visualization')
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    
    # Save the plot as an SVG file
    plt.savefig(output_viz_dir + '/UMAP_Visualization.svg', format='svg')
    plt.savefig(output_viz_dir + '/UMAP_Visualization.png', format='png')
    plt.show()

    # Save the trained model to a file
    joblib.dump(umap_model, output_viz_dir+'/umap_model.pkl')
    loaded_umap_model = joblib.load(output_viz_dir+'/umap_model.pkl')
    
    center_x = np.median(umap_result_df[umap_result_df['label']==0][0].values)
    center_y = np.median(umap_result_df[umap_result_df['label']==0][1].values)

    targeting_result_files = glob.glob(output_dir+'/targeting_results/*.csv')
    fp = open(output_viz_dir+'/target_genes.csv', 'w')
    fp.write('%s,%s\n' % ('Gene', 'Sample'))

    targeting_result_files = glob.glob(output_dir+'/targeting_results/*.csv')
    fp = open(output_viz_dir+'/target_genes.csv', 'w')
    fp.write('%s,%s\n' % ('Gene', 'Sample'))

    texts = []  # Initialize a list to hold the text objects for adjust_text

    for each_result in targeting_result_files:
        basename = os.path.basename(each_result).split('.')[0].strip()        
        df2 = pd.read_csv(each_result, index_col=0)
        df2 = df2.loc[target_reactions].T
        umap_result_sim = loaded_umap_model.transform(df2)
        sample_name = basename.replace('MOMA_target_results_','')
        
        sample_data = umap_result_df.loc[[sample_name]]
        sample_x = sample_data[0]
        sample_y = sample_data[1]
        
        tmp_df = df2
        tmp_result = loaded_umap_model.transform(tmp_df)
        tmp_df = pd.DataFrame(tmp_result, columns=['X', 'Y'], index=tmp_df.index)
        
#         # Save UMAP results to a CSV file
#         umap_csv_filename = output_viz_dir + '/UMAP_Results_' + basename + '.csv'
#         tmp_df.to_csv(umap_csv_filename)

        distances = np.sqrt((tmp_df['X'] - center_x)**2 + (tmp_df['Y'] - center_y)**2)
        closest_indices = np.argsort(distances)[:10]

        # 기존 코드는 생략하고 직접적으로 수정이 필요한 부분만 나타냈습니다.
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.scatter(umap_result[labels == 0, 0], umap_result[labels == 0, 1], color='skyblue', alpha=0.5)

        ax.scatter(sample_x, sample_y, color='green', marker='+')
        ax.scatter(center_x, center_y, color='blue', marker='+')
        ax.scatter(tmp_result[:, 0], tmp_result[:, 1], color='gray', alpha=0.05)
        ax.scatter(tmp_result[closest_indices, 0], tmp_result[closest_indices, 1], color='red', alpha=0.5, s=50)

        texts = []  # 텍스트 객체를 저장할 리스트를 초기화합니다.
        for i in closest_indices:
            texts.append(ax.text(tmp_result[i, 0], tmp_result[i, 1], tmp_df.index[i], fontsize=9))

        ax.set_xlabel('UMAP_1')
        ax.set_ylabel('UMAP_2')
        ax.set_title('Simulation Results')
        # ax.legend()

        # 그래프를 저장하기 전에 그림을 그립니다. 이것은 모든 텍스트 객체들이 그림과 연결되어 있음을 보장합니다.
        plt.draw()

        # adjust_text 함수를 호출할 때, 현재 축(ax) 인자를 전달합니다.
        adjust_text(texts, ax=ax)

        plt.savefig(output_viz_dir + '/UMAP_%s.svg' % basename, format='svg')
        plt.savefig(output_viz_dir + '/UMAP_%s.png' % basename, format='png')
        plt.close(fig)  # 사용이 끝난 후 그림을 닫습니다.

        target_genes = tmp_df.index[closest_indices].values
        for gene in target_genes:
            fp.write('%s,%s\n' % (gene, basename))
    fp.close()
    return

def visualize_other_results(output_dir, output_viz_dir):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    file_path = output_dir+'/Up_regulated_pathways.csv'
    df = pd.read_csv(file_path, index_col=0)
    df['-log10(FDR corrected P)'] = -np.log10(df['FDR corrected P'])
    df_sorted = df.sort_values('-log10(FDR corrected P)', ascending=False)
    plt.rcParams['font.size'] = 12
    plt.style.use('ggplot')
    plt.figure(figsize=(10, 6))
    sns.barplot(x=df_sorted.index, y='-log10(FDR corrected P)', data=df_sorted, palette='viridis')
    plt.xticks(rotation=90)
    plt.xlabel('')
    plt.ylabel('-log10(FDR corrected P)')
    plt.title('Up regulated pathways')
    plt.tight_layout()
    plt.savefig(output_viz_dir + '/Enrich_up_regulated_pathways.svg', format='svg')
    plt.savefig(output_viz_dir + '/Enrich_up_regulated_pathways.png', format='png')
    plt.show()
    
    file_path = output_dir+'/Down_regulated_pathways.csv'
    df = pd.read_csv(file_path, index_col=0)
    df['-log10(FDR corrected P)'] = -np.log10(df['FDR corrected P'])
    df_sorted = df.sort_values('-log10(FDR corrected P)', ascending=False)
    plt.rcParams['font.size'] = 12
    plt.style.use('ggplot')
    plt.figure(figsize=(10, 6))
    sns.barplot(x=df_sorted.index, y='-log10(FDR corrected P)', data=df_sorted, palette='viridis')
    plt.xticks(rotation=90)
    plt.xlabel('')
    plt.ylabel('-log10(FDR corrected P)')
    plt.title('Down regulated pathways')
    plt.tight_layout()
    plt.savefig(output_viz_dir + '/Enrich_down_regulated_pathways.svg', format='svg')
    plt.savefig(output_viz_dir + '/Enrich_down_regulated_pathways.png', format='png')
    plt.show()
    
    file_path = output_dir+'/condition1/functional_model_information.csv'
    df = pd.read_csv(file_path, index_col=0)
    df = df[['No. of genes', 'No. of metabolites', 'No. of reactions']]
    plt.style.use('ggplot')
    fig, ax = plt.subplots(figsize=(10, 6))
    df.plot(kind='bar', ax=ax, color=['skyblue', 'lightgreen', 'lightcoral'])
    plt.xticks(rotation=0)
    plt.xlabel('')
    plt.ylabel('Count')
    plt.title('Model Information (control)')

    for p in ax.patches:
        ax.annotate(str(int(p.get_height())), (p.get_x() * 1.005, p.get_height() * 1.005), color='black')
    plt.tight_layout()
    plt.savefig(output_viz_dir + '/Condition1_model_info.svg', format='svg')
    plt.savefig(output_viz_dir + '/Condition1_model_info.png', format='png')
    plt.show()
    
    file_path = output_dir+'/condition2/functional_model_information.csv'
    df = pd.read_csv(file_path, index_col=0)
    df = df[['No. of genes', 'No. of metabolites', 'No. of reactions']]
    plt.style.use('ggplot')

    fig, ax = plt.subplots(figsize=(10, 6))
    df.plot(kind='bar', ax=ax, color=['skyblue', 'lightgreen', 'lightcoral'])
    plt.xticks(rotation=0)
    plt.xlabel('')
    plt.ylabel('Count')
    plt.title('Model Information (test)')

    for p in ax.patches:
        ax.annotate(str(int(p.get_height())), (p.get_x() * 1.005, p.get_height() * 1.005), color='black')
    plt.tight_layout()
    plt.savefig(output_viz_dir + '/Condition2_model_info.svg', format='svg')
    plt.savefig(output_viz_dir + '/Condition2_model_info.png', format='png')
    plt.show()
    
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
    output_viz_dir = '%s/viz/'%(output_dir)
    if not os.path.isdir(output_viz_dir):
        os.mkdir(output_viz_dir)

    targeting_result_dir = output_dir + '/targeting_results/'
    if not os.path.isdir(targeting_result_dir):
        os.mkdir(targeting_result_dir)
        
    ## Load data    
    metabolic_task_file = pkg_resources.resource_filename('tf_metabolism', 'data/MetabolicTasks.csv')
    medium_file = pkg_resources.resource_filename('tf_metabolism', 'data/RPMI1640_medium.txt')
    present_metabolite_file = pkg_resources.resource_filename('tf_metabolism', 'data/essential_metabolites.txt')
    essential_reaction_file = pkg_resources.resource_filename('tf_metabolism', 'data/essential_reactions.txt')
    trrust = pkg_resources.resource_filename('tf_metabolism', 'data/TRRUST_v2_ensembl.tsv')
    
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
    
    target_genes_list1 = []
    target_genes_list2 = []
    
    for each_row, each_df in omics1_df.iterrows():
        min_exp = np.min(np.abs(each_df.values))
        if min_exp >= 0.5:
            target_genes_list1.append(int(each_row))
            
    for each_row, each_df in omics2_df.iterrows():
        min_exp = np.min(np.abs(each_df.values))
        if min_exp >= 0.5:
            target_genes_list2.append(int(each_row))
    
    target_genes2 = list(set(target_genes_list1) & set(target_genes_list2))
    omics1_df_tmp = omics1_df.loc[target_genes2]
    omics2_df_tmp = omics2_df.loc[target_genes2]

    comparison_result_df = statistical_comparison.two_grouped_data_comparison(omics1_df_tmp, omics2_df_tmp, related_sample_flag, 0.05)
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
    
    run_targeting_simulation(output_dir, targeting_result_dir)
    visualize_targeting_results(cobra_model, output_dir, output_viz_dir)
    visualize_other_results(output_dir, output_viz_dir)
    
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
    
if __name__ == '__main__':
    main()
