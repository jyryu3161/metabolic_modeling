import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection


def pathway_enrichment_analysis(cobra_model, selected_reactions):
    pathway_enrichment_result_info = {}
    all_reactions = [each_reaction.id for each_reaction in cobra_model.reactions]

    pathway_info = {}
    for each_reaction in cobra_model.reactions:
        subsystem = each_reaction.subsystem
        if subsystem not in pathway_info:
            pathway_info[subsystem] = [each_reaction.id]
        else:
            pathway_info[subsystem].append(each_reaction.id)

    p_value_list = []
    pathway_list = []
    odd_ratio_list = []
    data_information = {}
    for each_pathway in pathway_info:
        pathway_reactions = pathway_info[each_pathway]
        no_pathway_reactions = set(all_reactions).difference(pathway_reactions)
        not_selected_reactions = set(all_reactions).difference(set(selected_reactions))

        not_selected_reactions_no_pathway = set(no_pathway_reactions) & set(not_selected_reactions)
        not_selected_reactions_pathway = set(pathway_reactions) & set(not_selected_reactions)

        selected_reactions_no_pathway = set(no_pathway_reactions) & set(selected_reactions)
        selected_reactions_pathway = set(pathway_reactions) & set(selected_reactions)

        oddsratio, p_value = stats.fisher_exact(
            [[len(not_selected_reactions_no_pathway), len(not_selected_reactions_pathway)],
             [len(selected_reactions_no_pathway), len(selected_reactions_pathway)]])
        p_value_list.append(p_value)
        odd_ratio_list.append(oddsratio)
        pathway_list.append(each_pathway)
        data_information[each_pathway] = {}
        data_information[each_pathway]['No. of selected reactions in pathway'] = len(selected_reactions_pathway)
        data_information[each_pathway]['No. of all reactions in pathway'] = len(pathway_reactions)
    rejected, p_value_corrected = fdrcorrection(p_value_list)
    for i in range(len(p_value_corrected)):
        pathway = pathway_list[i]
        pathway_enrichment_result_info[pathway] = {}
        pathway_enrichment_result_info[pathway]['FDR corrected P'] = p_value_corrected[i]
        pathway_enrichment_result_info[pathway]['P'] = p_value_list[i]
        pathway_enrichment_result_info[pathway]['ODDS RATIO'] = odd_ratio_list[i]
        pathway_enrichment_result_info[pathway]['No. of selected reactions in pathway'] = data_information[pathway][
            'No. of selected reactions in pathway']
        pathway_enrichment_result_info[pathway]['No. of all reactions in pathway'] = data_information[pathway][
            'No. of all reactions in pathway']
    pathway_enrichment_result_df = pd.DataFrame.from_dict(pathway_enrichment_result_info)
    return pathway_enrichment_result_df.T


def tf_enrichment_analysis(tf_info, selected_genes):
    tf_enrichment_result_info = {}
    
    all_genes = []
    for tf in tf_info:
        genes = tf_info[tf]
        all_genes += genes
    all_genes = list(set(all_genes))
    
    p_value_list = []
    tf_list = []
    odd_ratio_list = []
    data_information = {}
    for each_tf in tf_info:
        tf_genes = tf_info[each_tf]
        no_tf_genes = set(all_genes).difference(tf_genes)
        not_selected_genes = set(all_genes).difference(set(selected_genes))
        
        not_selected_genes_no_tf = set(no_tf_genes) & set(not_selected_genes)
        not_selected_genes_tf = set(tf_genes) & set(not_selected_genes)
        
        selected_genes_no_tf = set(no_tf_genes) & set(selected_genes)
        selected_genes_tf = set(tf_genes) & set(selected_genes)
        
        oddsratio, p_value = stats.fisher_exact(
            [[len(not_selected_genes_no_tf), len(not_selected_genes_tf)],
             [len(selected_genes_no_tf), len(selected_genes_tf)]])
        
        p_value_list.append(p_value)
        odd_ratio_list.append(oddsratio)
        tf_list.append(each_tf)
        
        data_information[each_tf] = {}
        data_information[each_tf]['No. of selected genes in tf'] = len(selected_genes_tf)
        data_information[each_tf]['No. of all genes in tf'] = len(tf_genes)
        
    rejected, p_value_corrected = fdrcorrection(p_value_list)
    for i in range(len(p_value_corrected)):
        tf = tf_list[i]
        tf_enrichment_result_info[tf] = {}
        tf_enrichment_result_info[tf]['FDR corrected P'] = p_value_corrected[i]
        tf_enrichment_result_info[tf]['P'] = p_value_list[i]
        tf_enrichment_result_info[tf]['No. of selected genes in tf'] = data_information[tf][
            'No. of selected genes in tf']
        tf_enrichment_result_info[tf]['No. of all genes in tf'] = data_information[tf][
            'No. of all genes in tf']
    tf_enrichment_result_df = pd.DataFrame.from_dict(tf_enrichment_result_info)
    return tf_enrichment_result_df.T
