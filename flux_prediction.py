from tf_metabolism.metabolic_simulation import flux_prediction
from tf_metabolism.metabolic_simulation import Simulator

from tf_metabolism.omics_integration import LAD
from cobra.io import read_sbml_model

import pandas as pd
import os

def calculate_flux(output_dir, omics_file, model_file, minimum_biomass=0.01):    
    cobra_model = read_sbml_model(model_file)
    
    basename = os.path.basename(omics_file).split('.')[0].strip()
    expression_score = flux_prediction.read_expressoin_data(omics_file)

    reaction_weight_info = {}
    reaction_weights = flux_prediction.calculate_reaction_score(cobra_model, expression_score)
    reaction_weight_info[basename] = reaction_weights
    
    obj = LAD.LAD()
    obj.load_cobra_model(cobra_model)
    flux_constraints = {}
    flux_constraints['biomass_reaction'] = [minimum_biomass, 1000.0]
    
    model_reactions = [each_reaction.id for each_reaction in cobra_model.reactions]
    new_reaction_weight_info = {}
    for each_reaction in reaction_weights:
        if each_reaction in model_reactions:
            new_reaction_weight_info[each_reaction] = reaction_weights[each_reaction]/10000.0
    solution_status, objective_value, predicted_flux = obj.run_LP_fitting(opt_flux=new_reaction_weight_info, flux_constraints=flux_constraints)
    return predicted_flux


if __name__ == '__main__':
    omics_file = './input_data2/CAF-62T_TPM.csv'
    model_file = './input_data2/functional_CAF-62T_TPM.xml'
    output_dir = './output_test/'
    
    # omics 데이터로부터 flux를 예측하는 방법
    flux_dist = calculate_flux(output_dir, omics_file, model_file, minimum_biomass=0.1)  
    
    # 이미 계산된 flux를 읽는 방법
#     flux_dist2 = {}
#     flux_file = './input_data2/Flux_prediction_CAF-62T_TPM.csv'
#     with open(flux_file, 'r') as fp:
#         for line in fp:
#             sptlist = line.strip().split(',')
#             reaction = sptlist[0].strip()
#             flux = float(sptlist[1].strip())
#             flux_dist2[reaction] = flux
    
    # flux를 기반으로 perturbation 시뮬레이션 예
    cobra_model = read_sbml_model(model_file)
    
    obj = Simulator.Simulator()
    obj.load_cobra_model(cobra_model)
    
    for each_gene in cobra_model.genes:
        flux_constraints = {}
        target_reactions = []
        for each_reaction in each_gene.reactions:
            target_reactions.append(each_reaction.id)
            flux_constraints[each_reaction.id] = [0.0, 0.0]        
        
        status, obj_value, perturbed_flux_dist = obj.run_MOMA(wild_flux=flux_dist, flux_constraints=flux_constraints)
        print(perturbed_flux_dist)
        print(each_gene, target_reactions, perturbed_flux_dist['biomass_reaction'])
        break
        
        
        