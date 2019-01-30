from scipy import stats
import pandas as pd
import random
import math
import calculate
import time
from copy import deepcopy
from sklearn.metrics import roc_auc_score

def main():

    method_id = "1"

    # prepare used file
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output= pd.read_csv("mapping_sample_to_class_full.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    rows_to_read_file_pathway = 1329
    file_ref_name = "accession_number_to_entrez_id.csv"
    file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # list to collect sample name
    list_sample_train_relapse_name = ['GSM37011', 'GSM36931', 'GSM36881', 'GSM36941', 'GSM36952', 'GSM37004', 'GSM37029', 'GSM36835', 'GSM36902', 'GSM36996', 'GSM36813', 'GSM36858', 'GSM36955', 'GSM36784', 'GSM37053', 'GSM37052', 'GSM36860', 'GSM36967', 'GSM36939', 'GSM36792', 'GSM36954', 'GSM37035', 'GSM36918', 'GSM36885', 'GSM36976', 'GSM37020', 'GSM36898', 'GSM37038', 'GSM37006', 'GSM36870', 'GSM36950', 'GSM37022', 'GSM36986', 'GSM36826', 'GSM37027', 'GSM36838', 'GSM36927', 'GSM36960', 'GSM36800', 'GSM37005', 'GSM37008', 'GSM37028', 'GSM37040', 'GSM36957', 'GSM36872', 'GSM36920', 'GSM36956', 'GSM36908', 'GSM36811', 'GSM36972', 'GSM36923', 'GSM36911', 'GSM37049', 'GSM37026', 'GSM36969', 'GSM36897', 'GSM36985', 'GSM36874', 'GSM36926', 'GSM36877', 'GSM36875', 'GSM36974', 'GSM36949', 'GSM36862', 'GSM36971', 'GSM36943', 'GSM36989', 'GSM37031', 'GSM37036', 'GSM37002', 'GSM36797', 'GSM37037', 'GSM36815', 'GSM37003', 'GSM37018', 'GSM37007', 'GSM37039', 'GSM37030', 'GSM37051', 'GSM36839', 'GSM37023', 'GSM36973', 'GSM36999', 'GSM36924', 'GSM37041', 'GSM36903', 'GSM36937', 'GSM36818', 'GSM36946', 'GSM37050', 'GSM36778', 'GSM37058', 'GSM36998', 'GSM36789', 'GSM36879', 'GSM37013']
    list_sample_train_no_relapse_name = ['GSM37033', 'GSM36840', 'GSM36803', 'GSM36915', 'GSM36861', 'GSM36798', 'GSM36831', 'GSM36808', 'GSM36799', 'GSM36896', 'GSM36992', 'GSM36873', 'GSM36794', 'GSM36970', 'GSM36906', 'GSM36940', 'GSM36790', 'GSM37060', 'GSM36833', 'GSM37046', 'GSM36938', 'GSM36930', 'GSM36810', 'GSM36787', 'GSM36979', 'GSM36823', 'GSM37059', 'GSM36953', 'GSM36850', 'GSM36851', 'GSM36780', 'GSM36824', 'GSM36995', 'GSM36962', 'GSM36977', 'GSM36910', 'GSM36944', 'GSM37019', 'GSM36889', 'GSM36901', 'GSM36917', 'GSM36865', 'GSM36987', 'GSM37009', 'GSM36804', 'GSM37015', 'GSM36869', 'GSM36981', 'GSM36864', 'GSM37017', 'GSM36959', 'GSM36966', 'GSM36867', 'GSM36842', 'GSM36922', 'GSM36834', 'GSM37061', 'GSM37047', 'GSM36822', 'GSM36963', 'GSM37016', 'GSM37021', 'GSM36980', 'GSM36855', 'GSM36827', 'GSM36791', 'GSM36900', 'GSM36887', 'GSM36891', 'GSM36866', 'GSM36781', 'GSM36856', 'GSM36868', 'GSM36821', 'GSM36845', 'GSM36849', 'GSM36837', 'GSM36848', 'GSM36894', 'GSM37055', 'GSM36929', 'GSM36847', 'GSM37057', 'GSM37044', 'GSM36786', 'GSM37056', 'GSM36991', 'GSM36945', 'GSM36828', 'GSM36916', 'GSM36893', 'GSM36909', 'GSM36802', 'GSM36975', 'GSM36832', 'GSM36825', 'GSM37048', 'GSM36884', 'GSM36859', 'GSM36853', 'GSM36936', 'GSM36793', 'GSM36948', 'GSM36904', 'GSM37014', 'GSM36912', 'GSM36779', 'GSM36841', 'GSM36777', 'GSM36836', 'GSM36914', 'GSM36978', 'GSM36830', 'GSM37000', 'GSM36934', 'GSM36857', 'GSM36958', 'GSM36961', 'GSM36878', 'GSM36816', 'GSM36796', 'GSM36820', 'GSM36882', 'GSM36852', 'GSM36795', 'GSM36876', 'GSM36801', 'GSM36846', 'GSM37062', 'GSM36907', 'GSM36883', 'GSM36817', 'GSM36805', 'GSM36783', 'GSM36984', 'GSM37012', 'GSM36925', 'GSM36785', 'GSM36806', 'GSM37045', 'GSM36880', 'GSM36886', 'GSM36988', 'GSM36890', 'GSM36935', 'GSM36892', 'GSM36993', 'GSM36951', 'GSM36982', 'GSM36844', 'GSM37043', 'GSM36871', 'GSM37054', 'GSM36899', 'GSM36990', 'GSM36863', 'GSM36809', 'GSM36807', 'GSM36921', 'GSM36843', 'GSM37034']
    list_samples_test_relapse_name = ['GSM36888', 'GSM37042', 'GSM36983', 'GSM36994', 'GSM37001', 'GSM36905', 'GSM36928', 'GSM36947', 'GSM36997', 'GSM36964', 'GSM36814']
    list_samples_test_no_relapse_name = ['GSM36782', 'GSM36812', 'GSM36829', 'GSM36919', 'GSM36933', 'GSM37010', 'GSM36968', 'GSM36913', 'GSM36942', 'GSM37025', 'GSM36895', 'GSM37024', 'GSM36965', 'GSM36788', 'GSM37032', 'GSM36819', 'GSM36854', 'GSM36932']
    
    feature_set = ['REACTOME_METABOLISM_OF_CARBOHYDRATES', 'REACTOME_DNA_STRAND_ELONGATION', 'BIOCARTA_EPHA4_PATHWAY', 'NABA_COLLAGENS', 'ST_ERK1_ERK2_MAPK_PATHWAY', 'REACTOME_APC_CDC20_MEDIATED_DEGRADATION_OF_NEK2A', 'REACTOME_S_PHASE', 'REACTOME_REMOVAL_OF_THE_FLAP_INTERMEDIATE_FROM_THE_C_STRAND', 'REACTOME_FORMATION_OF_THE_HIV1_EARLY_ELONGATION_COMPLEX', 'REACTOME_EXTENSION_OF_TELOMERES', 'REACTOME_PHOSPHORYLATION_OF_THE_APC_C', 'REACTOME_UNWINDING_OF_DNA', 'REACTOME_RAF_MAP_KINASE_CASCADE', 'REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS', 'SIG_REGULATION_OF_THE_ACTIN_CYTOSKELETON_BY_RHO_GTPASES']

    # create list of all samples in trainig data in this fold
    list_train_all_samples = []
    list_train_all_samples.extend(list_sample_train_relapse_name)
    list_train_all_samples.extend(list_sample_train_no_relapse_name)

    list_test_all_samples = []
    list_test_all_samples.extend(list_samples_test_relapse_name)
    list_test_all_samples.extend(list_samples_test_no_relapse_name)

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_training_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)
    
    # get list of pathway name
    list_pathway_name = []
    for i in range(0, rows_to_read_file_pathway):
        pathway_name = []
        pathway_name.append(i)
        pathway_name.append(file_pathway.loc[i, "PATHWAY_NAME"])
        list_pathway_name.append(pathway_name)
    
    # list all gene expression to calculate mean and sd
    print("\n Gathering gene expressions are in progress ... ")
    list_train_all_gene_expression = []
    for line_index in range(0, row_to_read):
        for column in file_training_input.loc[line_index, list_train_all_samples]:
            list_train_all_gene_expression.append(column)

    mean_all_gene_expression_train = calculate.mean(list_train_all_gene_expression)
    sd_all_gene_expression_train = calculate.sd(list_train_all_gene_expression)
    max_all_gene_expression_train = max(list_train_all_gene_expression)
    min_all_gene_expression_train = min(list_train_all_gene_expression)

    # create samples
    samples_relapse = {}
    samples_no_relapse = {}
    samples_all = {}

    # class 'relapse'
    for element_index in range(0, len(list_sample_train_relapse_name)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
        print(str(len(list_sample_train_relapse_name) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_train_relapse_name[element_index]
        # pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

        if (method_id is "1"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression_train, sd_of_data = sd_all_gene_expression_train, \
                        method = "z_score")
        elif (method_id is "2"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                        method = "narrow_scaling")            
        elif (method_id is "3"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                        method = "wide_scaling")    

        sample.append(sample_name)
        sample.append(pathways)
        samples_relapse[element_index] = sample

    for element_index in range(0, len(list_sample_train_no_relapse_name)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " no-relapse is in progress ...")
        print(str(len(list_sample_train_no_relapse_name) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_train_no_relapse_name[element_index]
        # pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

        if (method_id is "1"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression_train, sd_of_data = sd_all_gene_expression_train, \
                        method = "z_score")
        elif (method_id is "2"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                        method = "narrow_scaling")            
        elif (method_id is "3"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                        method = "wide_scaling")    

        sample.append(sample_name)
        sample.append(pathways)
        samples_no_relapse[element_index] = sample
    
    for element_index in range(0, len(list_test_all_samples)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " - all is in progress ...")
        print(str(len(list_test_all_samples) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_test_all_samples[element_index]
        # pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway)

        if (method_id is "1"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression_train, sd_of_data = sd_all_gene_expression_train, \
                        method = "z_score")
        elif (method_id is "2"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                        method = "narrow_scaling")            
        elif (method_id is "3"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression_train, min_of_data = min_all_gene_expression_train, \
                        method = "wide_scaling")    

        sample.append(sample_name)
        sample.append(pathways)
        samples_all[element_index] = sample
    
    # calculate activity score
    samples_training_relapse_pathway_activity = {}
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
    for samples_index in range(0, len(samples_relapse)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_relapse[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_relapse[samples_index][1][pathway_index][1])):
                # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                # print(gene_index)
                gene_expression = samples_relapse[samples_index][1][pathway_index][1][gene_index][1]
                list_gene_expression_in_pathway.append(gene_expression)

            # data to collect as pathway activity
            pathway_name = samples_relapse[samples_index][1][pathway_index][0]
            pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            pathway.append(pathway_name)
            pathway.append(pathway_activity)
            list_pathway.append(pathway)
        
        sample_name = samples_relapse[samples_index][0]
        sample.append(sample_name)
        sample.append(list_pathway)
        samples_training_relapse_pathway_activity[samples_index] = sample
    
    samples_training_no_relapse_pathway_activity = {}
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
    for samples_index in range(0, len(samples_no_relapse)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_no_relapse[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_no_relapse[samples_index][1][pathway_index][1])):
                # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                # print(gene_index)
                gene_expression = samples_no_relapse[samples_index][1][pathway_index][1][gene_index][1]
                list_gene_expression_in_pathway.append(gene_expression)

            # data to collect as pathway activity
            pathway_name = samples_no_relapse[samples_index][1][pathway_index][0]
            pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            pathway.append(pathway_name)
            pathway.append(pathway_activity)
            list_pathway.append(pathway)
        
        sample_name = samples_no_relapse[samples_index][0]
        sample.append(sample_name)
        sample.append(list_pathway)
        samples_training_no_relapse_pathway_activity[samples_index] = sample

    samples_test_all_pathway_activity = {}
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
    for samples_index in range(0, len(samples_all)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_all[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_all[samples_index][1][pathway_index][1])):
                # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                # print(gene_index)
                gene_expression = samples_all[samples_index][1][pathway_index][1][gene_index][1]
                list_gene_expression_in_pathway.append(gene_expression)

            # data to collect as pathway activity
            pathway_name = samples_all[samples_index][1][pathway_index][0]
            pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            pathway.append(pathway_name)
            pathway.append(pathway_activity)
            list_pathway.append(pathway)
        
        sample_name = samples_all[samples_index][0]
        sample.append(sample_name)
        sample.append(list_pathway)
        samples_test_all_pathway_activity[samples_index] = sample
    
    # create list of sample with pahtway activity in relation to the feature set
    list_sample_relapse_pathway_expression_eval = []
    for sample_index in range(0, len(samples_training_relapse_pathway_activity)):
        list_pathway_activity = []
        for pathway_index in range(0, len(samples_training_relapse_pathway_activity[sample_index][1])):
            pathway_name = samples_training_relapse_pathway_activity[sample_index][1][pathway_index][0]
            # print(pathway_name)
            pathway_activity = samples_training_relapse_pathway_activity[sample_index][1][pathway_index][1]

            if (pathway_name in feature_set):
                list_pathway_activity.append(pathway_activity)

        list_sample_relapse_pathway_expression_eval.append(list_pathway_activity)

    list_sample_no_relapse_pathway_expression_eval = []
    for sample_index in range(0, len(samples_training_no_relapse_pathway_activity)):
        list_pathway_activity = []
        for pathway_index in range(0, len(samples_training_no_relapse_pathway_activity[sample_index][1])):
            pathway_name = samples_training_no_relapse_pathway_activity[sample_index][1][pathway_index][0]
            # print(pathway_name)
            pathway_activity = samples_training_no_relapse_pathway_activity[sample_index][1][pathway_index][1]

            if (pathway_name in feature_set):
                list_pathway_activity.append(pathway_activity)

        list_sample_no_relapse_pathway_expression_eval.append(list_pathway_activity)
    
    list_sample_all_pathway_expression_eval = []
    for sample_index in range(0, len(samples_test_all_pathway_activity)):
        list_pathway_activity = []
        for pathway_index in range(0, len(samples_test_all_pathway_activity[sample_index][1])):
            pathway_name = samples_test_all_pathway_activity[sample_index][1][pathway_index][0]
            # print(pathway_name)
            pathway_activity = samples_test_all_pathway_activity[sample_index][1][pathway_index][1]

            if (pathway_name in feature_set):
                list_pathway_activity.append(pathway_activity)

        list_sample_all_pathway_expression_eval.append(list_pathway_activity)

    # create list of desired output sorted according to the order of samples in 'list_test_all_samples'
    file_desired_outputs_eval = file_training_output.loc[file_training_output['GEO asscession number'].isin(list_test_all_samples)]
    file_desired_outputs_eval['pathway_id'] = file_desired_outputs_eval['GEO asscession number'].apply(lambda name: list_test_all_samples.index(name)) 
    file_desired_outputs_eval = file_desired_outputs_eval.sort_values(by = ['pathway_id'])
    file_desired_outputs_eval.drop(columns = 'pathway_id', inplace = True)

    list_desired_outputs_eval = []
    for element in file_desired_outputs_eval.loc[:, 'relapse (1=True)']:
        list_desired_outputs_eval.append(element)
    
    # find output using lda
    list_actual_outputs_eval = calculate.lda(list_sample_all_pathway_expression_eval, list_sample_relapse_pathway_expression_eval, list_sample_no_relapse_pathway_expression_eval)

    # calculate auc score 
    auc_score = roc_auc_score(list_desired_outputs_eval, list_actual_outputs_eval)    

    print(" #### SUMMARY ####")
    print(" Size of feature set : " + str(len(feature_set)))
    print(" Feture set : " )
    print(feature_set)
    print(" Actual outputs : " + str(list_actual_outputs_eval))
    print(" Desired outputs : " + str(list_desired_outputs_eval))
    print(" AUC score : " + str(auc_score))
if __name__ == '__main__':
    main()

   