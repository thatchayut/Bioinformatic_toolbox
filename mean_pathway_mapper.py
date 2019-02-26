from scipy import stats
import pandas as pd
import random
import math
import calculate
import time
import xlsxwriter
from copy import deepcopy
from sklearn.metrics import roc_auc_score

def main(): 
    # prepare data
    # default row_to_read = 22283
    row_to_read_file_input = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read_file_input)
    file_training_output= pd.read_csv("mapping_sample_to_class_full.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway = 1329
    rows_to_read_file_pathway = 1329
    file_ref_name = "accession_number_to_entrez_id.csv"
    file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read_file_input):
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

    # consider non-relapse and relapse (not in specific period of time)
    sample_relapse = file_training_output.loc[file_training_output['relapse (1=True)'].isin(['1'])]
    sample_no_relapse = file_training_output.loc[file_training_output['relapse (1=True)'].isin(['0'])]

    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)

    # get output file's name
    file_name = input("Name of output file : ")

    # initial file to be written as an output
    workbook = xlsxwriter.Workbook(file_name + ".xlsx")
    worksheet = workbook.add_worksheet()

    print("Process : Creating collection to collect samples and their genes' expression")
    # create dictionary used to collect pathways of each sample
    samples_relapse = {}
    samples_no_relapse = {}

    # get all pathways of all samples in class 'relapse'
    for element_index in range(0, len(list_sample_relapse)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " relapse is in progress ...")
        print(str(len(list_sample_relapse) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_relapse[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

        sample.append(sample_name)
        sample.append(pathways)
        samples_relapse[element_index] = sample

    for element_index in range(0, len(list_sample_no_relapse)):
        print()
        print("Creating pathways for sample " + str(element_index + 1) + " non-relapse is in progress ...")
        print(str(len(list_sample_no_relapse) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_no_relapse[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, normalize = False)

        sample.append(sample_name)
        sample.append(pathways)
        samples_no_relapse[element_index] = sample
    
    print("Process : Creating collections of samples with their pathways' activity ...")
    # create collections of samples with their pathways
    # data will be collected in this format
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
    samples_relapse_pathway_activity = {}
    samples_no_relapse_pathway_activity = {}
    for samples_index in range(0, len(samples_relapse)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_relapse[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_relapse[samples_index][1][pathway_index][1])):
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
        samples_relapse_pathway_activity[samples_index] = sample

    for samples_index in range(0, len(samples_no_relapse)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_no_relapse[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_no_relapse[samples_index][1][pathway_index][1])):
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
        samples_no_relapse_pathway_activity[samples_index] = sample

    # merge both classes together
    samples_all_pathway_activity = {}
    for sample_index in range(0, len(samples_relapse_pathway_activity)):
        samples_all_pathway_activity[sample_index] = samples_relapse_pathway_activity[sample_index]
    
    next_index = len(samples_all_pathway_activity)

    for sample_index in range(0, len(samples_no_relapse_pathway_activity)):
        index_to_add = next_index + sample_index

        samples_all_pathway_activity[index_to_add] = samples_no_relapse_pathway_activity[sample_index]

    # write to file
    print("Process : Write to an output file ...")
    worksheet.write(0, 0, "PATHWAY_NAME")
    
    # write pathways' name
    for pathway_index in range(0, len(list_pathway_name)):
        pathway_name = list_pathway_name[pathway_index][1]
        worksheet.write(pathway_index + 1, 0, pathway_name)

    # write samples' name
    for sample_index in range(0, len(samples_all_pathway_activity)):
        sample_name = samples_all_pathway_activity[sample_index][0]

        worksheet.write(0, sample_index + 1, sample_name)
    
    # write pathways' activity of each sample
    for sample_index in range(0, len(samples_all_pathway_activity)):
        for pathway_index in range(0, len(samples_all_pathway_activity[sample_index][1])):
            pathway_activity = samples_all_pathway_activity[sample_index][1][pathway_index][1]

            row_to_write = pathway_index + 1
            col_to_write = sample_index + 1

            worksheet.write(row_to_write, col_to_write, pathway_activity)


if __name__ == "__main__":
    main()