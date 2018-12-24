#!/usr/bin/python
import pandas as pd
from scipy import stats
import math

def main():
    # read files to be used
    # number of all genes in GSE2034-22071 is 22283 genes
    row_to_read = 22283
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output = pd.read_csv("mapping_sample_to_class.csv", usecols = ['GEO asscession number', 'relapses within 5 years (1 = yes, 0=no)'])
    
    # get number of ranked gene to be shown
    while True:
        number_of_ranked_gene = input("Number of ranked feature: ")
        if ((number_of_ranked_gene.isnumeric() == False) or (int(number_of_ranked_gene) > row_to_read) or (int(number_of_ranked_gene) <= 0)):
            print("Invalid input...")
        else:
            break

    # separate data into 2 classes
    sample_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['1'])]
    sample_no_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['0'])]
    # print(sample_no_relapse)
    
    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)
    # print(list_sample_relapse)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_training_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)
    # print(list_gene_name)
        

    # get gene expression of all sample in the same class
    # print(file_training_input.loc[:, list_sample_relapse])
    list_gene_exp_relapse = []
    for i in range(0, row_to_read):
        gene_exp_relapse = []
        for column in  file_training_input.loc[i, list_sample_relapse]:
            gene_exp_relapse.append(column)
        list_gene_exp_relapse.append(gene_exp_relapse)
    # print(list_gene_exp_relapse)

    # gene expression for each gene from samples with no relapse within 5 years
    # print(file_training_input.loc[:, list_sample_no_relapse])
    list_gene_exp_no_relapse = []
    for i in range(0, row_to_read):
        gene_exp_no_relapse = []
        for column in  file_training_input.loc[i, list_sample_no_relapse]:
            gene_exp_no_relapse.append(column)
        list_gene_exp_no_relapse.append(gene_exp_no_relapse)
    # print(list_gene_exp_no_relapse)

    # conducting t-test
    ttest_result = []
    for i in range(0, row_to_read):      
        score = []
        # print(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False))
        # get absolute magnitude of t-test value
        abs_ttest_value = math.fabs(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[0])
        p_value = stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[1]
        # add element with this format (gene_order_id, ttest_value)
        score.append(i)
        score.append(abs_ttest_value)
        ttest_result.append(score)
    # ranking elements using their t-test value in descending order
    ttest_result.sort(key=lambda x: x[1], reverse=True)
    # print(ttest_result)

    # create list of ranked gene
    ranked_gene = []
    for i in range(0, len(ttest_result)):
        gene_order_id = ttest_result[i][0]
        # print(gene_order_id)
        # print(list_gene_name[gene_order_id][1])
        ranked_gene.append(list_gene_name[gene_order_id][1])
    # print(ranked_gene)

    # show top ranked feature
    for i in range(0, int(number_of_ranked_gene)):
        print(ranked_gene[i] + " => " + "t-test value : " + str(ttest_result[i][1]))

if __name__ == '__main__':
    main()