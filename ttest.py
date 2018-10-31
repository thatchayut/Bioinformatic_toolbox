#!/usr/bin/python
import pandas as pd
from scipy import stats
import math

def main():
    # read files to be used
    row_to_read = 10
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_training_output = pd.read_csv("mapping_sample_to_class.csv", usecols = ['GEO asscession number', 'relapses within 5 years (1 = yes, 0=no)'])
    
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
        print(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False))
        # get absolute magnitude of t-test value
        abs_ttest_value = math.fabs(stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[0])
        p_value = stats.ttest_ind(list_gene_exp_relapse[i], list_gene_exp_no_relapse[i], equal_var = False)[1]
        # add element with this format (gene_order_id, ttest_value)
        score.append(i)
        score.append(abs_ttest_value)
        ttest_result.append(score)
    print(ttest_result)
        
if __name__ == '__main__':
    main()