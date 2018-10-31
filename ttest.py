#!/usr/bin/python
import pandas as pd

def main():
    # read files to be used
    file_training_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = 10)
    file_training_output = pd.read_csv("mapping_sample_to_class.csv", usecols = ['GEO asscession number', 'relapses within 5 years (1 = yes, 0=no)'])
    
    # separate data into 2 classes
    sample_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['1'])]
    sample_no_relapse = file_training_output.loc[file_training_output['relapses within 5 years (1 = yes, 0=no)'].isin(['0'])]
    print(sample_no_relapse)
    
    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)
    print(list_sample_relapse)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)
    # print(list_sample_no_relapse)
    # print(sample_no_relapse.loc[sample_no_relapse['GEO asscession number']])
if __name__ == '__main__':
    main()