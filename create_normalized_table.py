import xlwt
import pandas as pd
import calculate

def main():
    # prepare data
    # row_to_read = 22283
    row_to_read = 22283
    file_input = pd.read_csv("GSE2034-22071 (edited).csv", nrows = row_to_read)
    file_output= pd.read_csv("mapping_sample_to_class_full.csv", usecols = ['GEO asscession number', 'relapse (1=True)'])

    # files to be used to get pathways and their gene expression
    # default rows_to_read_file_pathway = 1329
    rows_to_read_file_pathway = 1329
    file_ref_name = "accession_number_to_entrez_id.csv"
    file_to_convert_name = "GSE2034-22071 (edited).csv"
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)

    # get gene order id with its name
    list_gene_name = []
    for i in range(0, row_to_read):
        # add element with this format (gene_order_id, gene_name)
        gene_name = []
        gene_name.append(i)
        gene_name.append(file_input.loc[i, "ID_REF"])
        list_gene_name.append(gene_name)
    
    # get list of pathway name
    list_pathway_name = []
    for i in range(0, rows_to_read_file_pathway):
        pathway_name = []
        pathway_name.append(i)
        pathway_name.append(file_pathway.loc[i, "PATHWAY_NAME"])
        list_pathway_name.append(pathway_name)
    
    # consider non-relapse and relapse (not in specific period of time)
    sample_relapse = file_output.loc[file_output['relapse (1=True)'].isin(['1'])]
    sample_no_relapse = file_output.loc[file_output['relapse (1=True)'].isin(['0'])]

    # add GEO asscession number to each list
    list_sample_relapse = []
    for element in sample_relapse.loc[:, 'GEO asscession number']:
        list_sample_relapse.append(element)

    list_sample_no_relapse = []
    for element in sample_no_relapse.loc[:, 'GEO asscession number']:
        list_sample_no_relapse.append(element)

    # create list contains all samples
    list_sample_all = []
    list_sample_all.extend(list_sample_relapse)
    list_sample_all.extend(list_sample_no_relapse)
    print(" ## Size of data ##")
    print(" Number of samples with relapse : " + str(len(list_sample_relapse))) 
    print(" Number of samples with non-relapse : " + str(len(list_sample_no_relapse)))
    print(" Number of all samples : " + str(len(list_sample_all)))

    # list all gene expression to calculate mean and sd
    print("\n Gathering gene expressions are in progress ... ")
    list_all_gene_expression = []
    for line_index in range(0, row_to_read):
        for column in file_input.iloc[line_index, 1 : ]:
            list_all_gene_expression.append(column)

    mean_all_gene_expression = calculate.mean(list_all_gene_expression)
    sd_all_gene_expression = calculate.sd(list_all_gene_expression)
    max_all_gene_expression = max(list_all_gene_expression)
    min_all_gene_expression = min(list_all_gene_expression)

    print()
    print(" Mean of all gene expression : " + str(mean_all_gene_expression))
    print(" SD of all gene expression : " + str(sd_all_gene_expression))
    print(" Max of all gene expression : " + str(max_all_gene_expression))
    print(" Min of all gene expression : " + str(min_all_gene_expression))

# get gene expression of each pathway of each sample in testing set
    samples_with_activity = {}
    for element_index in range(0, 3):
        print()
        print("Creating pathways for samples_testing_all #" + str(element_index + 1) + "  is in progress ...")
        print(str(len(list_sample_all) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_all[element_index]
        pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression, min_of_data = min_all_gene_expression, \
                    method = "wide_scaling")

        sample.append(sample_name)
        sample.append(pathways)
        samples_with_activity[element_index] = sample 



if __name__ == '__main__':
    main()