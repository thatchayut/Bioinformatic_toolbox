#!/usr/bin/python
import pandas as pd

def main():
    # prepare files to be used
    cols_to_read_file_to_convert = ["ID_REF", "GSM36777"]
    rows_to_read_file_pathway = 1
    file_ref = pd.read_csv("accession_number_to_entrez_id.csv")
    file_to_convert = pd.read_csv("GSE2034-22071 (edited).csv", usecols = cols_to_read_file_to_convert)
    file_pathway = pd.read_csv("c2.cp.v6.2.entrez.gmt.csv", nrows = rows_to_read_file_pathway)

    # list all probe id
    list_probe_id = []
    for element in file_to_convert.loc[:, 'ID_REF']:
        list_probe_id.append(element)

    num_of_probe_id = len(list_probe_id)
    count_not_found = 0
    # scan probe id by its entrez id
    list_entrez_id = []
    print("Scanning in progress ...")
    print()
    for index in range(0, num_of_probe_id):
        entrez_id = file_ref.loc[file_ref['From'].isin([list_probe_id[index]])]
        # add to list if entrez_id is found
        if (entrez_id.empty):
            list_entrez_id.append("none")
            count_not_found += 1
        else:
            list_entrez_id.append(entrez_id.iloc[0][1])  
    num_of_available_data = num_of_probe_id - count_not_found

    print(file_pathway)

    # create dictionary to collect each pathway
    pathways = {}
    for i in range(0, rows_to_read_file_pathway):
        key = i
        pathway = []
        pathway_name = file_pathway.iloc[i, 0]
        list_gene_expression = []
        for element in file_pathway.iloc[i, 2:-1]:
            if (str(element).isnumeric()):
                # check if genes in each pathway have their expression value
                list_gene_same_entrez = []
                for i in range(0, len(list_entrez_id)):
                    if (element == list_entrez_id[i]):
                        # WARNING : NEED to adjust column in file_to_convert.iloc[i, 1] for use further
                        list_gene_same_entrez.append(file_to_convert.iloc[i, 1])
                    else:
                        list_gene_same_entrez.append(0.0)

                # calculate genes expression of the same entrez id by using their average
                num_of_same_entrez = len(list_gene_same_entrez)
                avg_gene_expression = (sum(list_gene_same_entrez) / num_of_same_entrez)
                avg_gene_expression = round(avg_gene_expression, 7)
                # print(list_gene_same_entrez)
                list_gene_expression.append(avg_gene_expression)
        # combine pathway name with its gene expressions
        pathway.append(pathway_name)
        pathway.append(list_gene_expression)
        
        pathways[key] = pathway
                    



if __name__ == '__main__':
    main()