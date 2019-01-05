#!/usr/bin/python
import pandas as pd

def main():
    # prepare files to be used
    cols_to_read_file_to_convert = ["ID_REF", "GSM36777"]
    rows_to_read_file_pathway = 1329
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
                list_gene_name_with_expression = []
                # check if genes in each pathway have their expression value
                list_gene_same_entrez = []
                for list_entrez_id_index in range(0, len(list_entrez_id)):
                    if (element == list_entrez_id[list_entrez_id_index]):
                        # WARNING : NEED to adjust column in file_to_convert.iloc[i, 1] for use further
                        list_gene_same_entrez.append(file_to_convert.iloc[list_entrez_id_index, 1])
                # if gene expression is not found, assume it as 'zero'
                if not list_gene_same_entrez:
                    list_gene_same_entrez.append(0.0)
                
                # calculate genes expression of the same entrez id by using their average
                num_of_same_entrez = len(list_gene_same_entrez)
                avg_gene_expression = (sum(list_gene_same_entrez) / num_of_same_entrez)
                avg_gene_expression = round(avg_gene_expression, 7)
                
                # create a gene list in format [gene entrez id, gene expression]
                list_gene_name_with_expression.append(element)
                list_gene_name_with_expression.append(avg_gene_expression)
                list_gene_expression.append(list_gene_name_with_expression)
                # list_gene_expression.append(element)
                # list_gene_expression.append(avg_gene_expression)

        # combine pathway name with its gene expressions
        pathway.append(pathway_name)
        pathway.append(list_gene_expression)
        
        pathways[key] = pathway
    print(pathways[1])
    # print(len(pathways))
                    



if __name__ == '__main__':
    main()