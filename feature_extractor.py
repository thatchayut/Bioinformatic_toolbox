import xlsxwriter
import pandas as pd

def main():

    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Feature Extractor")
    print(" # Purpose : This program is used to extract features into member genes.")
    print(" # You will be asked to provide related files and required information about them including ")
    print(" #   [1] A file contains pathways and their member genes")
    print(" #   [2] number of rows of the file contaning pathways and their member genes")
    print(" #   [3] A file contains mapping between gene probe IDs and gene entrez IDs")  
    print(" #   [4] Feature set to be extracted")
    print(" # These files must follow a required format shown in file_format.pdf")
    print(" #")
    print(" # You will be asked for the name of an output file.")
    print(" # An output file will be created in directory 'result'.")
    print("------------------------------------------------------------------------------------------------------------------------")
    print()

    # get a pathway reference file
    print(" # Enter a pathway reference file : ")
    # file_pathway_name = input(" Enter file : ")
    file_pathway_name = "c2.cp.v6.2.entrez.gmt.csv"
    print()

    print(" # Enter number of rows of this file to be read : ")
    # rows_to_read_file_pathway = input(" Number of rows : ")
    # rows_to_read_file_pathway = int(rows_to_read_file_pathway)
    rows_to_read_file_pathway = 1329
    print()

    # get file mapping between gene entrez id and gene probe id
    file_ref_name = (" # Enter a file mapping between gene entrez id and gene probe id : ")
    # file_ref_name = input(" Enter file : ")
    file_ref_name = "accession_number_to_entrez_id.csv"

    # get feature set to extract
    print(" # Enter a feature set to be extracted : ")
    feature_set = input(" Feature set : ")
    print()

    # get output file's name
    file_name = input(" # Name of output file : ")
    print()

    # prepare text file for results to be written in
    result_file = open("./result/" + str(file_name) + ".txt", "w+")

    # select option :
    print(" # Extract to [1] Gene entrez id [2] Gene Probe id : ")
    extrct_choice = input(" Enter your choice : ")
    print()

    # get proper format of feature set
    feature_set = feature_set.replace(" ", '')
    feature_set = feature_set.replace("[", '')
    feature_set = feature_set.replace("]", '')
    feature_set = feature_set.replace("'", '')

    # split feature set into features
    feature_set = feature_set.split(",")

    # get data from file
    file_pathway = pd.read_csv(file_pathway_name, nrows = rows_to_read_file_pathway)
    file_ref = pd.read_csv(file_ref_name)

    # mapping feature to pathway
    list_gene_entrez_id = []
    list_gene_probe_id = []

    for feature_index in range(0, len(feature_set)):
        feature = feature_set[feature_index]

        for pathway_index in range(0, rows_to_read_file_pathway):
            pathway_name = file_pathway.iloc[pathway_index, 0]
            
            if (pathway_name == feature):
                # get gene entrez id of genes in this pathway
                for element in file_pathway.iloc[pathway_index, 2:-1]:
                    if (str(element).isnumeric()):
                        if (str(element) not in list_gene_entrez_id):
                            list_gene_entrez_id.append(str(element))

    if (extrct_choice == "2"):        
        # mapping gene entrez id to gene probe id 
        for gene_entrez_index in range(0, len(list_gene_entrez_id)):
            gene_entrez_id = list_gene_entrez_id[gene_entrez_index]

            gene_probe_id_dataframe = file_ref.loc[file_ref['To'].isin([gene_entrez_id])]
            
            for index, row in gene_probe_id_dataframe.iterrows():
                gene_probe_id = row['From']
                if (gene_probe_id not in list_gene_probe_id):
                    list_gene_probe_id.append(gene_probe_id)

    # write to file
    print(" Process : Write to an output file ...")
    result_file.write("Gene probe ids from feature set : \n")
    result_file.write(str(feature_set) + "\n")
    result_file.write("\n")

    if (extrct_choice == "1"):
        result_file.write("Gene entrez id : \n")
        for element in list_gene_entrez_id:
            result_file.write(str(element) + "\n")

    elif (extrct_choice == "2"):
        result_file.write("Gene probe id : \n")
        for element in list_gene_probe_id:
            result_file.write(str(element) + "\n")
    print(" Process : Done ...")

if __name__ == "__main__":
    main()