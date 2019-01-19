import xlwt
import xlsxwriter
import pandas as pd
import calculate

def main():

    # ask for scaling method
    print("\n Scaling methods: [1] z-score // [2] narrow scaling (range [0,1]) // [3] wide scaling (range [-1,1])") 
    method_id = None
    while True:
        method_id = input(" Enter method id : ")
        if (method_id not in ["1", "2", "3"]):
            print(" WARNING : Invalid method id")
        else:
            break
    # ask for output file name
    file_name = input(" Enter file name : ")

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
    # print(list_pathway_name)

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
    samples_with_gene_expression = {}
    for element_index in range(0, len(list_sample_all)):
        print()
        print("Creating pathways for samples_with_gene_expression #" + str(element_index + 1) + "  is in progress ...")
        print(str(len(list_sample_all) - (element_index + 1)) + " samples left")
        print()

        sample = []
        sample_name = list_sample_all[element_index]

        if (method_id is "1"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, mean_of_data = mean_all_gene_expression, sd_of_data = sd_all_gene_expression, \
                        method = "z_score")
        elif (method_id is "2"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression, min_of_data = min_all_gene_expression, \
                        method = "narrow_scaling")            
        elif (method_id is "3"):
            pathways = calculate.getPathway(file_ref_name, file_to_convert_name, file_pathway_name, sample_name, rows_to_read_file_pathway, max_of_data = max_all_gene_expression, min_of_data = min_all_gene_expression, \
                        method = "wide_scaling")       

        sample.append(sample_name)
        sample.append(pathways)
        samples_with_gene_expression[element_index] = sample 

        # print(samples_with_gene_expression[0])

    # calculate pathway activity using 'Mean'
    samples_with_pathway_activity = {}
    # { GSM1234, {0: ['KEGG_GLYCOLYSIS_GLUCONEOGENESIS', [[55902, 0.0], [2645, 0.0], ...}}
    for samples_index in range(0, len(samples_with_gene_expression)):
        sample = []
        list_pathway = []
        for pathway_index in range(0, len(samples_with_gene_expression[samples_index][1])):
            list_gene_expression_in_pathway = []
            pathway = []
            for gene_index in range(0, len(samples_with_gene_expression[samples_index][1][pathway_index][1])):
                # print(sample_relapse[samples_index][1][pathway_index][gene_index][1])
                # print(gene_index)
                gene_expression = samples_with_gene_expression[samples_index][1][pathway_index][1][gene_index][1]
                list_gene_expression_in_pathway.append(gene_expression)

            # data to collect as pathway activity
            pathway_name = samples_with_gene_expression[samples_index][1][pathway_index][0]
            pathway_activity = calculate.mean(list_gene_expression_in_pathway)

            pathway.append(pathway_name)
            pathway.append(pathway_activity)

            list_pathway.append(pathway)
            # # add to list only if it is in feature set
            # if (pathway_name in list_top_ranked_pathways):
            #     list_pathway.append(pathway)
        
        sample_name = samples_with_gene_expression[samples_index][0]
        sample.append(sample_name)
        sample.append(list_pathway)
        samples_with_pathway_activity[samples_index] = sample    
    
    # print(samples_with_pathway_activity[0])

    # write to file
    #initiate worksheet
    # workbook = xlwt.Workbook()
    # worksheet = workbook.add_sheet(file_name)
    # style = xlwt.easyxf("align : horiz right")
    workbook = xlsxwriter.Workbook(file_name + ".xlsx")
    worksheet = workbook.add_worksheet()

    print("Writing in progress ...")
    print()
    worksheet.write(0, 0, "PATHWAY_NAME")

    # write pathway name
    for i in range(0, len(list_pathway_name)):
        worksheet.write(i + 1, 0, str(list_pathway_name[i][1]))

    # write data of each sample
    num_of_total_line_in_file = len(list_pathway_name) + 1
    for sample_index in range(0, len(samples_with_pathway_activity)):
        for line_index in range(0, num_of_total_line_in_file):
            # write column name (sample id)
            if (line_index == 0):
                sample_name = str(samples_with_pathway_activity[sample_index][0])
                worksheet.write(line_index, sample_index + 1, sample_name)
            else:
                pathway_activity = str(samples_with_pathway_activity[sample_index][1][line_index - 1][1])
                worksheet.write(line_index, sample_index + 1, pathway_activity)
    # print(samples_with_pathway_activity[0])
    # workbook.save(file_name + ".csv")
    workbook.close()
    print(" Creating file is done ...")

if __name__ == '__main__':
    main()