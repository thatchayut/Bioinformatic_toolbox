import xlsxwriter
import pandas as pd

def main():

    # ask for output file name
    file_name = input(" Enter output file's name : ")

    # prepare files
    # prepare file containing probe ids and gene expressions
    file_gene_expression_name = "GSE3494_GPL96.csv"
    # default row_to_read_file_gene_expression = 22283
    row_to_read_file_gene_expression = 22283
    file_gene_expression  = pd.read_csv(file_gene_expression_name, nrows = row_to_read_file_gene_expression)

    # prepare file containing mapping between ID and Sample ID
    file_mapping_sample_name = "mapping_id_to_gse_3494.csv"
    # default row_to_read_file_mapping_sample = 251
    row_to_read_file_mapping_sample = 251
    file_mapping_sample = pd.read_csv(file_mapping_sample_name, nrows = row_to_read_file_mapping_sample)
    
    # prepare file of sample status
    file_sample_status_name = "GSE3494_sample_status.csv"
    row_to_read_file_sample_status = 251
    file_sample_status = pd.read_csv(file_sample_status_name, usecols = ["INDEX (ID)", "p53 seq mut status (p53+=mutant; p53-=wt)"])

    # Create list of mapping betwenn ID and Sample ID
    list_mapping_sample = []
    for row_index in range(0, row_to_read_file_mapping_sample):
        list_each_sample = []
        id_x = file_mapping_sample.loc[row_index, "INDEX (ID)"]
        sample_id = file_mapping_sample.loc[row_index, "SAMPLE_ID"]

        list_each_sample.append(id_x)
        list_each_sample.append(sample_id)

        list_mapping_sample.append(list_each_sample)
    
    # create list of sample status
    list_sample_status = []
    for row_index in range(0, row_to_read_file_sample_status):
        list_each_sample = []
        id_x = file_sample_status.loc[row_index, "INDEX (ID)"]
        p53_status = file_sample_status.loc[row_index, "p53 seq mut status (p53+=mutant; p53-=wt)"]

        list_each_sample.append(id_x)
        list_each_sample.append(p53_status)

        list_sample_status.append(list_each_sample)

    # convert id_x in list_sample_status to sample_id
    list_sample_status_converted = []
    for mapping_index in range(0, len(list_mapping_sample)):
        for sample_index in range(0, len(list_sample_status)):
            idx_mapping = list_mapping_sample[mapping_index][0]
            idx_sample_status = list_sample_status[sample_index][0]

            if (idx_mapping == idx_sample_status):
                sample_status = []
                sample_id = list_mapping_sample[mapping_index][1]
                p53_status = list_sample_status[sample_index][1]

                sample_status.append(sample_id)
                sample_status.append(p53_status)

                list_sample_status_converted.append(sample_status)
    
    # write to file
    workbook = xlsxwriter.Workbook(file_name + ".xlsx")
    worksheet = workbook.add_worksheet()

    print(" Writing in progress ...")
    print()

    worksheet.write(0, 0, "SAMPLE_ID")
    worksheet.write(0, 1, "p53 seq mut status (p53+=mutant; p53-=wt)")
    worksheet.write(0, 2, "relapse (1=True)")

    for row_index in range(0, len(list_sample_status_converted)):
        sample_id = list_sample_status_converted[row_index][0]
        p53_status = list_sample_status_converted[row_index][1]

        relapse_status = None
        if (p53_status == "p53+"):
            relapse_status = 1
        elif (p53_status == "p53-"):
            relapse_status = 0

        worksheet.write(row_index + 1, 0, sample_id)
        worksheet.write(row_index + 1, 1, p53_status)
        worksheet.write(row_index + 1, 2, relapse_status)

    workbook.close()
    print(" Creating file is done ...")

if __name__ == "__main__":
    main()