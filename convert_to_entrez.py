import xlwt
import pandas as pd

def main():
    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Gene Probe Ids to Gene Entrez IDs Converter")
    print(" # Purpose : This program is used to convert gene probe id to gene entrez id")
    print(" # You will be asked to provide related files and required information about them including ")
    print(" #   [1] A file contains pathways and their member genes")
    print(" #   [2] A file contains mapping between gene probe IDs and gene entrez IDs")  
    print(" #   [3] A file contains mapping between samples and their gene expression")  
    print(" # These files must follow a required format shown in file_format.pdf")
    print(" #")
    print(" # You will be asked for the name of an output file.")
    print(" # An output file will be created in directory 'result'.")
    print("------------------------------------------------------------------------------------------------------------------------")
    print()
    # ask for file name
    file_name = input("Enter file name : ")

    #initiate worksheet
    workbook = xlwt.Workbook()
    worksheet = workbook.add_sheet('worksheet')
    style = xlwt.easyxf("align : horiz right")

    # get file mapping between gene probe id  and gene entrez id
    file_ref_name = (" # Enter a file mapping between gene probe id and gene entrez id : ")
    file_ref_name = input(" Enter file : ")

    # get file mapping between samples and their gene expression
    file_to_convert_name = (" # Enter a file mapping between samples and their gene expression : ")
    file_to_convert_name = input(" Enter file : ")

    # prepare files to be used
    cols_to_read = ["ID_REF"]
    file_ref = pd.read_csv(file_ref_name)
    file_to_convert = pd.read_csv(file_to_convert_name, usecols = cols_to_read)


    # list all probe id
    list_probe_id = []
    for element in file_to_convert.loc[:, 'ID_REF']:
        list_probe_id.append(element)
    
    num_of_probe_id = len(list_probe_id)
    count_not_found = 0
    # scan and replace probe id by its entrez id
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
    
    print(list_entrez_id)
    print()
    print("Avalable data : " + str(num_of_available_data))
    print("Number of non-existed entrez id : " + str(count_not_found))

    # write to file
    print("Writing in progress ...")
    print()
    worksheet.write(0, 0, "ENTREZ_ID")
    worksheet.write(0, 1, "ID_REF")
    for i in range(0, num_of_probe_id):
        if list_entrez_id[i] != "none":
            worksheet.write(i + 1, 0, int(list_entrez_id[i]), style)
            worksheet.write(i + 1, 1, (list_probe_id[i]), style)
        else:
            worksheet.write(i + 1, 0, str(list_entrez_id[i]), style)
            worksheet.write(i + 1, 1, (list_probe_id[i]), style)
    print("Done ...")
    workbook.save(file_name + ".csv")

if __name__ == '__main__':
    main()