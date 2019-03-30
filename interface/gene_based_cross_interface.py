import add_ons

def main():
    print()
    print("------------------------------------------------------------------------------------------------------------------------")
    print(" # Method : Gene-Based Classification")
    print(" # Experiment : Cross-Dataset")
    print(" # This method requires 2 datasets.")
    print(" # You will be asked to provide related files and required information about them including ")
    print(" #   [1] A file contains mapping between gene probe IDs and samples of the first dataset")
    print(" #   [2] Number of rows of the file containing mapping between gene probe IDs and samples of the first dataset to be read")
    print(" #   [3] A file contains mapping between samples and their class of the first dataset")  
    print(" #   [4] A file contains mapping between gene probe IDs and samples of the second dataset")
    print(" #   [5] Number of rows of the file contains mapping between gene probe IDs and samples of the second dataset to be read")
    print(" #   [6] A file contains mapping samples and their class of the second dataset")
    print(" # These files must follow a required format shown in file_format.pdf")
    print(" #")
    print(" # You will be asked to provide required information to conduct an experiment including")
    print(" #   [1] Number of epochs")
    print(" #   [2] Number of folds")
    print(" #   [3] Number of top-ranked feature")
    print(" #")
    print(" # You will be asked for the name of an output file.")
    print("------------------------------------------------------------------------------------------------------------------------")
    print()

    # prepare variables
    file_gene_first_dataset_name = None
    file_gene_second_dataset_name = None

    row_to_read_file_gene_first_dataset = None
    row_to_read_file_gene_second_dataset = None
    
    file_output_first_dataset_name = None
    file_output_second_dataset_name = None

    epoch = None
    num_of_folds = None
    number_of_ranked_gene = None

    file_name = None


    print(" # Enter required information about the first dataset ")
    print(" 1. Enter name of a file containing mapping between probes IDs and samples of the first dataset ")
    file_gene_first_dataset_name = add_ons.checkFileValid()
    print()

    print(" 2. Enter number of rows of this file to be read ")
    while True:
        row_to_read_file_gene_first_dataset = input(" Number of rows : ")
        if (int(row_to_read_file_gene_first_dataset) < 1):
            print("WARNING : Number of rows cannot be lower than 1.")
        else:
            break
    print()

    print(" 3. Enter name of file containing mapping between samples and their class of the first dataset")
    file_output_first_dataset_name = add_ons.checkFileValid()
    print()

    print(" # Enter required information about the second dataset ")
    print(" 1. Enter name of a file containing mapping between probes IDs and samples of the second dataset ")
    file_gene_second_dataset_name = add_ons.checkFileValid()
    print()

    print(" 2. Enter number of rows of this file to be read ")
    while True:
        row_to_read_file_gene_second_dataset = input(" Number of rows : ")
        if (int(row_to_read_file_gene_second_dataset) < 1):
            print("WARNING : Number of rows cannot be lower than 1.")
        else:
            break
    print()

    print(" 3. Enter name of a containing mapping between samples and their class of the second dataset")
    file_output_second_dataset_name = add_ons.checkFileValid()
    print()

    print(" # Enter required information to conduct an experiment")
    print(" 1. Enter number of epochs ")
    while True:
        epoch = input(" Epochs : ")

        if (epoch.isnumeric() == False):
            print(" WARNING : Number of epochs must be numeric.")
        elif (int(epoch) <= 0):
            print(" WARINING : Number of epochs must be greater than 0.")
        else:
            break
    print()

    print(" 2. Enter number of folds ")
    while True:
        num_of_folds = input("Number of folds: ")
        if (num_of_folds.isnumeric() == False):
            print("WARNING : Invalid input must be numeric")
        
        # these conditions are not available in mock-up
        # elif(int(num_of_folds) > len(list_sample_relapse_second_dataset)):
        #     print("WARNING : Number of folds exceeds the size of the 1st dataset")
        # elif(int(num_of_folds) > len(list_sample_relapse_second_dataset)):
        #     print("WARNING : Number of folds exceeds the size of the 2nd dataset")

        elif(int(num_of_folds) <= 1):
            print("WARNING : Number of folds cannot lower than or equal to 1")
        else:
            break
    num_of_folds = int(num_of_folds)    
    print()

    print(" 3. Enter number of top-ranked features")
    while True:
        number_of_ranked_gene = input("Number of top-ranked features: ")
        if(number_of_ranked_gene.isnumeric() == False):
            print(" WARNING : Number of top-ranked features must be numeric.")

        # these conditions are not available in mock-up
        # elif(int(number_of_ranked_gene) > row_to_read_file_gene_first_dataset):
        #     print(" WARINING : Number of top-ranked features must not exceed available genes from the first file.")

        elif(int(number_of_ranked_gene) <= 0):
            print(" WARNING : Number of top-ranked features must not be lower than or equal to 0.")    
        else:
            break
    print()

    file_name = input(" # Enter name of an output file : ")

    


if __name__ == "__main()__":
    main()