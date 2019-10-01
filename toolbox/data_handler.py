import os

# def checkFileValid():
#     check_file_type = False
#     file_name = None
#     while check_file_type is False:
#         input_value = input(" Enter input file : ")
#         split_name = input_value.split(".")

#         if (split_name[len(split_name) - 1] != 'csv'):
#             print(" WARNING : File must be a comma-separated file (.csv) .")
#             continue
#         else:
#             file_name = input_value
#             check_file_type = True

#     return file_name

# def getFile():
#     file_name = None
#     check_file_exist = False

#     while check_file_exist is False:
#         try :
#             input_file_name = checkFileValid() 
#             if os.stat(input_file_name).st_size < 0: 
#                 print(" WARNING : " + str(input_file_name) + " is empty.") 
#             else:
#                 file_name = input_file_name
#                 check_file_exist = True
#         except OSError:
#             print(" WARNING : " + str(input_file_name) + " is not found in this directory.") 
    
#     return file_name

def checkFileValid(file_name):
    split_name = file_name.split(".")

    if (split_name[len(split_name) - 1] != 'csv'):
        print(" WARNING : File must be a comma-separated file (.csv) . ")
        return False
    else:
        return True


def validateFile(file_name):
    try:
        validation_flag = checkFileValid(file_name)
        if (validation_flag is True):
            if os.stat(file_name).st_size < 0: 
                print(" WARNING : " + str(file_name) + " is empty.") 
                return False
            else:
                return True 
        else : 
            print(" Provide a correct file in config.py ...")    
            return False 
    except OSError:
        print(" WARNING : " + str(file_name) + " is not found in this directory.") 
        print(" Provide a correct file in config.py ...")    
        return False 

def validateRowToRead(row_to_read):
    if(row_to_read.isnumeric() == False):
        print(" WARNING : Number of rows must be numeric.")
        print(" Provide correct data in config.py ...")   
        return False
    elif (int(row_to_read) < 1):
        print(" WARNING : Number of rows cannot be lower than 1.")
        print(" Provide correct data in config.py ...")   
        return False
    else:
        return True

def validateEpochs(epochs):
    if (epochs.isnumeric() == False):
        print(" WARNING : Number of epochs must be numeric.")
        print(" Provide correct data in config.py ...") 
        return False 
    elif (int(epochs) <= 0):
        print(" WARINING : Number of epochs must be greater than 0.")
        print(" Provide correct data in config.py ...") 
        return False 
    else:
        return True

def validateNumofFolds(num_of_folds, list_sample_relapse, list_sample_no_relapse):
    if (num_of_folds.isnumeric() == False):
        print(" WARNING : Number of folds must be numeric")
        print(" Provide correct data in config.py ...") 
        return False
    elif(int(num_of_folds) > len(list_sample_relapse)):
        print("WARNING : Number of folds exceeds the size of samples in class relapse")
        print(" Provide correct data in config.py ...") 
        return False
    elif(int(num_of_folds) > len(list_sample_no_relapse)):
        print("WARNING : Number of folds exceeds the size of samples in class non-relapse")
        print(" Provide correct data in config.py ...") 
        return False

    elif(int(num_of_folds) <= 1):
        print(" WARNING : Number of folds cannot lower than or equal to 1")
        print(" Provide correct data in config.py ...") 
        return False
    else:
        return True

def validateNumofRankedGenes(num_of_ranked_genes, row_to_read):
    if (num_of_ranked_genes.isnumeric() == False):
        print(" WARNING : Number of top-ranked features must be numeric.")
        print(" Provide correct data in config.py ...")
        return False
    elif(int(num_of_ranked_genes) > row_to_read):
        print(" WARINING : Number of top-ranked features must not exceed available genes from the first file.")
        print(" Provide correct data in config.py ...")
        return False
    elif (int(num_of_ranked_genes) <= 0):
        print(" WARNING : Number of top-ranked features must not be lower than or equal to 0.")    
        print(" Provide correct data in config.py ...")
        return False
    else:
        return True
