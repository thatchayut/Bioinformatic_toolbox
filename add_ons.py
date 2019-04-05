import os

def checkFileValid():
    check_file_type = False
    file_name = None
    while check_file_type is False:
        input_value = input(" Enter input file : ")
        split_name = input_value.split(".")

        if (split_name[len(split_name) - 1] != 'csv'):
            print(" WARNING : File must be a comma-separated file (.csv) .")
            continue
        else:
            file_name = input_value
            check_file_type = True

    return file_name

def getFile():
    file_name = None
    check_file_exist = False

    while check_file_exist is False:
        try :
            input_file_name = checkFileValid() # mod_line
            if os.stat(input_file_name).st_size < 0: 
                print(" WARNING : " + str(input_file_name) + " is empty.") 
            else:
                file_name = input_file_name
                check_file_exist = True
        except OSError:
            print(" WARNING : " + str(input_file_name) + " is not found in this directory.") # mod_line
    
    return file_name