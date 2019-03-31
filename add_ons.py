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