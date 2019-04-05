import random
import time
import numpy as np
import random
import os
import add_ons


def getFile():
    file_name = None
    check_file_exist = False

    while check_file_exist is False:
        try :
            input_file_name = add_ons.checkFileValid() # mod_line
            if os.stat(input_file_name).st_size < 0: 
                print(" WARNING : " + str(input_file_name) + " is empty.") 
            else:
                print(" VALID FILE")
                file_name = input_file_name
                check_file_exist = True
        except OSError:
            print(" WARNING : " + str(input_file_name) + " is not found.") # mod_line
    
    return file_name


def main():
    print("########################################################################################################################")
    print("########################################                                 ###############################################")
    print("######################################   TOOLBOX FOR BIOINFORMATIC TOOLS   #############################################")
    print("########################################                                 ###############################################")
    print("########################################################################################################################")
    print("#                                                                                                                      #")
    print("# This toolbox includes 4 classification methods including                                                             #")
    print("#    [1] Gene-Based Classification                                                                                     #")
    print("#    [2] Pathway-Based Classification Based on                                                                         #")
    print("#        [2.1] Mean                                                                                                    #")
    print("#        [2.2] Condition-Responsive Genes (CORGs)                                                                      #")
    print("#        [2.3] Log-Likelihood Ratio                                                                                    #")
    print("#                                                                                                                      #")
    print("# Each of these methods is divided into 2 types of the experiment which are                                            #")
    print("#    [1] Within Dataset Experiment                                                                                     #")
    print("#    [2] Cross-Dataset Experiment                                                                                      #")
    print("#                                                                                                                      #")
    print("# All files provided to this toolbox must strictly follow the format shown in file_format.pdf.                         #")
    print("########################################################################################################################")
    print()


    # filename = input(" Test file name : ") 

    # result_file = open("./result/" +str(filename) + ".txt", "w+")
    # result_file.write("testestestestes")
    print(" #### Enter name of the first file  ") # add_line
    file_name = getFile()
    print("file_name : " + str(file_name))


if __name__ == "__main__":
    main()