import os
import config as cf
import sys
import data_handler
import lib_gene_based

def main():

    # test_file = os.path.abspath("C:\Users\thatchayut\Desktop\LinuxEnvironment\classification")

    # while True :
    #     rep = input('Try for another file ?  y/n ?')
    #     if rep == 'n':
    #         sys.exit(0)
    
    # print("aaaa")

    lib_gene_based.gene_based()
   

    # flag = data_handler.validateFile(file_training_input_name)
    # print(flag)

if __name__ == "__main__":
    main()