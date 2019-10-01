import os
import config as cf
import sys
import files_handler

def main():

    # test_file = os.path.abspath("C:\Users\thatchayut\Desktop\LinuxEnvironment\classification")

    # while True :
    #     rep = input('Try for another file ?  y/n ?')
    #     if rep == 'n':
    #         sys.exit(0)
    
    # print("aaaa")

    file_training_input_name = cf.gene_based["file_training_input_name"]

    flag = files_handler.validateFile(file_training_input_name)
    print(flag)

if __name__ == "__main__":
    main()