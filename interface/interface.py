import gene_based_interface
import gene_based_cross_interface
import mean_interface
import mean_cross_interface
import pac_interface
import pac_cross_interface
import llr_interface
import llr_cross_interface

def main():
    print("########################################################################################################################")
    print("#                                                                                                                      #")
    print("#                                                                                                                      #")
    print("#                                                                                                                      #")
    print("#                                                                                                                      #")
    print("#                                                                                                                      #")
    print("########################################################################################################################")
    print()

    print(" 1. Select Method [1] Gene-Based  [2] Mean-Based  [3] CORGs-Based [4] LLR-Based ")
    method_number = None
    list_available_method_number = ["1", "2", "3", "4"]
    while True:
        method_number = input(" Enter method number : ")
        if (method_number in list_available_method_number):
            break
        else:
            print(" WARNING : Available choices are " + str(list_available_method_number[0] + " to " + str(list_available_method_number[len(list_available_method_number) - 1]) + "."))
    
    print()
    print(" 2. Select Type of Experiment [1] Within Dataset Experiment [2] Cross-Dataset Experiment")
    experiment_number = None
    list_available_experiment_number = ["1", "2"]
    while True:
        experiment_number = input(" Enter experiment number :")
        if (experiment_number in list_available_experiment_number):
            break
        else:
            print(" WARNING : Available choices are " + str(list_available_experiment_number[0] + " to " + str(list_available_experiment_number[len(list_available_experiment_number) - 1]) + "."))

    # gene_based_cross_interface.main()
    gene_based_interface.main()


if __name__ == "__main__":
    main()