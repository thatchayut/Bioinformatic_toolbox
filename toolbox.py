import gene_based
import gene_based_cross
import mean
import mean_cross
import pac
import pac_cross
import llr
import llr_cross

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

    if (method_number == "1"):
        if (experiment_number == "1"):
            gene_based.main()
        elif (experiment_number == "2"):
            gene_based_cross.main()

    elif (method_number == "2"):
        if (experiment_number == "1"):
            mean.main()
        elif (experiment_number == "2"):
            mean_cross.main()

    elif (method_number == "3"):
        if (experiment_number == "1"):
            pac.main()
        elif (experiment_number == "2"):
            pac_cross.main()
            
    elif (method_number == "4"):
        if (experiment_number == "1"):
            llr.main()
        elif (experiment_number == "2"):
            llr_cross.main()


if __name__ == "__main__":
    main()