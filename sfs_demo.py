import calculate
import random
import time 
from copy import deepcopy
from sklearn.metrics import roc_auc_score

def main():
    num_of_pathways = 10
    num_of_samples_each_class = 10

    samples_relapse = {}
    samples_no_relapse = {}
    samples_test = {}
    sample_test_evaluation = {}

    list_desired_output = []

    list_desired_output_evaluation = []

    # create demo samples
    for sample_index in range(0, num_of_samples_each_class):
        list_relapse_pathway_activity = []    
        # create list of pathway activity
        random.seed(time.time())
        for i in range(0, num_of_pathways):
            pathway_activity = random.randint(0, 75)
            list_relapse_pathway_activity.append(float(pathway_activity))
        samples_relapse[sample_index] = list_relapse_pathway_activity

    for sample_index in range(0, num_of_samples_each_class):
        list_non_relapse_pathway_activity = []
        random.seed(time.time())
        for i in range(0, num_of_pathways):
            pathway_activity = random.randint(25, 100)
            list_non_relapse_pathway_activity.append(float(pathway_activity))
        samples_no_relapse[sample_index] = list_non_relapse_pathway_activity

    for sample_index in range(0, num_of_samples_each_class):
        list_test_pathway_activity = []
        random.seed(time.time())

        class_options = random.randint(0,1)

        if (class_options == 1):
            list_desired_output.append(class_options)
            for i in range(0, num_of_pathways):
                pathway_activity = random.randint(0, 75)
                list_test_pathway_activity.append(float(pathway_activity))   
            samples_test[sample_index] = list_test_pathway_activity     
        elif (class_options == 0):
            list_desired_output.append(class_options)
            for i in range(0, num_of_pathways):
                pathway_activity = random.randint(25, 100)
                list_test_pathway_activity.append(float(pathway_activity))   
            samples_test[sample_index] = list_test_pathway_activity   

    for sample_index in range(0, num_of_samples_each_class):
        list_test_pathway_activity = []
        random.seed(time.time())

        class_options = random.randint(0,1)

        if (class_options == 1):
            list_desired_output_evaluation.append(class_options)
            for i in range(0, num_of_pathways):
                pathway_activity = random.randint(0, 75)
                list_test_pathway_activity.append(float(pathway_activity))   
            sample_test_evaluation[sample_index] = list_test_pathway_activity     
        elif (class_options == 0):
            list_desired_output_evaluation.append(class_options)
            for i in range(0, num_of_pathways):
                pathway_activity = random.randint(25, 100)
                list_test_pathway_activity.append(float(pathway_activity))   
            sample_test_evaluation[sample_index] = list_test_pathway_activity   

    # create demo pathway name
    list_pathway_name = []
    # create list of pathway name 
    for i in range(0, num_of_pathways):
        pathway_name = "pathway_" + str(i)
        list_pathway_name.append(pathway_name)

    print(" #### Data preparation ####")
    print(" Pathway name : " + str(list_pathway_name))
    print(" List desired output : " + str(list_desired_output))

    # Sequential Forward Selection
    check_finish = False
    count_iteration = 0
    max_auc_score_over_all_features = 0
    feature_set_final = []
    check_improve_auc = True

    list_pathway_selected = []
    while (check_finish == False):

        if (check_improve_auc == True):
            max_auc_in_consider = 0
            list_pathway = []
            for pathway_index in range(0, num_of_pathways):
                list_pathway_to_consider = deepcopy(list_pathway_selected)

                # if (len(list_pathway_to_consider) == 1) and (pathway_index not in list_pathway_selected):
                if (pathway_index not in list_pathway_selected):
                    list_pathway_to_consider.extend([pathway_index])

                    # collect pathway activity of each class to be used in lda 
                    input_relapse_to_test = []
                    for sample_index in range(0, len(samples_relapse)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index in range(0, len(samples_relapse[sample_index])):
                            if (pathway_index in list_pathway_to_consider):
                                pathway_activity_to_test = samples_relapse[sample_index][pathway_index]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_relapse_to_test.append(list_pathway_each_sample_to_test)
                    # print(" input_relapse_to_test : " + str(input_relapse_to_test))
                    # print()

                    input_no_relapse_to_test = []
                    for sample_index in range(0, len(samples_no_relapse)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index in range(0, len(samples_no_relapse[sample_index])):
                            if (pathway_index in list_pathway_to_consider):
                                pathway_activity_to_test = samples_no_relapse[sample_index][pathway_index]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_no_relapse_to_test.append(list_pathway_each_sample_to_test)
                    # print(" input_no_relapse_to_test : " + str(input_no_relapse_to_test))
                    # print()
                    
                    input_test_to_test = []
                    for sample_index in range(0, len(samples_test)):
                        list_pathway_each_sample_to_test = []
                        for pathway_index in range(0, len(samples_test[sample_index])):
                            if (pathway_index in list_pathway_to_consider):
                                pathway_activity_to_test = samples_test[sample_index][pathway_index]
                                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
                        input_test_to_test.append(list_pathway_each_sample_to_test)
                    # print(" input_test_to_test : " + str(input_test_to_test))

                    list_actual_output = calculate.lda(input_test_to_test, input_relapse_to_test, input_no_relapse_to_test)
                    auc_score = roc_auc_score(list_desired_output, list_actual_output)

                    print(" list_actual_output : " + str(list_actual_output))
                    print(" list_desired_output : " + str(list_desired_output))
                    print(" AUC score : " + str(auc_score))

                    if (auc_score >= max_auc_in_consider):
                        max_auc_in_consider = auc_score
                        list_pathway = deepcopy(list_pathway_to_consider)
                    print("list_pathway_to_consider : " + str(list_pathway_to_consider))
                    print("max_auc_in_consider : " + str(max_auc_in_consider))

            if (max_auc_in_consider >= max_auc_score_over_all_features):
                max_auc_score_over_all_features = max_auc_in_consider
                list_pathway_selected = deepcopy(list_pathway)
            else:
                check_improve_auc = False
        else:
            check_finish = True

    # convert pathway index to pathway name
    list_feature_set = []
    for i in range(0, len(list_pathway_selected)):
        list_feature_set.append(list_pathway_name[list_pathway_selected[i]])


    print()
    print(" #### Summary ####")
    print("list_pathway_selected : " + str(list_pathway_selected))
    print("list_feature_set : " + str(list_feature_set))
    print("max_auc_score_over_all_features : " + str(max_auc_score_over_all_features))
    print()

    print(" #### Testing ####")
    # collect pathway activity of each class to be used in lda 
    input_relapse_evaluation = []
    for sample_index in range(0, len(samples_relapse)):
        list_pathway_each_sample_to_test = []
        for pathway_index in range(0, len(samples_relapse[sample_index])):
            if (pathway_index in list_pathway_selected):
                pathway_activity_to_test = samples_relapse[sample_index][pathway_index]
                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
        input_relapse_evaluation.append(list_pathway_each_sample_to_test)
    # print(" input_relapse_evaluation : " + str(input_relapse_evaluation))
    print()

    input_no_relapse_evaluation = []
    for sample_index in range(0, len(samples_no_relapse)):
        list_pathway_each_sample_to_test = []
        for pathway_index in range(0, len(samples_no_relapse[sample_index])):
            if (pathway_index in list_pathway_selected):
                pathway_activity_to_test = samples_no_relapse[sample_index][pathway_index]
                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
        input_no_relapse_evaluation.append(list_pathway_each_sample_to_test)
    
    input_test_evaluation = []
    for sample_index in range(0, len(sample_test_evaluation)):
        list_pathway_each_sample_to_test = []
        for pathway_index in range(0, len(sample_test_evaluation[sample_index])):
            if (pathway_index in list_pathway_selected):
                pathway_activity_to_test = sample_test_evaluation[sample_index][pathway_index]
                list_pathway_each_sample_to_test.append(pathway_activity_to_test)
        input_test_evaluation.append(list_pathway_each_sample_to_test)
    
    list_actual_output_evaluation = calculate.lda(input_test_evaluation, input_relapse_evaluation, input_no_relapse_evaluation)

    auc_score_evaluation = roc_auc_score(list_desired_output_evaluation, list_actual_output_evaluation)

    print(" #### Evaluation ####")
    print(" list_actual_output_evaluation : " + str(list_actual_output_evaluation))
    print(" list_desired_output_evaluation : " + str(list_desired_output_evaluation))
    print(" auc_score_evaluation : " + str(auc_score_evaluation))
    print(" list_feature_set : " +str(list_feature_set))


    # version : 1
    # while (check_finish == False):
    #     if (count_iteration >= num_of_pathways):
    #         check_finish = True
    #     else:
    #         max_auc_score = 0
    #         pathway_index_in_list = None

    #         pathway_to_test = []
    #         pathway_selected = []
    #         for index in range(count_iteration, num_of_pathways):
    #             pathway_to_test.extend([index])

    #             # skip the case which has only 1 pathway in the list
    #             # if (len(pathway_to_test) <= 1):
    #             #     continue
    #             # else:
    #             print(" **************************************************")
    #             print(" pathway_to_test : " + str(pathway_to_test))
    #             print()
    #             # collect pathway activity of each class to be used in lda 
    #             input_relapse_to_test = []
    #             for sample_index in range(0, len(samples_relapse)):
    #                 list_pathway_each_sample_to_test = []
    #                 for pathway_index in range(0, len(samples_relapse[sample_index])):
    #                     if (pathway_index in pathway_to_test):
    #                         pathway_activity_to_test = samples_relapse[sample_index][pathway_index]
    #                         list_pathway_each_sample_to_test.append(pathway_activity_to_test)
    #                 input_relapse_to_test.append(list_pathway_each_sample_to_test)
    #             # print(" input_relapse_to_test : " + str(input_relapse_to_test))
    #             # print()

    #             input_no_relapse_to_test = []
    #             for sample_index in range(0, len(samples_no_relapse)):
    #                 list_pathway_each_sample_to_test = []
    #                 for pathway_index in range(0, len(samples_no_relapse[sample_index])):
    #                     if (pathway_index in pathway_to_test):
    #                         pathway_activity_to_test = samples_no_relapse[sample_index][pathway_index]
    #                         list_pathway_each_sample_to_test.append(pathway_activity_to_test)
    #                 input_no_relapse_to_test.append(list_pathway_each_sample_to_test)
    #             # print(" input_no_relapse_to_test : " + str(input_no_relapse_to_test))
    #             # print()
                
    #             input_test_to_test = []
    #             for sample_index in range(0, len(samples_test)):
    #                 list_pathway_each_sample_to_test = []
    #                 for pathway_index in range(0, len(samples_test[sample_index])):
    #                     if (pathway_index in pathway_to_test):
    #                         pathway_activity_to_test = samples_test[sample_index][pathway_index]
    #                         list_pathway_each_sample_to_test.append(pathway_activity_to_test)
    #                 input_test_to_test.append(list_pathway_each_sample_to_test)
    #             # print(" input_test_to_test : " + str(input_test_to_test))

    #             list_actual_output = calculate.lda(input_test_to_test, input_relapse_to_test, input_no_relapse_to_test)
    #             auc_score = roc_auc_score(list_desired_output, list_actual_output)

    #             print(" list_actual_output : " + str(list_actual_output))
    #             print(" list_desired_output : " + str(list_desired_output))
    #             print(" AUC score : " + str(auc_score))

    #             # add to feature set if it gives maximum auc
    #             if (auc_score > max_auc_score):
    #                 max_auc_score = auc_score
    #                 pathway_selected.extend([index])

    #         print()
    #         print(" pathway_selected : " + str(pathway_selected))
    #         print(" max_auc_score : " + str(max_auc_score))

    #         # track maximum auc
    #         if (max_auc_score > max_auc_score_over_all_features):
    #             max_auc_score_over_all_features = max_auc_score
            
    #         print()
    #         print("max_auc_score_over_all_features : " + str(max_auc_score_over_all_features))
                

    #         count_iteration += 1


if __name__ == '__main__':
    main()