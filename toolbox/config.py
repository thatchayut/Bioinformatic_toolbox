# data required for using gene_based.py
gene_based = {

    # required information about the first dataset
    # a file containing mapping between probes IDs and sample
    "file_training_input_name" : "GSE2034-220712 (edited).csv",

    # number of rows of this file to be read
    "row_to_read" : "22071", 

    # name of a file containing mapping between samples and their class
    "file_training_output_name" : "mapping_sample_to_class_gse2034.csv",

    # required information to conduct an experiment
    # number of epochs
    "epoch" : "5",

    # number of folds
    "num_of_folds" : "10",

    # number of top-ranked features
    "number_of_ranked_gene" : "10"
}