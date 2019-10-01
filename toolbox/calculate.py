import copy
import numpy as np
import pandas as pd
import math
from copy import deepcopy
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

def checkEqualListSize(list_1, list_2):
    check_valid = False
    num_of_chunks = None
    if (len(list_1) == len(list_2)):
        check_valid = True
        num_of_chunks = len(list_1)
    else:
        print(" WARNING : Number of chunks in the first dataset is not equal to Number of chunks in the second dataset.")
    return check_valid, num_of_chunks