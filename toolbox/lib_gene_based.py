from scipy import stats
from sklearn.metrics import roc_auc_score
from copy import deepcopy

import math
import calculate
import random
import time
import add_ons
import pandas as pd
import numpy as np

def gene_based():
    # record start time
    start_time = time.time()
    
    