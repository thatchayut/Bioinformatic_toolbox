import random
import time
import numpy as np
import random

def main():
    start_time = time.time()
    for i in range(0, 10000000):
        i +=1    
    end_time = time.time()
    print("start_time : " + str(start_time))
    print("end_time : " + str(end_time))
    time_elapse = end_time - start_time
    print("time_elapse second : " + str(time_elapse))
    time_elapse_min = time_elapse / 60
    print("time_elapse minute : " +str(time_elapse_min))
if __name__ == '__main__':
    main()