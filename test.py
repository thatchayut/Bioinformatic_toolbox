import random
import time
import numpy as np
import random

def main():
    matrix = np.matrix('1 2; 3 4 ; 5 6')
    print(matrix)
    print()
    result = 2 * matrix
    print(result)

    print(result.shape)

    noise = (0.00001 * (np.random.rand(result.shape[0], result.shape[1])))
    print(noise)
    print(type(noise))
    matrix_noise = np.matrix(noise)
    print(matrix_noise)
    print(type(matrix_noise))

    print(result + noise)
if __name__ == '__main__':
    main()