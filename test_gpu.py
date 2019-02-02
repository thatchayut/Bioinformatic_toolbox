import numpy as np 
from timeit import default_timer as timer 
from numba import vectorize

@vectorize(["float32(float32, float32)"], target="cuda")
def vAdd(a,b):
    return a + b

def main():
    N = 32000000

    a = np.ones(N, dtype=np.float32)
    b = np.ones(N, dtype=np.float32)
    c = np.ones(N, dtype=np.float32)

    start = timer()
    c = vAdd(a, b)
    vAdd_time = timer() - start

    print(vAdd_time)

if __name__ == "__main__":
    main()