import random
import time
def main():
    random.seed(1)
    for i in range (0, 100):
        # random.seed(time.time())
        a = random.randint(0,4)
        print(a)
        print("test")

if __name__ == '__main__':
    main()