import numpy as np
from functools import reduce
from multiprocessing import Pool
from itertools import product

def add(x, y):
    return x + y

def prod(x, y):
    return x * y

def dot_product(i, j):
    return reduce(add, map(prod, A[i, :], B[:, j]))

n = 100
np.random.seed(123)
A = np.random.randint(-10, 10, size=(n, n))
B = np.random.randint(-10, 10, size=(n, n))


if __name__ == "__main__":
    with Pool(processes=5) as pool:
        result = pool.starmap(dot_product, product(np.arange(n), np.arange(n)))
      
    print(np.reshape(list(result), (n, n)))
    