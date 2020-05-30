import numpy as np
from functools import reduce
from multiprocessing import Pool
from itertools import product
from math import sqrt

def add(x, y):
    return x + y

def subtract(x, y):
    return x - y

def square(x):
    return x*x

def sum_division(a, b, c, d):
    """Computes (a + b) / (c + d)"""
    return (a + b) / (c + d)

def mean(somelist):
    return reduce(add, somelist) / len(somelist)

def sd(somelist):
    return sqrt(reduce(add, map(square, list(map(subtract, somelist, [mean(somelist)]*len(somelist))))) / (len(somelist) - 1))


def dot_product(i, j):
    return reduce(add, map(prod, A[i, :], B[:, j]))

def prod(x, y):
    return x * y


n = 10
np.random.seed(123)
A = np.random.randint(-10, 10, size=(n, n))
B = np.random.randint(-10, 10, size=(n, n))


if __name__ == "__main__":
    # Store here all the futures
    futures = []
    
    with Pool(processes=5) as pool:
        # These two operations require only one argument so map_async is enough
        futures.append(pool.apply_async(add, [1, 2]))
        futures.append(pool.apply_async(subtract, [1, 2]))
        futures.append(pool.apply_async(sum_division, [10, 2, 3, 4]))
        futures.append(pool.apply_async(square, [20])) 
        futures.append(pool.apply_async(mean, [[1, 2, 3]]))
        futures.append(pool.apply_async(sd, [[1, 2, 3]]))
        # Map + Asynchronous function
        futures.append(pool.starmap_async(dot_product, product(np.arange(n), np.arange(n)), chunksize=50))
        
        ready = [False] * len(futures)
        while not all(ready):
            for i, future in enumerate(futures):
                if future.ready():
                    ready.pop(i)
                    futures.pop(i)
                    if future.successful():
                        print("Process finished. Result {result}".format(result=future.get()))
                    else:
                        print("Process finished, but unsuccessful.")
                        try:
                            future.get()
                        except Exception as exception:
                            print("Error: {exception_type} : {exception_raised}".format(exception_type=type(exception), exception_raised=exception))
                    
                else:
                    next
