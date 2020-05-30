from scoop import futures
from functools import reduce
import numpy as np


def prod_modulo(a, b, c, d):
    """Computes (a + b) / (c + d)"""
    return (a % 7) * (b % 7) * (c % 7) * (d % 7)

def subtract(x, y):
    """Computes x - y"""
    return x - y

if __name__ == "__main__":

    a = np.random.randint(-10, 10, size=(100, ))
    b = np.random.randint(-10, 10, size=(100, ))
    c = np.random.randint(-10, 10, size=(100, ))
    d = np.random.randint(-10, 10, size=(100, ))

    out = futures.mapReduce(prod_modulo, subtract, a, b, c, d)

    print("Result: ", out)