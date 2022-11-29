from multiprocessing import Pool
import os

def f(x):
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    return x*x


if __name__ == '__main__':
    with Pool(2) as p:
        print(p.map(f, [1, 2, 3]))
