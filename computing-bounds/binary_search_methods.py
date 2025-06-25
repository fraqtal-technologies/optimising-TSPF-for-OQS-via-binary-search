import numpy as np

def binarySearch(low,high,erorr_function,epsilon,**kwargs):
    high = high
    low = 1
    #mid = int(np.ceil((high+low)/2))
    mid = (high+low)//2
    error = erorr_function(**kwargs)

    while (low<high and error!=epsilon):
        error = erorr_function(**kwargs)
        if error > epsilon:
            low = mid+1
            mid = (high+low)//2
            #mid = int(np.ceil((high+low)/2))
        elif error < epsilon:
            high = mid-1
            mid = (high+low)//2
            #mid = int(np.ceil((high+low)/2))
    #return mid+1
    return mid