import numpy as np

def first_order_error(t,M,Omega,N):
    
    error = (t*Omega*M)**3
    error = error/(N**2)
    error = error*np.exp((t*Omega*M)/N)
    
    return error

def second_order_error(t,M,Omega,N):
    
    error = (t*Omega)**3
    error *= M**2
    error = error/(N**2)
    error = error*np.exp((t*Omega*M)/N)
    
    return error

