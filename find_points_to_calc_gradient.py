import numpy as np

def find_ppoints(x, y):
    """ todo"""
    x_new = []
    y_new = []
    for i in range(len(x)-1):
        if (y[i+1] - y[i])**2 > 0:
            a_inv_sqr = (x[i+1] - x[i])**2/(y[i+1] - y[i])**2 +1 
            a = 1/np.sqrt(a_inv_sqr)
            b = np.sqrt(1 - a**2)
            x_new.append(x[i] + a)
            y_new.append(y[i] + b)
        else:
            b = (y[i+1] -y[i])/abs(y[i+1] -y[i])
            x_new.append(x[i])
            y_new.append(y[i] + b)
    return x_new, y_new

    