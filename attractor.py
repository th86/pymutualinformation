import numpy as np
import pandas as pd
from mutualinfo import *

#from attractor import *
#df = pd.DataFrame(np.random.uniform(0,10,size=(30, 5)))
#findAttractor(df, df.iloc[:,0].values) 

def findAttractor(df, vec, a=5, bin=6, so=3, maxIter=100, epsilon=1E-14, verbose=True):
    c = 0
    mi = np.zeros(len(df.columns))

    for idx in range(len(df)):
        mi[idx] = mi2D(df.iloc[:,idx].values, vec, bin, so)

    premi = mi.copy()
    w = np.abs(mi)**a / sum(np.abs(mi)**a)
    w[mi < 0 ] = 0
    metagene = np.dot(df, w)
 
    while c < maxIter:
        for idx in range(len(df.columns)):
            mi[idx] = mi2D(df.iloc[:,idx].values, metagene, bin, so)

        delta = sum((mi - premi)**2)

        if verbose == True:
            print("\r{}, {}.".format(c, np.round(delta,4)), end="", flush=True)

        if delta < epsilon:
            break

        premi = mi.copy()
        w = np.abs(mi)**a / sum(np.abs(mi)**a)
        w[mi < 0] = 0
        metagene = np.dot(df, w)
    
        c = c + 1

    if c >= maxIter:
        return np.nan
    return mi 
