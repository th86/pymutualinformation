import numpy as np
import pandas as pd
from mutualinfo import *

#from attractor import *
#df = pd.DataFrame(np.random.uniform(0,10,size=(10, 5)))
#findAttractor(df, df.iloc[1,:].values)

def findAttractor(df, vec, a=5, bin=6, so=3, maxIter=100, epsilon=1E-14, verbose=True):
    c = 0
    mi = np.zeros(len(df))

    for idx in range(len(df)):
        mi[idx] = mi2D(df.iloc[idx,:].values, vec, bin, so)

    premi = mi.copy()
    w = np.abs(mi)**a / sum(np.abs(mi)**a)
    w[mi < 0 ] = 0
    metagene = np.dot(w, df)
 
    while c < maxIter:
        for idx in range(len(df)):
            mi[idx] = mi2D(df.iloc[idx,:].values, metagene, bin, so)

        delta = sum((mi - premi)**2)

        if verbose == True:
            print(c, ":", np.round(delta,4) )

        if delta < epsilon:
            break

        premi = mi.copy()
        w = np.abs(mi)**a / sum(np.abs(mi)**a)
        w[mi < 0] = 0
        metagene = np.dot(w, df)
    
        c = c + 1

    if c >= maxIter:
        return NaN
    return mi 