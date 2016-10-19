from __future__ import division #for floating number division in python2
import numpy as np

#Taihsien Ou Yang, 2016
#References:
#Carsten O Daub et al., Estimating mutual information using B-spline functions an improved similarity measure for analysing gene expression data
#BMC Bioinformatics, 2004
#Weiyi Cheng et al., https://github.com/weiyi-bitw/cafr, 2012

#Cox-de Boor recursion formula
def basisFunction( i, p, t, kVector, numBins):
  if p==1:
    
    if ((t>=kVector[i]) and (t< kVector[i+1]) and (kVector[i] < kVector[i-1])):
      return 1;
    elif abs(t-kVector[i+1])<1e-10 and (i+1 == numBins):
      return 1;
    else:
      return 0;
  else:

    d1 = kVector[i+p-1] - kVector[i]; #u(i+p)-ui
    n1 = t-kVector[i];                #u-ui
    d2 = kVector[i+p] - kVector[i+1]; #u(i+p+1)-u(i+1)
    n2 = kVector[i+p] - t;            #u(i+p+1)-u
    
    if d1<1e-10 and d2<1e-10: 
      e1=0;
      e2=0;
    elif d1<1e-10:
      e1=0
      e2=n2/d2*basisFunction(i+1, p-1, t, kVector, numBins); #N(i+1,p-1)(u)
    elif d2<1e-10:
      e1=n1/d1*basisFunction(i, p-1, t, kVector, numBins); #N(i,p-1)(u)
      e2=0;
    else:
      e1=n1/d1*basisFunction(i, p-1, t, kVector, numBins); #N(i,p-1)(u)
      e2=n2/d2*basisFunction(i+1, p-1, t, kVector, numBins); #N(i+1,p-1)(u)

    if e1+e1<0 :
      return 0;

  return e1+e2;

#Convert the values into B-Spline weights
def findWeights(x, knots, numBins, splineOrder): #only takes values within 0-1
  numSamples = len(x)
  weights = np.zeros(numSamples*numBins)
  for curSample in range(0,(numSamples-1)):
    for curBin in range(0,(numBins-1)):
      weights[curBin*numSamples+curSample]=basisFunction(curBin, splineOrder, x[curSample], knots, numBins);

  return weights;

#1-D Entropy
def entropy1D(weights, numSamples, numBins ):
  H = 0;
  for i in range(0, (numBins-1)): #compute the entropy of every bin using the B-spline-based probabilities.
    p = 0;
    for j in range(0,(numSamples-1)): #go through every sample
      p = p + weights[i*numSamples + j];
    if p>0:
      H=H-(p/numSamples)*np.log2(p/numSamples);

  return H;

#2-D Entropy 
def entropy2D(weightsx, weightsy, numSamples, numBins):
  H = 0;
  for ix in range(0, (numBins-1)):
    for iy in range(0, (numBins-1)):
      p = 0;
      for j in range(0,(numSamples-1)):
        p = p + weightsx[ix*numSamples + j] * weightsy[iy*numSamples + j];
      if p>0:
        H=H-(p/numSamples)*np.log2(p/numSamples);

  return H;

#Mutual Information
def mi2D(x,y,numBins,splineOrder):
  numSamples=len(x);
  u=np.zeros(numBins+splineOrder);

  #Put default weights in u, so that u[0:(bin+so-1)]= _____/-----
  nInternalPoints = numBins - splineOrder; 
  for i in range(splineOrder, (splineOrder+nInternalPoints-1) ):
    u[i]=(i-splineOrder+1)/(nInternalPoints+1);
  for i in range(splineOrder+nInternalPoints, (numBins+splineOrder-1) ):
    u[i]=1;

  wx=findWeights(x, u, numBins, splineOrder);
  wy=findWeights(y, u, numBins, splineOrder);

  mi= entropy1D(wx, numSamples, numBins ) + entropy1D(wy, numSamples, numBins ) - entropy2D(wx ,wy, numSamples, numBins);

  return mi;