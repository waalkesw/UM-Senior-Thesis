#! /usr/bin/python
#
import math
import os
import time
import sys
import numpy as np
from math import *
import pyfits
import pylab as pl
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
############################################################

def chi2(fmodel,fobs,ferr,fileout):
    Q = len(fobs)
    m = 3
    if Q <= m:
        m = 0 # free parameters
    else:
        m = 3
        
    tmp,hdr = pyfits.getdata("/Users/willwaalkes/Desktop/HCN_Research/Radex-Results/"+fmodel[0]+".fits", 0,header=True)
    #print "/Users/willwaalkes/Desktop/Radex/data/fits-radex-models-lines/"+fmodel[0]+".fits"
    #print tmp.shape
    ncd = tmp.shape[0]   
    model = np.zeros(shape=(len(fobs),tmp.shape[0],tmp.shape[1],tmp.shape[2]))
    obs = fobs
    err = ferr
    
    gastemp = np.linspace(60,220,30)
    gasdens = np.logspace(2,8,30)
    cdens = np.logspace(14,15,30)
    ncd = len(cdens)
    ntmp = len(gastemp)
    ndens = len(gasdens)

    for i in range(len(fobs)):
        model[i,:,:,:] = pyfits.getdata("/Users/willwaalkes/Desktop/HCN_Research/Radex-Results/"+fmodel[i]+".fits", 0)

    chi2 = np.zeros(shape=(len(fobs),model.shape[1],model.shape[2],model.shape[3]))
    
    N = np.array(range(len(fobs)))
    n = np.array(range(len(fobs))) 
    T = np.array(range(len(fobs))) 
                    
    for iobs in range(len(fobs)):
        observed = obs[iobs]
        expected = model[iobs,:,:,:]
        sigma = err[iobs]**2
        chi2[iobs,:,:,:] = (observed-expected)**2/sigma
        #Best_Model = model[np.where()]
    chi2_sum = np.nansum(chi2,axis=0)/(Q-m)
    print chi2_sum.shape
    
    #NP.NANMIN IF NANS CARRY OVER        

    chi2_min = np.amin(chi2_sum)
    print chi2_min
    indices = np.where(chi2_sum == chi2_min)
    print indices
    N = cdens[indices[0]]
    T = gastemp[indices[1]]
    n = gasdens[indices[2]]
    print N,T,n
    #print chi2_min.shape   
                
    for i in range(ncd):      
        plt.imshow(chi2_sum[i,:,:],
                   extent=(gastemp[0],gastemp[-1],np.log10(gasdens[0]),np.log10(gasdens[-1])),
                    aspect='auto',cmap='coolwarm',
                    vmin=500,vmax=1000)
        plt.xlabel('Temp')
        plt.ylabel('Log (gasdens)')
        plt.title(i)
        plt.show()
        plt.clf()
      
    #PLOT THE MODEL VS OBSERVATION
    Freq = [531.7163,620.3040,708.8770,797.4333,885.9707,974.4872,1062.9807]
    plt.scatter(Freq,obs,s=7,marker='o',label = 'Observation',color='r')
    plt.errorbar(Freq,obs,yerr=err[:],fmt='none',ecolor='k')
    plt.scatter(Freq,model[:,indices[0],indices[1],indices[2]],s=7,marker='o',label = 'Model',color='b')
    plt.legend(loc='upper right')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Intensity (K*Km/s)')
    plt.title('Model vs Observations')
    plt.savefig(filename='Quality-of-fit.eps')
    plt.show()
    
    pyfits.writeto("/Users/willwaalkes/Desktop/HCN_Research/N-Best-Fit/"+fileout+"-chi2.fits", chi2_min, hdr, clobber=True)

    pyfits.writeto("/Users/willwaalkes/Desktop/HCN_Research/N-Best-Fit/"+fileout+".fits", N, hdr, clobber=True)
  #  pyfits.writeto("/Users/willwaalkes/Desktop/Radex/data/N-best-fit/"+fileout+"-interp.fits", N_interp, hdr, clobber=True)
 
############################################################

def HCN(sou,source):
    model = [sou+"-radex-hcn-531.7163",sou+"-radex-hcn-620.304",sou+"-radex-hcn-708.877",
    sou+"-radex-hcn-797.4333",sou+"-radex-hcn-885.9707",sou+"-radex-hcn-974.4872",
    sou+"-radex-hcn-1062.9807"]
    obs = [63.649,49.284,25.267,20.235,29.377,6.5874,5.925] #in K*km/s
    err = [0.361,0.346,0.272,0.254,0.354,0.301,0.563]    
    fileout = source+"-hcn"
    chi2(model,obs,err,fileout)

#############################################################

def H13CN(sou,source):
    model = [sou+"-radex-h13cn-517.9698",sou+"-radex-h13cn-604.2679",sou+"-radex-h13cn-690.2251"]
    obs = [2.8916,1.9734,1.7245]
    err = [0.076,0.07,0.115]
    fileout = source+"-h13cn"
    chi2(model,obs,err,fileout)
    
#############################################################
    
def HC15N(sou,source):
    model = [sou+"-radex-hc15n-516.2606"]
    obs = [1.1093]
    err = [0.1]
    fileout = source+"-N-hc15n"
    chi2(model,obs,err,fileout)
    
#############################################################

#Orion South

HCN("Orion-S","Orion-S")
#H13CN("Orion-S","Orion-S")
#HC15N("Orion-S","Orion-S")

#Orion KL

#HCN("Orion-KL","Orion-KL")
#H13CN("Orion-KL","Orion-KL")
#HC15N("Orion-KL","Orion-KL")

#NGC 6334

#HCN("NGC-6334","NGC-6334")
#H13CN("NGC-6334","NGC-6334")
#HC15N("NGC-6334","NGC-6334")