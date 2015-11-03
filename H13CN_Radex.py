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


#
maxiter = 100
debug   = False
radexpath = '~/Desktop/Radex/data/'
extension = '.dat'

pc = 3.08568025e18 # cm
pi = 3.14159265
k =  1.3806503e-16 # erg /K

tbg = 2.73   # cmb temp
dv = 0.4


fmt = "%7.2e cm^-2"

# Collision partners
#
# ph2co -> ph2, oh2
# oh2co -> ph2, oh2

moles = ['h13cn']
frequencies = [[517.9698,604.2679,690.5520]]

dir_input = "/Users/willwaalkes/Desktop/Radex/data/inp/"
dir_output = "/Users/willwaalkes/Desktop/Radex/data/out/"

############################################################

# this routine will write the input file for RADEX
def write_input(file,source,tkin,nh2,mole,cdmol,freqs):
    ne = 0.1*nh2   #electron density
    op = 3    #Ortho/Para ratio
    nph2 = nh2/(op+1)  #Total density is the sum of O and P densities
    noh2 = nph2*op   #getting O H2
    file.write(mole+'@xpol.dat\n') #THESE ALL WRITE A SINGLE LINE IN THE INPUT FILE
    file.write(dir_output+source+'-radex-'+mole+'.out\n') #
    file.write(str(freqs[0]-0.001)+' '+str(freqs[-1]+0.001)+'\n') 
    file.write(str(tkin)+'\n')
    file.write('2\n')
    file.write('oh2\n')
    file.write(str(noh2)+'\n')
    file.write('ph2\n')
    file.write(str(nph2)+'\n')
    file.write(str(tbg)+'\n')
    file.write(str(cdmol)+'\n')
    file.write(str(dv)+'\n')
 
# this reads the output created by RADEX and extracts the information you want: tkin,density,column density and flux 
def read_radex(results,freq):
    line  = results.readline()
    words = line.split()
    while (words[1] != "T(kin)"):
        line  = results.readline()
        words = line.split()
    tkin = words[-1] 
    line  = results.readline()
    words = line.split()
    dens = words[-1]
    line  = results.readline()
    words = line.split()
    while (words[1] != "Column"):
        line  = results.readline()
        words = line.split()
    cdmol = words[-1] 
    line  = results.readline()
    words = line.split()
    while (words[-1] != "FLUX"):
        line  = results.readline()
        words = line.split()
    line  = results.readline()
    line  = results.readline()
    words = line.split()
    ftmp  = float(words[4])
    while (freq!=ftmp):
        line  = results.readline()
        words = line.split()  
        ftmp  = float(words[4])
    TR = float(words[-5])
    #print freq,ftmp,TR
    return tkin,dens,cdmol,TR
 
# Begin of main program
def main(source):
    # Read Tdust and nH2 data 
    tgas = np.linspace()
    ngas = np.linspace()  # determine the H2 density from the H2 column density map. 
    cdens = np.linspace(1e16,1e19,50)
    cdmol =  tgas*0

    #tgas = tgas[0:10,0:10]

    file = open(dir_input+source+'-radex-'+mole+'.inp','w')
        for i in range(len(gastemp)):
            for j in range(len(gasdens)):
                for icd in range(ncd):
                    tgas = gastemp[i]
                    ngas = gasdens[j]
                    cd = cdens[icd]
                    write_input(file,source,tgas,ngas,mole,cd,freqs)
                    number_models += 1

                    if((i==tgas[-2]) & (icd==ncd-1)): 
                        file.write('0\n')
                        file.close()
                    else:
                        file.write('1\n')

        start = time.time()
         
        ## run the models
        print "Running ",number_models," models"
        os.system('radex < '+dir_input+source+'-radex-'+mole+'.inp > /dev/null')
      
        stop = time.time()
        duration = stop-start
        print "Run time = ",duration," seconds"

        #print "test"
        #os.sys.exit()
        # Read and save model results 

        for freq in freqs:

            tgas_radex = np.zeros(shape=(ncd,len(tgas)))
            ngas_radex = np.zeros(shape=(ncd,len(tgas)))
            cdmol_radex = np.zeros(shape=(ncd,len(tgas)))
            TR_radex = np.zeros(shape=(ncd,len(tgas)))

            print "Reading results for line "+mole+" "+str(freq)
            results = open(dir_output+source+'-radex-'+mole+'.out','r')
            for i in range(len(tgas)):
                for j in range(ncd):
                    tgas_radex[i,j],ngas_radex[i,j],cdmol_radex[i,j],TR_radex[i,j] = read_radex(results,freq)
                    #print tgas_radex[k,i,j],ngas_radex[k,i,j],cdmol_radex[k,i,j],TR_radex[k,i,j]

            #os.sys.exit()
                
            # update header            
            #hdr.update('NAXIS', 3)
            #hdr.update('NAXIS3', ncd)

            print TR_radex.shape

            #print hdr
            #os.sys.exit()

            # save results: FITS cube with the different column density maps
            pyfits.writeto("/Users/willwaalkes/Desktop/HCN_Research/Radex-Results/"+source+"-radex-"+mole+"-"+str(freq)+".fits", TR_radex,clobber=True)
            #pyfits.writeto("test.fits", TR_radex, clobber=True)

        results.close()

############################################################

main("Orion-S")
#main("Orion-KL")
#main("NGC-6334")
