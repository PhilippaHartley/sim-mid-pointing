import sys, numpy as np, os
from Timba.Apps import meqserver
from Timba.TDL import Compile
from Timba.TDL import TDLOptions
import run_HPBW
from run_HPBW import hdf2npy

i = 0
while os.path.exists('simulation%s'%i):
    i+=1
dirname = 'simulation%s/'%i
os.system('mkdir '+dirname)

print '******* making ideal visibilities *******'

os.system('casa -c pointing_errorsSKA.py make_ms %s'%dirname)

print '**** corrupting with pointing errors ****'


errs = [0.0,1.0,2.0,4.0,8.0,16.0,32.0,64.0,128.0,256.0]
for i in errs:
    #os.system('casa -c pointing_errorsSKA.py get_offsets %s %s'%(dirname,i))
    os.system('python run_HPBW_using_RandomNoise.py %s %s'%(i,dirname)) # would rather do this and residuals on the fly so will rearrange in future

print '*********** getting residuals ***********'

os.system('casa -c pointing_errorsSKA.py get_residuals %s'%dirname)



