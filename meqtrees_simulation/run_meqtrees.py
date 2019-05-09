import sys, numpy as np, os
from Timba.Apps import meqserver
from Timba.TDL import Compile
from Timba.TDL import TDLOptions
import run_HPBW
from run_HPBW import hdf2npy

print '******* making ideal visibilities *******'
os.system('casa -c pointing_errorsSKA.py make_ms')

print '**** corrupting with pointing errors ****'


errs = [0.0, 1.0,2.0]#,4.0,8.0,16.0,32.0,64.0,128.0,256.0]
for i in errs:
    os.system('python run_HPBW.py %s'%i ) # would rather do this and residuals on the fly so will rearrange

print '*********** getting residuals ***********'

os.system('casa -c pointing_errorsSKA.py get_residuals')



