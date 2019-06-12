import numpy as np
from matplotlib import pyplot as plt


from numpy import genfromtxt



for i in range(4):
    pes = genfromtxt('arltim%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,19].astype(np.float) 
    print pes
    maxabs = genfromtxt('arltim%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,21].astype(np.float)    
    maxabs_MT = genfromtxt('errors%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,0].astype(np.float)
    if i ==0:
        plt.scatter(pes, maxabs, marker = 'x', s= 50,facecolors = 'none',edgecolors='red', label = 'ARL')
        plt.scatter(pes, maxabs_MT, marker = 'x', s= 50,facecolors = 'none',edgecolors='blue',label = 'MeqTrees')
    else:
        plt.scatter(pes, maxabs, marker = 'x', s= 50,facecolors = 'none',edgecolors='red')
        plt.scatter(pes, maxabs_MT, marker = 'x', s= 50,facecolors = 'none',edgecolors='blue')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('sigma offset (arcsec)')
    plt.ylabel('maxabs (Jy)')
    plt.xlim(0.8,300)
plt.legend(loc = 'upper left')

plt.savefig('maxabs_both.png')
plt.clf()

for i in range(4):
    pes = genfromtxt('arltim%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,19].astype(np.float) 
    print pes
    rms = genfromtxt('arltim%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,22].astype(np.float)    
    rms_MT = genfromtxt('errors%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,1].astype(np.float)
    if i ==0:
        plt.scatter(pes,rms, marker = 'x',s= 50,facecolors = 'none',edgecolors='red',label = 'ARL')
        plt.scatter(pes, rms_MT, marker = 'x', s= 50,facecolors = 'none',edgecolors='blue',label = 'MeqTrees')
    else:
        plt.scatter(pes,rms, marker = 'x',s= 50,facecolors = 'none',edgecolors='red')
        plt.scatter(pes, rms_MT, marker = 'x', s= 50,facecolors = 'none',edgecolors='blue')

    plt.xscale('log')
    plt.yscale('log')    
    plt.xlabel('sigma offset (arcsec)')
    plt.ylabel('rms (Jy)')
    plt.xlim(0.8,300)
    plt.ylim(1e-6,1e-3)
plt.legend(loc = 'upper left')
plt.savefig('rms_both.png')
plt.clf()


for i in range(4):
    pes = genfromtxt('arltim%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,19].astype(np.float) 
    print pes
    medianabs = genfromtxt('arltim%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,23].astype(np.float)    
    medianabs_MT = genfromtxt('errors%s.csv'%str(i+1), delimiter=',', dtype = None)[1:,2].astype(np.float)

    if i ==0:
        plt.scatter(pes, medianabs,marker = 'x', s= 50, facecolors = 'none',edgecolors='red',label = 'ARL')
        plt.scatter(pes, medianabs_MT, marker = 'x', s= 50,facecolors = 'none',edgecolors='blue',label = 'MeqTrees')
    else:
        plt.scatter(pes, medianabs,marker = 'x', s= 50, facecolors = 'none',edgecolors='red')
        plt.scatter(pes, medianabs_MT, marker = 'x', s= 50,facecolors = 'none',edgecolors='blue')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('sigma offset (arcsec)')
    plt.ylabel('medianabs (Jy)')
    plt.xlim(0.8,300)
    plt.ylim(1e-7,1e-3)
plt.legend(loc = 'upper left')

plt.savefig('medianabs_both.png')


