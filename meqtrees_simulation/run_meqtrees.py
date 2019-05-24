import os

i = 0
while os.path.exists('simulation%s'%i):
    i+=1
dirname = 'simulation%s/' % i
os.system('mkdir %s' % dirname)

print '******* making ideal visibilities *******'

os.system('casa -c pointing_errorsSKA.py make_ms %s' % dirname)

print '**** corrupting with pointing errors ****'

os.system('tigger-convert %ssource_list.txt %ssource_list.lsm.html' % (dirname, dirname))

errs = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
for i in errs:
    os.system('python run_HPBW.py %s %s' % (i, dirname)) # would rather do this and residuals on the fly so will rearrange

print '*********** getting residuals ***********'

os.system('casa -c pointing_errorsSKA.py get_residuals %s' % dirname)
