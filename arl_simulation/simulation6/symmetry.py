import socket
import sys

import matplotlib as mpl
import numpy
import seqfile

mpl.use('Agg')

from matplotlib import pyplot as plt
from wrappers.serial.image.operations import show_image, import_image_from_fits, export_image_to_fits, \
    qa_image

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

pes = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]

results = []

for pe in pes:
    residual_image_plus = import_image_from_fits('simulation0/PE_%.1f_arcsec_arl.fits' % (pe))
    residual_image_minus = import_image_from_fits('simulation6/PE_%.1f_arcsec_arl.fits' % (pe))
    residual_image_plus.data += residual_image_minus.data
    
    print(qa_image(residual_image_plus, context='simulation6/PE_%.1f_arcsec_arl_difference.fits' % (pe)))
    show_image(residual_image_plus, title='simulation6/PE_%.1f_arcsec_arl_difference.fits' % (pe))
    plt.savefig('simulation6/PE_%.1f_arcsec_arl_difference.png' % (pe))
    plt.show()
    export_image_to_fits(residual_image_plus, 'simulation6/PE_%.1f_arcsec_arl_difference.fits' % (pe))
    
    nchan, npol, ny, nx = residual_image_minus.shape
    
    qa = qa_image(residual_image_plus)
    result = {}
    for field in ['maxabs', 'rms', 'medianabs']:
        result["onsource_" + field] = qa.data[field]
    result['onsource_abscentral'] = numpy.abs(residual_image_plus.data[0, 0, ny // 2, nx // 2])
    
    results.append(result)

import pprint

pp = pprint.PrettyPrinter()
pp.pprint(results)

plt.clf()
colors = ['b', 'r', 'g', 'y']
for ifield, field in enumerate(['onsource_maxabs', 'onsource_rms', 'onsource_medianabs', 'onsource_abscentral']):
    plt.loglog(pes, [result[field] for result in results], '-', label=field, color=colors[ifield])

plt.xlabel('Pointing error (arcsec)')
plt.ylabel('Error (Jy)')

title = 'simulation0 + simulation6'
plt.title(title)
plt.legend(fontsize='x-small')

plotfile = seqfile.findNextFile(prefix='simulation6/simulation0_plus_simulation6_%s_' % socket.gethostname(), suffix='.jpg')

print('Saving plot to %s' % plotfile)

plt.savefig(plotfile)
plt.show()
