import logging
import sys

import numpy

from processing_library.image.operations import create_empty_image_like
from wrappers.serial.image.operations import export_image_to_fits, import_image_from_fits

import matplotlib.pyplot as plt

log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(logging.StreamHandler(sys.stdout))
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

import pprint

pp = pprint.PrettyPrinter()

from scipy import interpolate

# x = np.arange(0, 10)
# y = np.exp(-x/3.0)
# f = interpolate.interp1d(x, y)
#
# xnew = np.arange(0,9, 0.1)
# ynew = f(xnew)   # use interpolation function returned by `interp1d`
# plt.plot(x, y, 'o', xnew, ynew, '-')
# plt.show()

elevations_in = numpy.array([15, 45, 90], dtype='float')
elevations_out = numpy.array([15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90], dtype='float')
elevations_out = numpy.arange(15.0, 90, 1.0)
default = 1
nchan = 1
npol = 4
ny = 1024
nx = 1024

array_in = numpy.zeros([nchan, npol, ny, ny, len(elevations_in)])
array_out = numpy.zeros([nchan, npol, ny, ny, len(elevations_out)])

im_in = "../Ku_{el:d}_11700_{type}.fits"
im_out = "Ku_{el:d}_11700_{type}_interpolated.fits"
im_diff_out = "Ku_{el:d}_11700_{type}_interpolated_difference.fits"

im_template = None

for type in ['real', 'imag']:
    for iel, el in enumerate(elevations_in):
        print("Reading elevation %s part elevation %.0f" % (type, el))
        im_in_file = im_in.format(el=int(el), type=type)
        im = import_image_from_fits(im_in_file)
        array_in[..., iel] = im.data
        if im_template is None:
            im_template = create_empty_image_like(im)
    
    f = interpolate.interp1d(elevations_in, array_in, axis=4, kind='quadratic')
    array_out = f(elevations_out)

    rms_vp = []
    max_vp = []
    min_vp = []
    rms_diff = []
    max_diff = []
    min_diff = []


    for iel, el in enumerate(elevations_out):
        print("Writing elevation %s part %.0f" % (type, el))
        im_template.data = array_out[..., iel]
        im_out_file = im_out.format(el=int(el), type=type)
        export_image_to_fits(im_template, im_out_file)
        rms_vp.append(numpy.std(im_template.data[0,0:1,...]))
        max_vp.append(numpy.max(im_template.data[0,0:1,...]))
        min_vp.append(numpy.min(im_template.data[0,0:1,...]))
        im_template.data -= array_in[..., default]
        im_diff_out_file = im_diff_out.format(el=int(el), type=type)
        export_image_to_fits(im_template, im_diff_out_file)
        rms_diff.append(numpy.std(im_template.data[0,0:1,...]))
        max_diff.append(numpy.max(im_template.data[0,0:1,...]))
        min_diff.append(numpy.min(im_template.data[0,0:1,...]))

    plt.clf()
    plt.plot(elevations_out, rms_vp, '-', color='r', label='VP rms')
    if type == 'imag':
        plt.plot(elevations_out, max_vp, '.', color='g', label='VP max')
    plt.plot(elevations_out, min_vp, '-', color='b', label='VP min')
    plt.plot(elevations_out, rms_diff, '.', color='r', label='VP diff rms')
    plt.plot(elevations_out, max_diff, '.', color='g', label='VP diff max')
    plt.plot(elevations_out, min_diff, '.', color='b', label='VP diff min')
    plt.xlabel('Elevation')
    plt.ylabel('Value')
    plt.title('Statistics in %s part of 11700MHz voltage pattern' % type)
    plt.legend()
    plt.savefig('%s_vp_statistics.png' % type)
    plt.show(block=False)