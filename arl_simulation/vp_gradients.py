"""Calculate gradients of voltage patterns. Also calculate expected noise level for standard sim.


"""

import logging
from data_models.parameters import arl_path
from processing_components.image.gradients import image_gradients
from processing_components.image.operations import export_image_to_fits, show_image, import_image_from_fits
import numpy

log = logging.getLogger(__name__)

real_vp = import_image_from_fits(arl_path('data/models/MID_GRASP_VP_real.fits'))
gradx, grady = image_gradients(real_vp)

gradxx, gradxy = image_gradients(gradx)
gradyx, gradyy = image_gradients(grady)

gradx.data *= real_vp.data
grady.data *= real_vp.data

gradxx.data *= real_vp.data
gradxy.data *= real_vp.data
gradyx.data *= real_vp.data
gradyy.data *= real_vp.data

import matplotlib.pyplot as plt

plt.clf()
fig, axs = plt.subplots(3, 2, gridspec_kw={'hspace': 0.3, 'wspace': 0})

axs[0,0].imshow(gradx.data[0,0])
axs[0,0].set_title("Gradx")
axs[0,0].axis('off')

axs[0,1].imshow(grady.data[0,0])
axs[0,1].set_title("Grady")
axs[0,1].axis('off')

axs[1,0].imshow(gradxx.data[0,0])
axs[1,0].set_title("Gradxx")
axs[1,0].axis('off')

axs[1,1].imshow(gradxy.data[0,0])
axs[1,1].set_title("Gradxy")
axs[1,1].axis('off')

axs[2,0].imshow(gradyx.data[0,0])
axs[2,0].set_title("Gradyx")
axs[2,0].axis('off')

axs[2,1].imshow(gradyy.data[0,0])
axs[2,1].set_title("Gradyy")
axs[2,1].axis('off')

plt.show()

plt.clf()
show_image(gradx, fig=fig, title='gradx')
plt.show()
plt.clf()
show_image(grady, title='grady')
plt.show()
plt.clf()
show_image(gradxx, title='gradxx')
plt.show()
plt.clf()
show_image(gradxy, title='gradxy')
plt.show()
plt.clf()
show_image(gradyx, title='gradyx')
plt.show()
plt.clf()
show_image(gradyy, title='gradyy')
plt.show()

export_image_to_fits(gradx, "MID_GRASP_gradients_gradx.fits" )
export_image_to_fits(grady, "MID_GRASP_gradients_grady.fits" )

export_image_to_fits(gradxx, "MID_GRASP_gradients_gradxx.fits" )
export_image_to_fits(gradxy, "MID_GRASP_gradients_gradxy.fits" )
export_image_to_fits(gradyx, "MID_GRASP_gradients_gradyx.fits" )
export_image_to_fits(gradyy, "MID_GRASP_gradients_gradyy.fits" )

cellsize = abs(real_vp.wcs.wcs.cdelt[0]) * numpy.pi / 180.0
print(cellsize)

nant= 197
nt = 65
s2r = numpy.pi / (3600.0 * 180)
pes = numpy.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0])

error = numpy.sqrt(2.0) * numpy.max(numpy.abs(gradx.data)) * pes * s2r / (cellsize * (nant - 1) * numpy.sqrt(nt))
print(error)

plt.clf()
plt.loglog(pes, error, '-', color='r')
plt.xlabel('Pointing error (arcsec)')
plt.ylabel('Error (Jy)')
plt.title('Predicted error due to pointing')
plt.savefig("prediction.png")
plt.show()
