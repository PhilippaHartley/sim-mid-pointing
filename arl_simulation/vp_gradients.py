"""Unit tests for testing support


"""

import logging
from data_models.parameters import arl_path
from processing_components.image.gradients import image_gradients
from processing_components.image.operations import export_image_to_fits, show_image, import_image_from_fits

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

exit()

plt.clf()
show_image(gradx, fig=fig, ax=axs[0,0], title='gradx')
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
