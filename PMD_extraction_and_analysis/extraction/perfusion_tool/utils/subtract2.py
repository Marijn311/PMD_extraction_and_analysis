import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
from matplotlib.widgets import Slider

path_new = r"C:\Users\20192010\OneDrive - TU Eindhoven\Junior Researcher TUe\data\test\SSR-INS-04711\20181028_5000670\perf_nii\svd\CBV_lc_original.nii.gz"
path_old = r"C:\Users\20192010\OneDrive - TU Eindhoven\Junior Researcher TUe\data\test\SSR-INS-04711\20181028_5000670\perf_nii\svd\CBV_original.nii.gz"

# Load the images
img_new = nib.load(path_new)
img_old = nib.load(path_old)

# Get the data
data_new = img_new.get_fdata()
data_old = img_old.get_fdata()


# print("data_new at index 120,120,10 ", data_new[120,120,10])
# print("data_old at index 120,120,10 ", data_old[120,120,10])
# print("data_diff at index 120,120,10 ", data_new[120,120,10] - data_old[120,120,10])
# print("this difference as a factor is ", (data_new[120,120,10] - data_old[120,120,10])/data_old[120,120,10])

# # do the same for 130,130,10
# print("data_new at index 130,130,10 ", data_new[130,130,10])
# print("data_old at index 130,130,10 ", data_old[130,130,10])
# print("data_diff at index 130,130,10 ", data_new[130,130,10] - data_old[130,130,10])
# print("this difference as a factor is ", (data_new[130,130,10] - data_old[130,130,10])/data_old[130,130,10])




# Subtract the images
data_diff = data_new - data_old


# Function to update the plot
def update(val):
    slice_idx = int(slider.val)
    ax.imshow(data_diff[:, :, slice_idx], cmap='gray')
    fig.canvas.draw_idle()

# Plot the initial slice
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
slice_idx = data_diff.shape[2] // 2
img = ax.imshow(data_diff[:, :, slice_idx], cmap='gray')

# Add a slider
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Slice', 0, data_diff.shape[2] - 1, valinit=slice_idx, valstep=1)
slider.on_changed(update)

plt.show()
