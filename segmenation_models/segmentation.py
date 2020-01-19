import segmentation_models as sm
from os import sys

sys.path.append('../')

from pharynx_redox import (
    pharynx_io as pio,
    image_processing as ip,
    experiment,
    plots,
    profile_processing as pp,
    utils,
    data_analysis as da,
    constants
)

import numpy as np
import matplotlib.pyplot as plt
from skimage import measure, transform

imgs = pio.load_tiff_from_disk('/home/niketrp/apfeld/pharynx_redox/data/segmentation_training_data/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47-410-0.tif')
segs = pio.load_tiff_from_disk('/home/niketrp/apfeld/pharynx_redox/data/segmentation_training_data/2017_02_22-HD233_SAY47/seg-2017_02_22-HD233_SAY47-410-0.tif')

def crop(s_img, fl_img, crop_amount=32):
    # Padding the images with 32 cells down and left
    s_img = np.pad(s_img, ((0,32),(32,0)), 'constant')
    fl_img = np.pad(fl_img, ((0,32),(32,0)), 'constant')
   
    # Getting the center of pharynx based on padded image
    l = measure.label(s_img)
    props = measure.regionprops(l)
    center = props[0].centroid
    
    # Cropping images around center
    upper_bound = int(center[0] - crop_amount)
    lower_bound = int(center[0] + crop_amount)
    left_bound = int(center[1] - crop_amount)
    right_bound = int(center[1] + crop_amount)
    s_img_cropped = s_img[upper_bound:lower_bound, left_bound:right_bound]
    f_img_cropped = fl_img[upper_bound:lower_bound, left_bound:right_bound]
    
    cropped_imgs = [s_img_cropped, f_img_cropped]
    return cropped_imgs

all_cropped = []
for i in range(0, len(segs) - 1):
    S = segs[i]
    fl = imgs[i]
    cropped = crop(S,fl)
    all_cropped.append(cropped)

plt.imshow(all_cropped[119][1])

BACKBONE = 'resnet34'
preprocess_input = sm.get_preprocessing(BACKBONE)

# load your data
# x_train, y_train, x_val, y_val = load_data(all_cropped)

# # preprocess input
# x_train = preprocess_input(x_train)
# x_val = preprocess_input(x_val)

# # define model
# model = sm.Unet(BACKBONE, encoder_weights='imagenet')
# model.compile(
#     'Adam',
#     loss=sm.losses.bce_jaccard_loss,
#     metrics=[sm.metrics.iou_score],
# )

# # fit model
# # if you use data generator use model.fit_generator(...) instead of model.fit(...)
# # more about `fit_generator` here: https://keras.io/models/sequential/#fit_generator
# model.fit(
#    x=x_train,
#    y=y_train,
#    batch_size=16,
#    epochs=100,
#    validation_data=(x_val, y_val),
# )