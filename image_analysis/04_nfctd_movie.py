#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import copy
from skimage.color import rgb2gray, rgb2hsv
from skimage import io
import scipy.ndimage as ndimage
from plantcv import plantcv as pcv
from assess_noChl import *

# sys.argv[1] is directory of images
# sys.argv[2] is the image to be tested
# sys.argv[3] is the output directory

# Generate overlay masks for all images via assess. Usually used for movies. 
# CST & SE
# 2024_05_07

def vslz_rslts(img, quant_inf, path):
        fig, (test, orig) = plt.subplots(ncols=2, nrows=1, gridspec_kw = {'wspace':0, 'hspace':0})
        plt.subplots_adjust(wspace=None, hspace=None)
        fig.tight_layout()
        #subplots_adjust(left=None, bottom=None, right=None, top=None)
        orig.imshow(img)
        orig.axis('off')
        orig.text(0, 80, "Infected", fontsize=10, color='Red', ma='center')
        test.imshow(quant_inf, cmap='inferno', alpha=1)
        test.axis('off')
        plt.savefig(path + "_appended.png", transparent=True, dpi=558.8, bbox_inches='tight', pad_inches=0)
        plt.close()


if __name__ == '__main__':
    flnm = sys.argv[1] + "/" + sys.argv[2] + ".png"
    image = io.imread(flnm)  # generate array of color-values (RGB) by pixel
    phenotype= NfctnPhntyp("", image)
    phenotype.qntfy_nfctn()
    path =  sys.argv[3] + "/" + sys.argv[2]
    vslz_rslts(image, phenotype.nfctd_pxls_scaled , path)
# python '/home/user1/Downloads/navautron/leaf_hsv_1_3_movie_inf.py' '/home/user1/Downloads/navautron/2022.05.12/control/cropped' '1008' '/home/user1/Downloads/navautron'

