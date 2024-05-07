#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
import copy
from skimage.color import rgb2gray, rgb2hsv
from plantcv import plantcv as pcv
from skimage import io
import scipy.ndimage as ndimage
from PIL import Image, ImageDraw, ImageFilter


# sys.argv[1] is directory of images
# sys.argv[2] is the image to be tested
# sys.argv[3] is the output directory

# Perform HSV-based thresholding on single images 
# modify line 57 to match observations
# Severin Einspanier
# 07.05.2024

class NfctnPhntyp:
    def __init__(self, name, image):
        self.nfst_name = name  # only required for infest
        self.img = image  # image from the infest.py script
        self.nfctn_area = 0  # area of infected pixels
        self.lflt_area = 0  # area of leaf pixels
        self.bckgrnd_scaled = 0  # background pixels scaled to 255
        self.lflt_masked_scaled = 0  # leaf pixels scaled to 255
        self.nfctd_pxls_scaled = 0  # infected pixels scaled to 255

    def qntfy_nfctn(self):
        # load images
        hsv = rgb2hsv(copy.copy(self.img))
        bckgrnd = rgb2gray(copy.copy(self.img))
        lflt = rgb2gray(copy.copy(self.img))
        nfctd_pxls = rgb2gray(copy.copy(self.img))
        
        # Background removed if blue and not Black (important to keep dark lesions in)
        
        bckgrnd_index = np.nonzero(((hsv[:, :, 0] > 190/360) & (hsv[:,:, 0] < 260/360)  & (hsv[:,:, 2] > .15)))
        
          
        bckgrnd[bckgrnd_index] = 1
        bckgrnd[bckgrnd < 1] = 0
        
        # Use background as a mask
        lflt[lflt > 0] = 1
        lflt_masked = np.multiply(copy.copy(bckgrnd), copy.copy(lflt))

        ## THE FOLLOWING PART is modified by SE. Idea: removing the shadows by increasing the background.
        # So the bg will be dilated to make it a little bit bigger resulting in overwriting the shadows. 
        
        lflt_masked = pcv.dilate(lflt_masked, ksize=9, i=1)
        lflt_masked = 1 - lflt_masked 

        # determine infected pixels

        nfctn_index = np.nonzero((hsv[:, :, 0] < 48/255)|((hsv[:,:,2]<.5)&(hsv[:,:,1]<.4))|(hsv[:,:,2]<.4))
        #

        nfctd_pxls[nfctn_index] = 1
        nfctd_pxls[nfctd_pxls < 1] = 0
        nfctd_pxls_masked = np.multiply(copy.copy(nfctd_pxls), copy.copy(lflt_masked))
        #import pdb; pdb.set_trace()
        inf_fltr = ndimage.generic_filter(nfctd_pxls_masked, sum, size = 3, mode = 'constant')
        nfctd_pxls_fltrd = copy.copy(nfctd_pxls_masked)
        nfctd_pxls_fltrd[inf_fltr < 5] = 0
 
        # Output
        bckgrnd = np.absolute(bckgrnd - 1)  # invert colors
        self.nfctn_area = np.sum(nfctd_pxls_masked)
        self.lflt_area = np.sum(lflt_masked)
        self.lflt_masked_scaled = lflt_masked.astype(np.uint8) * 255
        self.nfctd_pxls_scaled = nfctd_pxls_fltrd.astype(np.uint8) * 255
        return self.nfctd_pxls_scaled
        print(self.nfctn_area)
        print(self.lflt_area)

    def vslz_rslts(self, path):
        hght, wdth = self.img.shape[:2]
        ##
        lflt_rgb = np.array(copy.copy(self.img))
        lflt_rgb[:, :, 0] = copy.copy(self.lflt_masked_scaled)
        lflt_rgb[:, :, 1] = np.zeros((hght, wdth), dtype=np.uint8)
        lflt_rgb[:, :, 2] = np.zeros((hght, wdth), dtype=np.uint8)
        lflt_rgba = np.dstack((lflt_rgb, copy.copy(self.lflt_masked_scaled)))
        Image.fromarray(lflt_rgba).save(path + "_lflt.png")
        ##
        nfctn_rgb = np.array(copy.copy(self.img))
        nfctn_rgb[:, :, 0] = copy.copy(self.nfctd_pxls_scaled)
        nfctn_rgb[:, :, 1] = np.zeros((hght, wdth), dtype=np.uint8)
        nfctn_rgb[:, :, 2] = np.zeros((hght, wdth), dtype=np.uint8)
        nfctn_rgba = np.dstack((nfctn_rgb, copy.copy(self.nfctd_pxls_scaled)))
        Image.fromarray(nfctn_rgba).save(path + "_nfctn.png")


if __name__ == '__main__':
    print("runs")
    flnm = sys.argv[1] + "/" + sys.argv[2] + ".png"
    name = ""  # only required for infest
    image = io.imread(flnm)  # generate array of color-values (RGB) by pixel
    inf_phe = NfctnPhntyp(name, image)
    inf_phe.qntfy_nfctn()
    inf_phe.vslz_rslts(sys.argv[3] + "/" + sys.argv[2])
# python '/home/christr/Downloads/navautron/leaf_hsv_1_3.py' '/home/christr/Desktop/test' '0522' '/home/christr/Desktop/test'
# python3 '/home/se_residevo/Desktop/SEVERIN LAB/code/2022_08_22assess.py' '//home/se_residevo/Desktop/2022_09_14EVAL' '0007' '/home/se_residevo/Desktop/2022_09_14EVAL/output'



