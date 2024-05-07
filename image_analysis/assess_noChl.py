#!/usr/bin/python
import sys
import numpy as np
from skimage.color import rgb2hsv
from skimage import io
import scipy.ndimage as ndimage
from PIL import Image, ImageDraw
import copy

class InfectionPhenotype:
    def __init__(self, image):
        # Initialize class variables
        self.img = image
        self.infection_area = 0
        self.leaf_area = 0
        self.infected_pixels_scaled = 0
        self.leaf_masked_scaled = 0

    def quantify_infection(self):
        # Convert image to HSV color space
        hsv = rgb2hsv(copy.deepcopy(self.img))
        
        # Create arrays for background, leaf, and infected pixels
        background = np.ones_like(hsv[:,:,0])
        leaf = np.ones_like(hsv[:,:,0])
        infected_pixels = np.ones_like(hsv[:,:,0])

        # Determine background pixels based on HSV values
        background_indices = np.nonzero(
            ((hsv[:, :, 0] > 190/360) & (hsv[:,:, 0] < 260/360) & (hsv[:,:, 2] > .15)))
        background[background_indices] = 0

        # Mask the leaf pixels using the background
        leaf[leaf > 0] = 1
        leaf_masked = np.multiply(copy.deepcopy(background), copy.deepcopy(leaf))

        # Dilate the leaf mask to handle shadows
        leaf_masked = ndimage.binary_dilation(leaf_masked, iterations=1)
        leaf_masked = 1 - leaf_masked

        # Determine infected pixels based on HSV values
        infection_indices = np.nonzero(
            (hsv[:, :, 0] < 48/255) | ((hsv[:,:,2] < .5) & (hsv[:,:,1] < .4)) | (hsv[:,:,2] < .4))
        infected_pixels[infection_indices] = 0
        infected_pixels_masked = np.multiply(copy.deepcopy(infected_pixels), copy.deepcopy(leaf_masked))
        
        # Filter out small infected regions
        infection_filter = ndimage.generic_filter(infected_pixels_masked, np.sum, size=3, mode='constant')
        infected_pixels_filtered = copy.deepcopy(infected_pixels_masked)
        infected_pixels_filtered[infection_filter < 5] = 0

        # Calculate areas and scale the masked images
        self.leaf_area = np.sum(leaf_masked)
        self.infection_area = np.sum(infected_pixels_masked)
        self.leaf_masked_scaled = (leaf_masked * 255).astype(np.uint8)
        self.infected_pixels_scaled = (infected_pixels_filtered * 255).astype(np.uint8)
        return self.infected_pixels_scaled

    def visualize_results(self, path):
        # Visualize the results by overlaying the masks on the original image
        height, width = self.img.shape[:2]
        
        leaf_rgb = np.zeros((height, width, 3), dtype=np.uint8)
        leaf_rgb[:,:,0] = self.leaf_masked_scaled
        leaf_rgba = np.dstack((leaf_rgb, self.leaf_masked_scaled))
        Image.fromarray(leaf_rgba).save(path + "_leaf.png")
        
        infection_rgb = np.zeros((height, width, 3), dtype=np.uint8)
        infection_rgb[:,:,0] = self.infected_pixels_scaled
        infection_rgba = np.dstack((infection_rgb, self.infected_pixels_scaled))
        Image.fromarray(infection_rgba).save(path + "_infection.png")

if __name__ == '__main__':
    print("Running...")
    # Parse command line arguments
    image_path = sys.argv[1] + "/" + sys.argv[2] + ".png"
    output_directory = sys.argv[3] + "/" + sys.argv[2]

    # Read input image
    image = io.imread(image_path)

    # Instantiate the InfectionPhenotype class
    phenotype = InfectionPhenotype(image)

    # Quantify infection and visualize results
    phenotype.quantify_infection()
    phenotype.visualize_results(output_directory)
