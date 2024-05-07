# Parser for assess, takes all images and ROIs and performs 'assess.py' only there. 
# Severin Einspanier 
# 2024_05_07

import numpy as np
import sys
from skimage.color import rgb2gray
from skimage.measure import regionprops
from skimage import io
from scipy.cluster.vq import *
from assess_noChl import *
import os
import datetime
import argparse

class Panel:
    def __init__(self, image, path, N):
        # Initialize panel properties
        self.N = N
        self.path = path
        self.i_original = image
        self.i_rgb_thre = self.remove_background()
        self.leaf_stack = []
        self.i_file = rgb2gray(self.i_rgb_thre)
        self.Nx = 0
        self.Ny = 0
        self.exist_layout = self.test_layout()
        if not self.exist_layout:
            sys.exit()
        else:
            self.order_bb1()

    def test_layout(self):
        # Check if the grid layout file exists
        out = False
        if os.path.isfile(args.grid):
            out = True
        if not out:
            print("\nNo grid layout found!\nPlease create a grid layout")
            sys.exit()
        return out

    def get_mean(self):
        # Calculate mean and standard deviation of region areas
        tab = []
        for region in regionprops(self.label_image):
            tab.append(region.area)
        return np.mean(tab), np.std(tab)

    def get_layout(self):
        # Extract panel layout information from layout file
        self.m_print("Finding panel layout", 0)
        found = False
        layout, tab = [], []
        for fi in os.listdir(self.path):
            if fi.endswith(".thelay"):
                f = open(self.path+fi, 'r')
                for line in f.readlines():
                    tab = line.split()
                    self.Nx = len(tab)
                    for n in tab:
                        layout.append(n)
                    self.Ny += 1
                f.close()
                found = True
        if not found:
            print("Layout file not found!", 0)
            sys.exit(0)
        self.layout = layout

    def find_grid(self, inp):
        # Find grid of leaf regions based on coordinates
        dx, dy = [], []
        outx, outy = [], []
        for i in inp:
            dy.append([i[2], 0])
        data = np.vstack(dy)
        centroids, _ = kmeans(data, self.Ny)
        idx, _ = vq(data, centroids)
        for ii in range(0, self.Ny):
            mini = np.min(data[idx == ii, 0])
            maxi = np.max(data[idx == ii, 0])
            outy.append([mini, maxi])
        outy = sorted(outy[:], key=lambda s: s[1])
        dx = []
        for i in inp:
            dx.append([i[1], 0])
        data = np.vstack(dx)
        centroids, _ = kmeans(data, self.Nx)
        idx, _ = vq(data, centroids)
        for ii in range(0, self.Nx):
            mini = np.min(data[idx == ii, 0])
            maxi = np.max(data[idx == ii, 0])
            outx.append([mini, maxi])
        outx = sorted(outx[:], key=lambda s: s[1])
        out = []
        for y in range(0, len(outy)):
            for x in range(0, len(outx)):
                k = 0
                for k in range(0, len(inp)):
                    if (inp[k][1] >= outx[x][0]) and (inp[k][1] <= outx[x][1]) and (inp[k][2] >= outy[y][0]) \
                            and (inp[k][2] <= outy[y][1]):
                        out.append(inp[k])
                    else:
                        k += 1
        return out

    def order_bb1(self):
        # Order bounding boxes based on grid layout
        f = open(args.grid, 'r')
        for line in f:
            tab = line[:-1].split()
            x, y, w, h = float(tab[1]), float(tab[2]), float(tab[3]), float(tab[4])
            x2 = x+w
            y2 = y+h
            line = NfctnPhntyp(tab[0], self.i_rgb_thre[int(y):int(y2), int(x):int(x2)])
            self.leaf_stack.append(line)
        f.close()

    def m_print(self, inp, level):
        # Print debugging messages based on verbosity level
        if self.verbose > level:
            print(inp)

    def remove_background(self):
        # Remove background from original image
        temp = self.i_original
        return temp

def check_arg(path):
    # Check the range of image numbers in the directory
    import glob
    mstart, mstop = 0, 0
    file_list = glob.glob(path+"/*.png")
    file_list1 = [int(i.split('/')[-1].split('.png')[0]) for i in file_list]
    mstart = min(file_list1)
    mstop = max(file_list1)
    return mstart, mstop

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("mpath", help="Path to the directory containing pictures")
    parser.add_argument("grid", help="Filepath to the grid layout")
    parser.add_argument("-f", "--first", type=int, help="Number of the first picture", default=0)
    parser.add_argument("-l", "--last", type=int, help="Number of the last picture", default=0)
    parser.add_argument("-o", "--output", help="filepath of the output file")
    args = parser.parse_args()

    # Get image range
    start, stop = check_arg(args.mpath)
    if args.first != 0:
        start = args.first
    if args.last != 0:
        stop = args.last

    # Create or open output file
    if not os.path.isfile(args.output):
        f1 = open(args.output, "w")
    else:
        date_now = datetime.date.today().strftime("%B %d, %Y")
        f1 = open(args.output+"_"+date_now+".txt", "w")
    
    # Write headers to output file
    f1.write("ID\ttime\tlesion_area\tleaf_area\n")
    
    # Process each image
    for N in range(start, stop):
        # Zero padding for image number
        N = f'{N:04}'
        out = ""
        sys.stdout.write('\rImage '+args.mpath+"/"+str(N)+".png -> "+str(stop)+".png")
        sys.stdout.flush()
        try:
            # Read image
            image = io.imread(args.mpath+"/"+str(N)+".png")
            # Process panel
            p = Panel(image, args.grid, N)
            for lflt in p.leaf_stack:
                # Quantify infection and write to output file
                lflt.quantify_infection()
                out = lflt.nfst_name+"\t"+str(N)+"\t"+str(lflt.nfctn
