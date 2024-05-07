import sys
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import numpy as np

## Check Results - Panel from Fiji as overlay over ímage
## Severin Einspanier
## 07.05.2024


# filepath of grid layout
grid = sys.argv[1]
# directory of images
images = sys.argv[2]
# which image to align to grid
img = sys.argv[3]
# output directory for the panel image
outdir = sys.argv[4]

back = mpimg.imread(images + "/" + img + ".png")
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6, 6))
ax.imshow(back)

# copy/paste from FIJI measurements (bounding rectangle) into .csv
grid_f = open(grid, "r")
data = []
for line in grid_f.readlines():
   data.append(line.replace(",", "\t").split("\t"))

# Prompt the user to check if the first row has been removed
while True:
   first_row_removed = input("Has the first row of data already been removed? (y/n): ")
   if first_row_removed.lower() in ["y", "n"]:
       break
   else:
       print("Invalid input. Please enter 'y' or 'n'.")

# Remove the first row of data if necessary
if first_row_removed.lower() == "n":
   data = data[1:]

   
# Extract coordinates and text
for row in data:
    x, y, w, h = int(row[1]), int(row[2]), int(row[3]), int(row[4])
    text = row[0]

    # Create and add rectangle
    rect = mpatches.Rectangle((x, y), w, h, fill=False, edgecolor='red', linewidth=0.25)
    ax.add_patch(rect)

    # Add text centered within the rectangle
    ax.text(x + w / 2, y + h / 2, text, fontsize=2, color="black")

grid_f.close()

# Save the figure to the specified output directory
fig.savefig(outdir + "/" + img + "_panel.png", dpi=300)
print("FERTIG :)")

# python3 check_image_grid.py Results.csv raw/ 0876 .
