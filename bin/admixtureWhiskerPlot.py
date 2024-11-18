#!/bin/env python

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np

MEMBERSHIP_CUTOFF = 0.85

# Check if a filename was provided
if len(sys.argv) < 2:
    print("Please provide a file name.")
    sys.exit()
    
boxplotName = "boxplots.pdf"
if len(sys.argv) < 3:
    print("No output filename provided.  Saving as \"boxplots.pdf\"")
else:
    boxplotName = sys.argv[2]

# Load the data
data = pd.read_csv(sys.argv[1], sep=' ')

# Filter the data
data = data[data > MEMBERSHIP_CUTOFF]

# Create a box plot for each column
fig, ax1 = plt.subplots(figsize=(10,5))
data.boxplot(ax=ax1)

# Set the label for the left y-axis
ax1.set_ylabel('Probability of group membership')

# Get the number of data points in each column
num_datapoints = data.count()

# Create a secondary y-axis for the number of data points
ax2 = ax1.twinx()
ax2.plot(ax1.get_xticks(), num_datapoints, color='red')
ax2.set_ylabel('Number of data points')

# Print the number of data points above each box plot
for i, column in enumerate(data.columns):
    ax1.text(i + 1, ax1.get_ylim()[1], 'N = {}'.format(num_datapoints[i]), ha='center', va='bottom')

# Save the plot as a PDF
plt.tight_layout()
plt.savefig(boxplotName)

plt.close('all')
