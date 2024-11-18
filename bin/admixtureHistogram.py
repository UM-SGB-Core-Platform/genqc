#!/bin/env python

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import sys

# Check if a filename was provided
if len(sys.argv) < 2:
    print("Please provide a file name.")
    sys.exit()

# Load the data
data = pd.read_csv(sys.argv[1], sep=' ')

# Determine the number of rows and columns for the subplot grid
num_columns = len(data.columns)
num_rows = num_columns // 2
if num_columns % 2:
    num_rows += 1

# Create a figure and axes for the subplot grid
fig, axes = plt.subplots(num_rows, 2, figsize=(10, num_rows*5))

# Flatten the axes array
axes = axes.flatten()

# Loop through each column
for i, column in enumerate(data.columns):
    # Create a histogram for the column on the corresponding subplot
    data[column].hist(bins=50, ax=axes[i])
    axes[i].set_title('Histogram of column {}'.format(i))
    axes[i].set_xlabel(column)
    axes[i].set_ylabel('Frequency')

# Remove any unused subplots
for j in range(num_columns, len(axes)):
    fig.delaxes(axes[j])

# Save the plot as a PDF
plt.tight_layout()
plt.savefig('histograms.pdf')

plt.close('all')
