#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np

FILTER_VALUE = 0.95

# Read data from stdin
df = pd.read_csv(sys.stdin, delimiter='\s+')  # Change delimiter if needed

# Create a DataFrame to store the calculated statistics
result_df = pd.DataFrame(columns=['Column', 'Number of Entries', 'Mean', 'Standard Deviation', 'Sum', 'Minimum', 'Maximum'])

# Iterate through each column in the DataFrame
for column in range(df.shape[1]):  # Use the range of the number of columns
    # Filter values greater than or equal to FILTER_VALUE
    filtered_column = df.iloc[:, column][df.iloc[:, column] >= FILTER_VALUE]
    
    # Calculate the statistics for the filtered column
    num_entries = filtered_column.count()
    mean_value = filtered_column.mean()
    std_dev_value = filtered_column.std()
    sum_value = filtered_column.sum()
    min_value = filtered_column.min()
    max_value = filtered_column.max()

    # Add the results to the result DataFrame
    result_df = pd.concat([result_df, pd.DataFrame({
        'Column': [column + 1],
        'Number of Entries': [num_entries],
        'Mean': [mean_value],
        'Standard Deviation': [std_dev_value],
        'Sum': [sum_value],
        'Minimum': [min_value],
        'Maximum': [max_value]
    })], ignore_index=True)

# Print the result DataFrame or save it to an Excel file, as needed
#print(result_df)
# If you want to save it to an Excel file:
result_df.to_excel('admixture_statistics.xlsx', index=False)
