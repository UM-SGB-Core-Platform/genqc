import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Use Agg backend for Matplotlib (no display required)
plt.switch_backend('agg')

# Check if there is a column in the directory string
# If so, extract that column
# Returns two values from this function.
def manipulate_dir(dir_string):
    defaultCategory = "Pipeline"
    parts = dir_string.split(os.sep)
    if len(parts) < 3:
        return dir_string, defaultCategory  # Return the original string if there are not enough parts
    match = re.search(r'column(\d+)', parts[1])
    if match:
        return os.sep.join(parts[::2]), "Column" + match.group(1)  # Remove the 2nd part
    else:
        return os.sep.join(parts[1:]), defaultCategory  # Remove the 1st part

# Set the path to the parent directory containing subdirectories
parent_directory = "/home/umstua02/UofM_QC/projectDir/testRun"

# get all directories in the path, and one level of subdirectory
directories = [dirpath for dirpath, dirnames, filenames in os.walk(parent_directory) if dirpath.count(os.sep) - parent_directory.count(os.sep) < 4]
directories = [item for item in directories if '/work/' not in item]
print( directories )
# sort directories by modification time
sorted_directories = sorted(directories, key=os.path.getmtime)


# Define the desired order for the X-axis
desired_order = []

for directory in sorted_directories:
    if( directory.endswith != 'work/' ):
        desired_order.append(directory)
    
# Initialize empty data frame to store data
data = pd.DataFrame()

# Loop through each subdirectory in the desired order
for subdir in desired_order:
    directoryDepth = -3
    processdir = [subdir]
    
    for checkDir in processdir:
        # Create the full path to individuals.txt and variants.txt
        individuals_path = os.path.join(parent_directory, checkDir, "individuals.txt")
        variants_path = os.path.join(parent_directory, checkDir, "variants.txt")

        # Check if both files exist and are non-empty
        if (os.path.exists(individuals_path) and os.path.exists(variants_path)) :
            print( "Processing directory:", checkDir )
            individuals_value = pd.read_csv(individuals_path, header=None).squeeze("columns").astype(float)
            variants_value = pd.read_csv(variants_path, header=None).squeeze("columns").astype(float)

            # Check if both values are non-empty and numeric
            if not individuals_value.empty and not variants_value.empty:
                # Create a data frame with directory name, individuals value, and variants value
                temp_data = pd.DataFrame({
                    'Directory': [os.sep.join(checkDir.split(os.sep)[directoryDepth:])],
                    'Individuals': [individuals_value.mean()],  # Assuming mean for individuals values
                    'Variants': [variants_value.mean()]  # Assuming mean for variants values
                })
                temp_data['Directory'], temp_data['Category'] = zip(*temp_data['Directory'].apply(manipulate_dir))
                #temp_data['Directory'] = temp_data['Directory'].apply(manipulate_dir)

                # Append the temporary data frame to the main data frame
                data = pd.concat([data, temp_data], ignore_index=True)
            else:
                print("Skipping directory:", checkDir, "- Invalid or empty numeric files")
        else:
            print("Skipping directory:", checkDir, "- Missing files")

# Calculate the difference between Initial and Final for Variants
variants_difference = data['Variants'].iloc[-1] - data['Variants'].iloc[0]

# Print the data frame as an excel file
data.to_excel( '/home/umstua02/UofM_QC/pipelineStats.xlsx' )

# Plot for Variants using matplotlib

plt.figure(figsize=(10, 6))
sns.lineplot(x='Directory', y='Variants', hue='Category', data=data, markers=True, palette='colorblind', linewidth=1)
plt.subplots_adjust(bottom=0.4)
plt.title('Variants')
plt.xlabel('Directory')
plt.ylabel('Value')
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.annotate(f'Initial: {data["Variants"].iloc[0]}', xy=(1, data['Variants'].max()), xytext=(1, data['Variants'].max()), ha='left', va='bottom')
plt.annotate(f'Final: {data["Variants"].iloc[-1]}\nDifference: {variants_difference}', xy=(1, data['Variants'].iloc[-1]), xytext=(1, data['Variants'].iloc[-1]), ha='left', va='bottom')

plt.savefig("/home/umstua02/UofM_QC/plot_variants.pdf")
plt.show()
plt.close()


#plt.figure(figsize=(10, 6))
#plt.plot(data['Directory'], data['Variants'], color='blue', marker='o', label='Variants', linewidth=1)
#plt.subplots_adjust(bottom=0.4)
#plt.title('Variants')
#plt.xlabel('Directory')
#plt.ylabel('Value')
#plt.xticks(rotation=45, ha='right')
#plt.legend()
#plt.annotate(f'Initial: {data["Variants"].iloc[0]}', xy=(1, data['Variants'].max()), xytext=(1, data['Variants'].max()), ha='left', va='bottom')
#plt.annotate(f'Final: {data["Variants"].iloc[-1]}\nDifference: {variants_difference}', xy=(1, data['Variants'].iloc[-1]), xytext=(1, data['Variants'].iloc[-1]), ha='left', va='bottom')

#plt.savefig("/home/umstua02/UofM_QC/plot_variants.pdf")
#plt.show()

# Calculate the difference between Initial and Final for Individuals
individuals_difference = data['Individuals'].iloc[-1] - data['Individuals'].iloc[0]

# Plot for Individuals using matplotlib
plt.figure(figsize=(10, 6))
sns.lineplot(x='Directory', y='Individuals', hue='Category', data=data, markers=True, palette='colorblind', linewidth=1)
plt.subplots_adjust(bottom=0.4)
#plt.plot(data['Directory'], data['Individuals'], color='red', marker='o', label='Individuals', linewidth=1)
plt.title('Individuals')
plt.xlabel('Directory')
plt.ylabel('Value')
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.annotate(f'Initial: {data["Individuals"].iloc[0]}', xy=(1, data['Individuals'].max()), xytext=(1, data['Individuals'].max()), ha='left', va='bottom')
plt.annotate(f'Final: {data["Individuals"].iloc[-1]}\nDifference: {individuals_difference}', xy=(1, data['Individuals'].iloc[-1]), xytext=(1, data['Individuals'].iloc[-1]), ha='left', va='bottom')
plt.savefig("/home/umstua02/UofM_QC/plot_individuals.pdf")
plt.show()

# Save the Variants plot as a PDF file
#plt.savefig("/home/umstua02/UofM_QC/plot_variants.pdf")
# Close the plot to free up resources
#plt.close()

# Save the Individuals plot as a PDF file
#plt.savefig("/home/umstua02/UofM_QC/plot_individuals.pdf")
# Close the plot to free up resources
#plt.close()
