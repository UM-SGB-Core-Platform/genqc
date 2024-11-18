#!/bin/bash

# Switch Java versions for obtaining the nextflow executable
module load java/jdk13.0.1

# Get the nextflow executable into the current directory
# This command comes from the site https://www.nextflow.io

curl -s https://get.nextflow.io | bash

if [[ -z "${NXF_CACHE_DIR}" ]]; then
	echo Adding NXF_CACHE_DIR environment variable to your environment
	echo Please restart your shell \(log out and back in\)
	echo "export NXF_CACHE_DIR=\"$HOME/projects/def-mliu/NXF_CACHE_DIR\"" >> ~/.bashrc
fi




