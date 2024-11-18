# genqc

UofM_QC pipeline originally created by Matthew S., curated by SGB Core Platform.

## UofM_QC

Requires Docker & Nextflow to be installed on the computer.
Docker can be installed using Docker Desktop on desktop computers.  
    https://www.docker.com/products/docker-desktop/


Also requires the liftover container to be built on the local computer before being able to use it.
After docker is installed go to the `liftOver-Docker`` directory and run
`docker build . -t liftover`

## Running the QC Pipeline
Run the pipeline by invoking `nextflow` on the QC.nf script here.  You can run it from whichever location you wish as output files will be placed in your current directory.

Supply the following arguments:
`--dataDir='<pathToInputData>'`
The data dir is the path to where your input data can be found.
The data is a .bed, .bim and .fam file.
Make sure this ends with a /
`--dataName='<dataName>'`
The .bed, .bim and .fam file must share the same base name:
ie. data.bed, data.bim, data.fam
In this case, `--dataName=data` would be supplied by this argument.

--liftoverChain='<liftover.chain'>
The liftover chain file to use when performing liftover.
If your input dataset uses hg18 annotations and you want to go forward with hg38, use the hg18ToHg38.over.chain file from http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/

Todo in the future will be to just download the appropriate file for users.

Make sure the file referenced here is in the same location as your --dataDir directory!

## Example:

`nextflow QC.nf --dataDir=/home/matthew/HapMap3/ --dataName=hapmap3_r3_b36_fwd.consensus.qc.poly --liftoverChain=hg18ToHg38.over.chain`

## Running on UofM HPC

module load java/jdk13.0.1
curl -s https://get.nextflow.io | bash
module load singularity


## TODO: Explore Running Nextflow on DNANexus
https://documentation.dnanexus.com/user/running-apps-and-workflows/tools-list

