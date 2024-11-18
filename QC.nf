#!/usr/bin/env nextflow

params.dataDir = '/mnt/c/Users/Matthew/Downloads/HapMap3/' // Make sure you put a / at the end
params.dataName = 'hapmap3_r3_b36_fwd.consensus.qc.poly'

// Liftover parameters.  Should be in the same directory as the dataDir.  Might change that.
// Can be obtained from http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/

// TODO: Need a better way to find the liftover files
params.liftoverChainLocation = '/home/umstua02/UofM_QC/liftoverChainFiles/'

params.hg18File = 'hg18ToHg38.over.chain'
params.hg19File = 'hg19ToHg38.over.chain'
params.hg17File = 'hg17ToHg38.over.chain'
params.hg16File = 'hg16ToHg38.over.chain'

params.convertToHg38 = true

// TODO: Determine what annotation set is in use and determine which liftover to use.
params.convertFrom = params.hg19File

//plinkContainer = /*Singularity*/ 'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' //Docker: 'biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
plinkContainer = /*Singularity*/ 'plink1.9_v1.90b6.6-181012-1-deb_cv1.sif' //Docker: 'biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
    

// LiftOver requires a license.  It is free for non-commercial use, which is the intended purpose of this pipeline,
// but it means that distributing a public liftover container is not acceptable.  
// I have created a container for liftover with the executable.  Please have it available on the system you are using.
liftoverContainer = 'liftover.sif'

//admixtureContainer = 'admixture' 
params.admixtureCommand = 'admixture' //plain admixture crashes


// Sometimes plink is called plink.  Sometimes it's plink1.9.  Lets create a variable for it.
params.plinkCommand = 'plink1.9'

plinkOutputRegex = '\'(\\d+) variants and (\\d+) people pass filters and QC\''

getVariantsAndIndividualsFunction = """
    latest_log_file=\$(ls -t *.log | head -n 1)
    plinkResults=\$(grep -P $plinkOutputRegex \$latest_log_file)

    echo "\$plinkResults" | grep -oE '[0-9]+' | head -n 1 > variants.txt
    echo "\$plinkResults" | grep -oE '[0-9]+' | tail -n 1 > individuals.txt
"""

process PLINK_HG38_UPDATE {
    container plinkContainer
    //container liftoverContainer

    input:
    path pedFile
    path mapFile
    path oldBimFile
    path oldBedFile
    path oldFamFile
    path chainFile
    val workflowName

    output:
    path '*.bed', emit: bedFile
    path '*.bim', emit: bimFile
    path '*.fam', emit: famFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    // TODO publish log files
    //path '*.log', emit: logFiles

    script:
    // Nextflow and Plink do not work well together, because the file you provide to Plink is often a basename, not a full file.
    // Nextflow wants a full file to stage in to the workign directory.
    // Therefore, to merge the two, you need to tell Nextflow the exact files you want, while providing Plink the basename of those files.
    // So even if an input isn't referened, it is required so that Nextflow will stage that file.

    baseFile = pedFile.baseName
    // PLINK 1 binary is PLINK 1.9's preferred input format. In fact, PLINK 1.9 automatically converts most other formats to PLINK 1 binary before the main loading sequence
    // So we'll look for the binary files first.
    
    bedFile = baseFile + ".updated"

    """
    if [ -e ${oldBimFile} ] && [ -e ${oldBedFile} ] && [ -e ${oldFamFile} ]; then
        ${params.plinkCommand} --bfile ${baseFile} --make-bed --out $bedFile
    elif [ -e ${pedFile} ] && [ -e ${mapFile} ]; then
        ${params.plinkCommand} --file ${baseFile} --make-bed --out $bedFile
    fi

    $getVariantsAndIndividualsFunction
    """
}

process LIFTOVER_PLINK_HG38 {
    container liftoverContainer
    memory '24 GB'

    input:
    path mapFile
    path pedFile
    path chainFile
    val workflowName

    output:
    path "lifted.ped", emit: pedFile
    path "lifted.map", emit: mapFile
    path "lifted.bed.unlifted", emit: unlifted, optional: true

    publishDir "${workflowName}/LiftOver", mode: 'symlink'

    script:
    """
    #Turn chr 23 and 24 into X and Y.  If no changes are needed then this will just print the original map file.
    awk '{if (\$1 == "23") \$1 = "X"; else if (\$1 == "24") \$1 = "Y"; print}' $mapFile | sed 's/ /\\t/g' > chrUpdated.map
    ${workflow.projectDir}/bin/liftOverPlink.py -m chrUpdated.map -p $pedFile -o lifted -c $chainFile
    """
}

// TODO: This was taken out because the output wasn't properly being processed.
// Have opted for LIFTOVER_PLINK_HG38 instead.
process LIFTOVER_HG38 {
    container liftoverContainer
    memory '24 GB'

    input:
    path bedFile
    path bimFile
    path famFile
    path chainFile
    val workflowName

    output:
    path "*_liftover", emit: lifted
    path "*_unlifted", emit: unlifted
    path bedFile, emit: bedFile
    path bimFile, emit: bimFile
    path famFile, emit: famFile
    path "*.pos", emit: posFile
    path "*.snps", emit: snpsFile


    publishDir "${workflowName}/LiftOver", mode: 'symlink'

    script:
    baseFile = bimFile.baseName
    newFileLifted = baseFile + "_liftover"
    newFileUnlifted = baseFile + "_unlifted"

    """


    awk '{print "chr" \$1, \$4 -1, \$4, \$2 }' ${bimFile} | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${bedFile}.toLift
    liftOver ${bedFile}.toLift $chainFile $newFileLifted $newFileUnlifted
    awk '{print \$4}' $newFileLifted > ${newFileLifted}.snps
    awk '{print \$4, \$3}' $newFileLifted > ${newFileLifted}.pos
        """
}
// TODO: This isn't needed with the converion of liftover program to LIFTOVER_PLINK_HG38
process POST_LIFTOVER_PLINK{
    container plinkContainer
    time '24h'
    memory '96 GB'

    input:
    path bedFile
    path bimFile
    path famFile
    path snpsFile
    path posFile
    val workflowName

    output:
    path "hg38.bed", emit: bedFile
    path 'hg38.bim', emit: bimFile
    path 'hg38.fam', emit: famFile
    path 'hg38.hh', emit: hhFile, optional: true
    path 'hg38.log', emit: logFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed

    publishDir "${workflowName}/UpdatedReferences", mode: 'symlink'

    script:
    baseFile = bedFile.baseName
    """
    ${params.plinkCommand} --bfile $baseFile --extract $snpsFile --update-map $posFile --make-bed --out "hg38"

    $getVariantsAndIndividualsFunction
    """

}

process MARKER_MISSINGNESS{
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path "*.bed", emit: bedFile
    path '*.bim', emit: bimFile
    path '*.fam', emit: famFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column

    publishDir "${workflowName}/$column/MarkerMissingness", mode: 'symlink'

    script:
    def baseFile = plinkFam.baseName

    // First make sure that the ped file has 
    // the same basename as the bedfile. It might
    // have been renamed by liftOver.
    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --geno 0.05 --out "${baseFile}.geno" --make-bed

    $getVariantsAndIndividualsFunction
    """
}

process CONVERT_TO_PLINK_FORMAT{
    container plinkContainer

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName

    output:
    path '*.map', emit: mapFile
    path '*.ped', emit: pedFile

    script:
    baseFile = plinkFam.baseName
    """
    ${params.plinkCommand} --bfile $baseFile --recode --tab --out pedOut
    """
}

process CONVERT_TO_PLINK_BINARY_FORMAT{
    container plinkContainer

    input:
    path plinkMap
    path plinkPed
    val workflowName

    output:
    path '*.bim', emit: bimFile
    path '*.bed', emit: bedFile
    path '*.fam', emit: famFile

    script:
    baseFile = plinkMap.baseName
    """
    ${params.plinkCommand} --file $baseFile --allow-extra-chr --make-bed --out bedOut
    """
}

process PLINK_FILTER_BY_IDS{
    container plinkContainer

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    path filterList
    val workflowName

    output:
    path "columnsFiltered.bed", emit: bedFile
    path 'columnsFiltered.bim', emit: bimFile
    path 'columnsFiltered.fam', emit: famFile
    path 'columnsFiltered.log', emit: logFile
    path 'columnsFiltered.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val columnStripped, emit: column
    //stdout, emit: log

    publishDir "${workflowName}/${columnStripped}/PlinkFilterIDS", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName
    columnStripped = filterList.simpleName.replaceAll(/_.*/, '')
    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --keep $filterList --make-bed --out columnsFiltered

    $getVariantsAndIndividualsFunction
    """
}


// TODO: Revisit this command, might not be necessary.
process MENDELIAN_ERROR_RATE {
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName

    output:
    path "mendelianErrorFiltered.bed", emit: bedFile
    path 'mendelianErrorFiltered.bim', emit: bimFile
    path 'mendelianErrorFiltered.fam', emit: famFile
    path 'mendelianErrorFiltered.log', emit: logFile
    path 'mendelianErrorFiltered.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    //stdout, emit: log

    publishDir "${workflowName}/MendelianErrorRate", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    // After running the --mendel command, you will see the Mendelian error rate in the PLINK output.
    // You can then decide on a threshold for filtering markers based on the Mendelian error rate.

    // Use the --geno option to filter out markers with a Mendelian error rate greater than 1%. 
    // To do this, you will need to calculate the threshold based on the error rate. 
    // If the error rate is in the range of 0-1, you can calculate the threshold by taking 1% 
    // of the error rate and then use the --geno option. 
    // For example, if the Mendelian error rate is 0.02 (2%), you can use a --geno threshold of 0.02 to 
    // filter out markers with an error rate greater than 1%:

    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --mendel > mendel.txt
    error_rate=\$(grep 'Mendel error rates' mendel.txt | awk '{print \$NF}' | sed 's/,//')
    echo "Mendelian error rate: \$error_rate"
    
    threshold=0.01
    if (( \$(echo "\$error_rate \$threshold" | awk '{if (\$1 >= \$2) print 1;}') )); then
        echo "Mendelian error rate is greater than \$threshold (1%)."

        # Filter markers with error rate greater than 1%
        ${params.plinkCommand} --bfile $baseFile  --allow-extra-chr --geno \$threshold --make-bed --out mendelianErrorFiltered

        echo "Markers with high Mendelian error rate filtered and saved as mendelianErrorFiltered."
    else
        echo "Mendelian error rate is within acceptable limits."
    fi

    $getVariantsAndIndividualsFunction
    """
}

 process INDIVIDUALS_AUTOSOMAL_HETEROZYGOSITY {
    container plinkContainer

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path "plink.het", emit: hetCheck
    path "hetFiltered.bim", emit: bimFile
    path "hetFiltered.bed", emit: bedFile
    path "hetFiltered.fam", emit: famFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column

    publishDir "${workflowName}/$column/CheckHeterozygosity", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    ${params.plinkCommand} --bfile $baseFile  --allow-extra-chr --het

    # Calculate the mean and stddev of the het file's F score
    # Remove individuals who are 3 stdevs away from the mean.
    #!/bin/bash

    # File name
    file="plink.het"

    # Calculate mean and standard deviation
    mean=\$(awk '{ total += \$6 } END { print total/NR }' \$file)
    stddev=\$(awk -v mean=\$mean '{ total += (\$6 - mean)^2 } END { print sqrt(total/NR) }' \$file)

    # Print 1st and 2nd column where the 6th column is outside 3 standard deviations from the mean
    awk -v mean=\$mean -v stddev=\$stddev '(\$6 < mean - 3*stddev) || (\$6 > mean + 3*stddev) { print \$1, \$2 }' \$file > hetFiltered.txt

    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --remove hetFiltered.txt --make-bed --out hetFiltered
    
    $getVariantsAndIndividualsFunction
    """
}

process CHECK_SEX {
    container plinkContainer

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path "plink.sexcheck", emit: sexCheck
    path "sexFiltered.bim", emit: bimFile
    path "sexFiltered.bed", emit: bedFile
    path "sexFiltered.fam", emit: famFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column

    publishDir "${workflowName}/$column/CheckSex", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --check-sex

    # remove PROBLEM individuals
    awk '\$5 != "OK" {print \$1 " " \$2}' plink.sexcheck > toFilter.txt
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --remove toFilter.txt --make-bed --out sexFiltered

    $getVariantsAndIndividualsFunction
    """
}

process MAF_AND_MISSING {
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path "maf_missing.bed", emit: bedFile
    path 'maf_missing.bim', emit: bimFile
    path 'maf_missing.fam', emit: famFile
    path 'maf_missing.log', emit: logFile
    path 'maf_missing.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column
    //stdout, emit: log

    publishDir "${workflowName}/$column/MafMissing", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    # max-maf obtains the maximimum maf.  So this filters out all variants above that. Nice feature for this task.
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --max-maf 0.05 --make-bed --out maf_inverse
    # geno filters variants.  This is more lax than previous --geno runs so I called it 'lax'
    # Need to get a list of the Id's NOT removed from this filter.
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --geno 0.01 --make-bed --out missing_lax
    awk '{print \$2}' missing_lax.bim > missing_lax_ids.txt

    # Now go back to the original base file and exclude all the variants that weren't removed before.
    # Now we have all variants that should have failed the 0.01 geno filter.
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --exclude missing_lax_ids.txt --make-bed --out missing_lax_inverse
    awk '{print \$2}' missing_lax_inverse.bim > missing_lax_inverse.txt

    # Load the max-maf output file and remove the variants that failed the 0.01 geno filter.
    # Now we have a list of variants that both are below 0.05 from max-maf and belwo 0.01 from geno.
    ${params.plinkCommand} --bfile maf_inverse --allow-extra-chr --extract missing_lax_inverse.txt --make-bed --out missing_and_maf
    awk '{print \$2}' missing_and_maf.bim > missing_and_maf.txt

    # Go back to the original input and exclude the variants identified as meeting both criteria for maf and geno
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --exclude missing_and_maf.txt --make-bed --out maf_missing

    $getVariantsAndIndividualsFunction
    """
}

process MINOR_ALLELE_FREQUENCY{
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path "maf.bed", emit: bedFile
    path 'maf.bim', emit: bimFile
    path 'maf.fam', emit: famFile
    path 'maf.log', emit: logFile
    path 'maf.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column
    //stdout, emit: log

    publishDir "${workflowName}/$column/MinorAlleleFrequency", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr -maf 0.01 --make-bed --out maf

    $getVariantsAndIndividualsFunction
    """
}

process HARDY_WEINBERG_EQUILIBRIUM{
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path "hwe.bed", emit: bedFile
    path 'hwe.bim', emit: bimFile
    path 'hwe.fam', emit: famFile
    path 'hwe.log', emit: logFile
    path 'hwe.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column
    //stdout, emit: log

    publishDir "${workflowName}/$column/HardyWeinbergEquilibrium", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --hwe 0.0001 --make-bed --out hwe
    
    $getVariantsAndIndividualsFunction
    """
}

process FILTER_UNAFFTECTED_AND_NO_PARENTS{
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName

    output:
    path "familyAffectedFilter.bed", emit: bedFile
    path 'familyAffectedFilter.bim', emit: bimFile
    path 'familyAffectedFilter.fam', emit: famFile
    path 'familyAffectedFilter.log', emit: logFile
    path 'familyAffectedFilter.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    //stdout, emit: log

    publishDir "${workflowName}/FamilyFiltered", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    // Gotta check that this actually is working correctly.
    """
    ${params.plinkCommand} --bfile $baseFile --family --list-duplicate-vars ids-only suppress-first --out families

    awk '{print \$1}' $plinkFam | sort -u > all_families.txt

    awk 'NR==FNR{a[\$1];next} !(\$1 in a){print \$1}' families.dupvar all_families.txt > families_to_keep.txt

    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --family --keep-fam families_to_keep.txt --make-bed --out familyAffectedFilter

    $getVariantsAndIndividualsFunction
    """
}

process SAME_LOCATION_FILTER{
    container plinkContainer
    // This really only needs the bim file.  Considering
    // adding the others just to help with the workflow, but
    // lets leave it at the minimal needed for now.
    //path plinkBed
    input:
    path plinkBed
    path plinkBim
    path plinkFam
    val workflowName

    output:
    path 'updated_dataset.bed', emit: bedFile
    path 'updated_dataset.bim', emit: bimFile
    path 'updated_dataset.fam', emit: famFile
    path 'dupID.txt', emit: dupFile
    path 'duplicates.txt', emit: duplicateEvidence
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed

    publishDir "${workflowName}/SameLocations", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    # Escape all \$var references for bash script.
    #
    # Look at the .bim file, chr and position will be the same in some locations.
    # chr is 1st column
    # pos is the 4th column
    # This will output duplicate ID's which can then be filtered.

    # I need to feed these locations to plink to filter them.  Filtering the .bim file 
    # itself will result in a mismatched .bed file.

    # If no duplicates found, have an empty file on hand
    touch "duplicates.txt"

    awk '{
        key = \$1 " " \$4;
        if (seen[key] == 1) {
            print initialDup[key], \$2 > "duplicates.txt";
            print \$0, \$2 > "duplicates.txt";
            duplicates[key]++;
        } else if (seen[key] > 1) {
            print \$0, \$2 > "duplicates.txt";
            duplicates[key]++;
        } else {
            seen[key] = 1;
            initialDup[key] = \$0;
        }
    }' "$plinkBim"

    # TODO: Remove both duplicates, not just the first duplicate instance found.
    awk '{print \$2}' duplicates.txt > dupID.txt

    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --exclude dupID.txt --make-bed --out updated_dataset

    $getVariantsAndIndividualsFunction
    """
}

process AT_GC_FILTER{
    container plinkContainer

    input:
    path plinkBed
    path plinkBim
    path plinkFam
    val workflowName

    output:
    path 'indelsRemoved.bed', emit: bedFile
    path 'indelsRemoved.bim', emit: bimFile
    path 'indelsRemoved.fam', emit: famFile
    path 'AT_GC_locations.txt', emit: indelMarkers
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed

    publishDir "${workflowName}/AT_GC_Filter", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    # Find an I in the 5th or 6th column
    # print out those locations and use plink to filter them.
    awk '\$5 == "A" && \$6 == "T" || \$5 == "C" && \$6 == "G" || \$5 == "T" && \$6 == "A" || \$5 == "G" && \$6 == "C" {print \$2 " " \$5 " " \$6}' $plinkBim > AT_GC_locations.txt
    awk '{print \$1}' AT_GC_locations.txt > filterLocations.txt
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --exclude filterLocations.txt --make-bed --out indelsRemoved

    $getVariantsAndIndividualsFunction
    """
}

process INDEL_MARKER_FILTER{
    container plinkContainer

    input:
    path plinkBed
    path plinkBim
    path plinkFam
    val workflowName

    output:
    path 'indelsRemoved.bed', emit: bedFile
    path 'indelsRemoved.bim', emit: bimFile
    path 'indelsRemoved.fam', emit: famFile
    path 'indelLocations.txt', emit: indelFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed

    publishDir "${workflowName}/IndelLocations", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    # Find an I in the 5th or 6th column
    # print out those locations and use plink to filter them.

    awk '(\$5 == "I" || \$6 == "I") {print \$2}' $plinkBim > indelLocations.txt

    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --exclude indelLocations.txt --make-bed --out indelsRemoved

    $getVariantsAndIndividualsFunction
    """
}

process INDIVIDUALS_CNV_FAMILIES_OR_IDIOPATHIC_FILTER{
    container plinkContainer
    // This really only needs the bim file.  Considering
    // adding the others just to help with the workflow, but
    // lets leave it at the minimal needed for now.
    //path plinkBed
    input:
    path plinkBed
    path plinkBim
    path plinkFam
    val workflowName

    output:
    path 'updated_dataset.bed', emit: bedFile
    path 'updated_dataset.bim', emit: bimFile
    path 'updated_dataset.fam', emit: famFile
    path 'dup.txt', emit: dupFile
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed

    publishDir "${workflowName}/SameLocations", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    """
    # Escape all \$var references for bash script.
    #
    # Look at the .bim file, chr and position will be the same in some locations.
    # chr is 1st column
    # pos is the 4th column
    # This will output duplicate ID's which can then be filtered.
    # This will retain the original location though. 
    # I feel that duplicate locations mean there is likely a good SNP in that location,
    # just that with the updated genome assembly positions sometimes get merged.

    # I need to feed these locations to plink to filter them.  Filtering the .bim file 
    # itself will result in a mismatched .bed file.

    # print duplicates to dup.txt and print everything else to stdout
    awk '{
        key = \$1 " " \$4;
        if (seen[key] == 1) {
            print \$1, \$2, \$4 > "dup.txt";
            duplicate[key] = 1;
        } else {
            seen[key] = 1;
        }
    }' < $plinkBim
    awk '{print \$2}' dup.txt > dupID.txt

    #${baseDir}/markerFilterScripts/sameLocationsFilter.sh < $plinkBim > 'locationsFiltered.bim'
    #ln -s $plinkBed locationsFiltered.bed
    #ln -s $plinkFam locationsFiltered.fam
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --exclude dupID.txt --make-bed --out updated_dataset

    $getVariantsAndIndividualsFunction
    """
}

process INDIVIDUALS_MISSINGNESS{
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path 'individualMissing.bed', emit: bedFile
    path 'individualMissing.bim', emit: bimFile
    path 'individualMissing.fam', emit: famFile
    path 'individualMissing.log', emit: logFile
    path 'individualMissing.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column
    //stdout, emit: log

    publishDir "${workflowName}/$column/IndividualMissingness", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    // Gotta check that this actually is working correctly.
    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --mind 0.05 --make-bed --out individualMissing

    $getVariantsAndIndividualsFunction
    """
}

process INDIVIDUALS_MENDELIAN_ERROR{
    container plinkContainer


    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName
    val column

    output:
    path 'individualMendelianError.bed', emit: bedFile
    path 'individualMendelianError.bim', emit: bimFile
    path 'individualMendelianError.fam', emit: famFile
    path 'individualMendelianError.log', emit: logFile
    path 'individualMendelianError.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    val column, emit: column
    //stdout, emit: log

    publishDir "${workflowName}/$column/IndividualMendelianError", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    // Gotta check that this actually is working correctly.
    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --me 0.01 0.01 --make-bed --out individualMendelianError

    $getVariantsAndIndividualsFunction
    """
}

process ADMIXTURE_PRUNE{
    container plinkContainer

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val workflowName

    output:
    path 'ldPrune.bed', emit: bedFile
    path 'ldPrune.bim', emit: bimFile
    path 'ldPrune.fam', emit: famFile
    path 'ldPrune.log', emit: logFile
    path 'ldPrune.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    //stdout, emit: log

    publishDir "${workflowName}/AdmixturePrune", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName

    // Gotta check that this actually is working correctly.
    // Remove sex chromosomes
    // Filter linkage disequilibrium markers
    // filter founders.  Only non-founders.
    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --chr 1-22 --indep 50 5 1.25 --filter-founders --make-bed --out pruned
    ${params.plinkCommand} --bfile pruned --allow-extra-chr --extract pruned.prune.in --make-bed --out ldPrune 

    $getVariantsAndIndividualsFunction
    """
}

process ADMIXTURE_JOIN_REFERENCE{
    container plinkContainer
    errorStrategy 'retry'
    maxRetries 1

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    path referenceBedFile
    path referenceBimFile
    path referenceFamFile
    val workflowName

    output:
    path 'referenceJoined.bed', emit: bedFile
    path 'referenceJoined.bim', emit: bimFile
    path 'referenceJoined.fam', emit: famFile
    path 'referenceJoined.log', emit: logFile
    path 'referenceJoined.hh', emit: hhFile, optional: true
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed
    //stdout, emit: log

    publishDir "${workflowName}/AdmixtureReferencesJoined", mode: 'symlink'

    script:
    baseFile = plinkFam.baseName
    referenceBase = referenceFamFile.baseName
    doFlip = task.attempt // If this is a retry, then do the SNP flip.

    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --write-snplist 
    mv plink.snplist markers_base.txt
    ${params.plinkCommand} --bfile $referenceBase --allow-extra-chr --write-snplist
    mv plink.snplist markers_reference.txt
    awk 'NR==FNR{a[\$1];next} (\$1 in a)' markers_base.txt markers_reference.txt > common_markers.txt

    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --extract common_markers.txt --make-bed --out base_filtered
    ${params.plinkCommand} --bfile $referenceBase --allow-extra-chr --extract common_markers.txt --make-bed --out reference_filtered


    if [ $doFlip ]; then
        ${params.plinkCommand} --bfile base_filtered --allow-extra-chr --bmerge reference_filtered --make-bed --out referenceNeedFlip || echo "Flipping heck!"
        ${params.plinkCommand} --bfile reference_filtered --allow-extra-chr --flip "referenceNeedFlip-merge.missnp" --make-bed --out referenceFlipped
        ${params.plinkCommand} --bfile base_filtered --allow-extra-chr --bmerge referenceFlipped --make-bed --out referenceJoined
    else
        ${params.plinkCommand} --bfile base_filtered --allow-extra-chr --bmerge reference_filtered --make-bed --out referenceJoined
    fi

    $getVariantsAndIndividualsFunction
    """
}

process FILTER_KING {
    container plinkContainer
    //cpus 12
    //memory '32 GB'
    //time '3d'
    //clusterOptions '--ntasks=1'
    //queue = 'compute'
    

    input:
    //path plinkBed
    //path plinkBim // Needed even if not used in the script
    //path plinkFam // Needed even if not used in the script
    tuple path(plinkBed), path(plinkBim), path(plinkFam),path(kinFile), val (column)
    //path kinFile 
    val workflowName
    //val column

    output:
    path 'kinsFiltered.bed', emit: bedFile // Just outputting the inputs since nothing will be changed.
    path 'kinsFiltered.bim', emit: bimFile
    path 'kinsFiltered.fam', emit: famFile
    path 'kinsMissing.txt', emit: kinsMissingness
    path 'higherMissingKins.txt', emit: higherMissingKings
    val column, emit: column
    path 'variants.txt', emit: variantsPassed
    path 'individuals.txt', emit: individualsPassed

    publishDir "${workflowName}/$column/FilteredPostKing", mode: 'copy'

    script:
    baseFile = plinkFam.baseName

    """
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --missing

    # I think king.kin0 is the only output file created by King.
    # Put both individuals into a new file.  We'll know they're kins due to one following the other.
    tail -n +2 king.kin0 | awk '{print \$1, \$2 > "individualKins.txt"; print \$3, \$4 >> "individualKins.txt"}'

    # Match the first two columns in the indivualKins.txt file with the plink results, put that into a file so we have just the kins and their missing value.
    awk 'NR==FNR{a[\$1,\$2]=\$0;next} {if((\$1,\$2) in a) print \$1 " " \$2 " " \$6}' individualKins.txt plink.imiss > kinsMissing.txt

    # Now for each two lines (representing a pair of kins) print only the kin with the higher missingness value
    awk 'NR%2==1{line=\$0; max=\$3; next} {if(\$3>max){print \$0} else {print line}}' kinsMissing.txt > higherMissingKins.txt

    # Use Plink to filter (remove) the kin with the higher missing value, thereby leaving the kin with the lower missing value.
    ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --remove higherMissingKins.txt --make-bed --out kinsFiltered

    $getVariantsAndIndividualsFunction
    """
}

process KING {
    //cpus 12
    //memory '32 GB'
    //time '3d'
    //clusterOptions '--ntasks=1'
    //queue = 'compute'
    

    input:
    path plinkBed
    path plinkBim
    path plinkFam
    val workflowName
    val column

    output:
    //path plinkBed, emit: bedFile
    //path plinkBim, emit: bimFile
    //path plinkFam, emit: famFile
    tuple (path plinkBed), (path plinkBim), (path plinkFam), (path '*.kin*'), val (column), emit: kingOut, optional: true

    publishDir "${workflowName}/$column/King", mode: 'copy'

    script:
    baseFile = plinkFam.baseName

    """
    ${workflow.projectDir}/bin/king -b $plinkBed --related
    """
}


process ADMIXTURE{
    // Don't need this on Grex, they have admixture installed.
    //container admixtureContainer
    cpus 12
    memory '32 GB'
    time '3d'
    //clusterOptions '--ntasks=1'
    queue = 'compute'

    input:
    path plinkBed
    path plinkBim // Needed even if not used in the script
    path plinkFam // Needed even if not used in the script
    val k
    val workflowName

    output:
    //path 'individualMendelianError.bed', emit: bedFile
    //path 'individualMendelianError.bim', emit: bimFile
    //path 'individualMendelianError.fam', emit: famFile
    //path 'individualMendelianError.log', emit: logFile
    //path 'individualMendelianError.hh', emit: hhFile, optional: true
    path '*.out', emit: logFiles
    path '*.P', emit: pFile
    path '*.Q', emit: qFile    

    publishDir "${workflowName}/AdmixtureOutput", mode: 'copy'

    script:
    baseFile = plinkFam.baseName

    // Gotta check that this actually is working correctly.
    // 
    """
    ${params.admixtureCommand} --cv -j12 $plinkBed ${k} | tee log${k}.out
    """
}

process ADMIXTRUE_EXTRACT_FROM_GROUP{
    // Extracting the ids from the input dataset is difficult.
    // When we joined the two groups together to do Admixture, we have a dataset with both.
    // We want to go through the admixture output and pull ids from the specific categories (columns), but excude the reference dataset.

    input:
    path qFile
    path joinedFamFile
    path referenceFamFile
    path famFile
    val workflowName

    output:
    path "admixture_columns.txt", emit: admixtureColumnsOutput
    path "admixture_reference_columns.txt", emit: admistureReferenceColumnsOutput
    path "AdmixtureGroupMembershipBoxplot.pdf", emit: boxplot
    path "AdmixtureReferenceGroupMembershipBoxplot.pdf", emit: referenceBoxplot
    path "column*_ids_filtered.txt", emit: filteredIdsByColumn

    publishDir "${workflowName}/AdmixtureGroupMembership"

    script:

    """
    #  Pull the list of IDs that we want to extract.  These are not the reference IDs.
    awk '{print \$2}' $famFile > ids.txt
    # line_numbers.txt will contain line numbers from the joined family file that contain only the input dataset IDs.
    # We can then use this to exclude reference dataset and pull out our input dataset.
    awk 'NR==FNR{ids[\$0]; next} \$2 in ids {print FNR}' ids.txt $joinedFamFile > line_numbers.txt
    # This will get the Ids of the input dataset
    awk 'NR==FNR{ids[\$0]; next} \$2 in ids {print \$1 " " \$2}' ids.txt $joinedFamFile > line_numbers_ids.txt
    # This will remove lines from the admixture Q file that come from the reference dataset, leaving only the input data.
    awk 'NR==FNR{lines[\$0]; next} (FNR in lines)' line_numbers.txt $qFile > admixture_columns.txt

    # Repeat for the reference data.  Need to have this information to know if we have admixture representation 
    # This means we know enough about the population in a column to infer the ancenstry.
    # Note that line_numbers.txt is intended for short term use.
    awk '{print \$2}' $referenceFamFile > ids.txt
    awk 'NR==FNR{ids[\$0]; next} \$2 in ids {print FNR}' ids.txt $joinedFamFile > line_numbers.txt
    awk 'NR==FNR{ids[\$0]; next} \$2 in ids {print \$1 " " \$2}' ids.txt $joinedFamFile > line_numbers_reference_ids.txt
    awk 'NR==FNR{lines[\$0]; next} (FNR in lines)' line_numbers.txt $qFile > admixture_reference_columns.txt
    rm line_numbers.txt

    # Generate graphs
    ${workflow.projectDir}/bin/admixtureWhiskerPlot.py admixture_columns.txt AdmixtureReferenceGroupMembershipBoxplot.pdf
    ${workflow.projectDir}/bin/admixtureWhiskerPlot.py admixture_reference_columns.txt AdmixtureGroupMembershipBoxplot.pdf

    # Determine which columns go forward by counting the column ratio in the reference file.
    # Only columns with a population greater than 10% of the size of the population will go forward in the pipeline.
    # TODO: Make sure the population is at least 1000 individuals in order to continue with the pipeline.
    awk '{for(i=1;i<=NF;i++) if(\$i>0.85) count[i]++} END{for(i in count) print i "," count[i] }' admixture_reference_columns.txt > referenceColumnCount.txt
    awk '{for(i=1;i<=NF;i++) if(\$i>0.85) count[i]++} END{for(i in count) print i "," count[i] }' admixture_columns.txt > columnCount.txt

    line_count=\$(wc -l < admixture_reference_columns.txt)
    # This line is the filter for 10% of the size of the overall population being able to continue the pipeline.
    line_count=\$((line_count / 10))
    
    awk -v lc="\$line_count" -F, '\$2 > lc {print \$1}' referenceColumnCount.txt > passingColumns.txt

    # Now that I have the passing columns, I need to get the IDs that pass the filter and remove the rest.
    while IFS= read -r col; do
        awk -v column="\$col" '\$column < 0.85 {print NR " " \$column}' admixture_reference_columns.txt > lines_to_remove.txt
        awk 'NR==FNR {exclude[\$1]; next} !(FNR in exclude)' lines_to_remove.txt line_numbers_reference_ids.txt > column\${col}_ids_filtered.txt
    done < passingColumns.txt

    """
}

process ADMIXTURE_STATS{

    input:
    path qFile
    val workflowName

    output:
    path 'admixture_statistics.xlsx', emit: stats

    publishDir "${workflowName}/AdmixtureStats", mode: "copy"

    script:
    """
    ${workflow.projectDir}/bin/admixtureStats.py < $qFile
    """
}

// Run this workflow to update the annotations in the input dataset.
// Updating annotations requires the proper liftover file from the annotations in the dataset to the requested annotation version
workflow UPDATE_ANNOTATIONS{
    // Report same (two old locations mapped to the same new location) or no physical location after update

    take:
    pedFile
    mapFile
    bimFile
    bedFile
    famFile
    liftover

    main:
    workflowName = "UpdateAnnotations"

    PLINK_HG38_UPDATE( pedFile, mapFile, bimFile, bedFile, famFile, liftover, workflowName )
    CHECK_NUM_INDIVIDUALS( PLINK_HG38_UPDATE.out.bedFile, PLINK_HG38_UPDATE.out.bimFile, PLINK_HG38_UPDATE.out.famFile, workflowName )
    CONVERT_TO_PLINK_FORMAT( CHECK_NUM_INDIVIDUALS.out.bedFile, CHECK_NUM_INDIVIDUALS.out.bimFile, CHECK_NUM_INDIVIDUALS.out.famFile, workflowName )
    //LIFTOVER_HG38( CHECK_NUM_INDIVIDUALS.out.bedFile, CHECK_NUM_INDIVIDUALS.out.bimFile, CHECK_NUM_INDIVIDUALS.out.famFile, liftover )
    LIFTOVER_PLINK_HG38( CONVERT_TO_PLINK_FORMAT.out.mapFile, CONVERT_TO_PLINK_FORMAT.out.pedFile, liftover, workflowName)
    CONVERT_TO_PLINK_BINARY_FORMAT( LIFTOVER_PLINK_HG38.out.mapFile, LIFTOVER_PLINK_HG38.out.pedFile, workflowName )
    //POST_LIFTOVER_PLINK( CONVERT_TO_PLINK_BINARY_FORMAT.out.bedFile, CONVERT_TO_PLINK_BINARY_FORMAT.out.bimFile, CONVERT_TO_PLINK_BINARY_FORMAT.out.famFile, LIFTOVER_HG38.out.snpsFile, LIFTOVER_HG38.out.posFile)
    
    //CONVERT_TO_PLINK_BINARY_FORMAT.out[0].view()
    emit:
    bedFile = CONVERT_TO_PLINK_BINARY_FORMAT.out.bedFile
    bimFile = CONVERT_TO_PLINK_BINARY_FORMAT.out.bimFile
    famFile = CONVERT_TO_PLINK_BINARY_FORMAT.out.famFile
    
}

process CHECK_NUM_INDIVIDUALS{
    container plinkContainer
    // If there are more than 5000 individuals, that's too many to process in a timely manner.  Grab 5000 random individuals.
    input:
    path bimFile
    path bedFile
    path famFile
    val workflowName

    output:
    path 'checkedIndividuals.bed', emit: bedFile
    path 'checkedIndividuals.bim', emit: bimFile
    path 'checkedIndividuals.fam', emit: famFile

    script:
    baseFile = famFile.baseName

    """
    # Count the number of lines in the input file
    num_lines=\$(wc -l < "$famFile")

    if [ "\$num_lines" -ge 5000 ]; then
        shuf -n 5000 $famFile | awk '{print \$1, \$2}' > filtered.fam
        ${params.plinkCommand} --bfile $baseFile --allow-extra-chr --keep filtered.fam --make-bed --out checkedIndividuals
    else
        mv $bimFile checkedIndividuals.bim
        mv $bedFile checkedIndividuals.bed
        mv $famFile checkedIndividuals.fam
    fi

    """
    
}

workflow FILTER_INDIVIDUALS{
    // Individuals Workflow:
    // Missing rate >= 5%
    // Mendelian error rate >= 1%
    // Individuals from families with known CNV or non-idiopathic cases
    // Admixture
    // Duplicates or relatedness

    take:
    bedFile
    bimFile
    famFile

    main:
    workflowName = "FilterIndividuals"
    //FILTER_UNAFFTECTED_AND_NO_PARENTS( HARDY_WEINBERG_EQUILIBRIUM.out.bedFile, HARDY_WEINBERG_EQUILIBRIUM.out.bimFile, HARDY_WEINBERG_EQUILIBRIUM.out.famFile)
    INDIVIDUALS_MISSINGNESS(bedFile, bimFile, famFile, workflowName, "")
    INDIVIDUALS_MENDELIAN_ERROR( INDIVIDUALS_MISSINGNESS.out.bedFile, INDIVIDUALS_MISSINGNESS.out.bimFile, INDIVIDUALS_MISSINGNESS.out.famFile, workflowName, "" )
    //INDIVIDUALS_CNV_FAMILIES_OR_IDIOPATHIC()


    //FILTER_DUPLICATES_OR_RELATEDNESS()

    emit:
    bedFile = INDIVIDUALS_MENDELIAN_ERROR.out.bedFile
    bimFile = INDIVIDUALS_MENDELIAN_ERROR.out.bimFile
    famFile = INDIVIDUALS_MENDELIAN_ERROR.out.famFile
}

workflow ADMIXTURE_PIPELINE {
    take:
    bedFile
    bimFile
    famFile
    referenceBedFile
    referenceBimFile
    referenceFamFile
    liftover
    
    main:
    workflowName = "AdmixturePipeline"
    // Prune the dataset for Admixture.  The pruned dataset is discarded for further analysis
    //CONVERT_TO_PLINK_FORMAT( referenceBedFile, referenceBimFile, referenceFamFile )
    //LIFTOVER_PLINK_HG38( CONVERT_TO_PLINK_FORMAT.out.mapFile, CONVERT_TO_PLINK_FORMAT.out.pedFile, liftover )
    //CONVERT_TO_PLINK_BINARY_FORMAT( LIFTOVER_PLINK_HG38.out.mapFile, LIFTOVER_PLINK_HG38.out.pedFile )
    //ADMIXTURE_JOIN_REFERENCE( bedFile, bimFile, famFile, CONVERT_TO_PLINK_BINARY_FORMAT.out.bedFile, CONVERT_TO_PLINK_BINARY_FORMAT.out.bimFile, CONVERT_TO_PLINK_BINARY_FORMAT.out.famFile)
    ADMIXTURE_JOIN_REFERENCE( bedFile, bimFile, famFile, referenceBedFile, referenceBimFile, referenceFamFile, workflowName )
    ADMIXTURE_PRUNE( ADMIXTURE_JOIN_REFERENCE.out.bedFile, ADMIXTURE_JOIN_REFERENCE.out.bimFile, ADMIXTURE_JOIN_REFERENCE.out.famFile, workflowName )
    ADMIXTURE( ADMIXTURE_PRUNE.out.bedFile, ADMIXTURE_PRUNE.out.bimFile, ADMIXTURE_PRUNE.out.famFile, Channel.of(5), workflowName )

    // Extract the individuals from each group with a cutoff of 85% confidence in the categorization
    ADMIXTRUE_EXTRACT_FROM_GROUP( ADMIXTURE.out.qFile, ADMIXTURE_JOIN_REFERENCE.out.famFile, famFile, referenceFamFile, workflowName )

    //ADMIXTURE_STATS( ADMIXTURE.out.qFile )
    emit:
    idsByColumn = ADMIXTRUE_EXTRACT_FROM_GROUP.out.filteredIdsByColumn
    
}

workflow FILTER_MARKERS{
    // Markers Workflow:
    
    // Indel markers
    // AT, TA, CG and GC markers removed
    // Missing rate >= 5%
    // Mendelian error rate >= 1%

    take:
    bedFile
    bimFile
    famFile

    main:
    workflowName = "FilterMarkers"
    SAME_LOCATION_FILTER( bedFile, bimFile, famFile, workflowName)
    INDEL_MARKER_FILTER( SAME_LOCATION_FILTER.out.bedFile, SAME_LOCATION_FILTER.out.bimFile, SAME_LOCATION_FILTER.out.famFile, workflowName )
    AT_GC_FILTER( INDEL_MARKER_FILTER.out.bedFile, INDEL_MARKER_FILTER.out.bimFile, INDEL_MARKER_FILTER.out.famFile, workflowName )
    MARKER_MISSINGNESS( AT_GC_FILTER.out.bedFile, AT_GC_FILTER.out.bimFile, AT_GC_FILTER.out.famFile, workflowName, "" )
    MENDELIAN_ERROR_RATE( MARKER_MISSINGNESS.out.bedFile, MARKER_MISSINGNESS.out.bimFile, MARKER_MISSINGNESS.out.famFile, workflowName )
    //MINOR_ALLELE_FREQUENCY( MENDELIAN_ERROR_RATE.out.bedFile, MENDELIAN_ERROR_RATE.out.bimFile, MENDELIAN_ERROR_RATE.out.famFile )
    //HARDY_WEINBERG_EQUILIBRIUM( MINOR_ALLELE_FREQUENCY.out.bedFile, MINOR_ALLELE_FREQUENCY.out.bimFile, MINOR_ALLELE_FREQUENCY.out.famFile )

    emit:
    bedFile = MENDELIAN_ERROR_RATE.out.bedFile
    bimFile = MENDELIAN_ERROR_RATE.out.bimFile
    famFile = MENDELIAN_ERROR_RATE.out.famFile
}

workflow POST_ADMIXTURE_PIPELINE {
    take:
    bedFile
    bimFile
    famFile
    admixtureIdsByColumn
    //TODO: Update the publishDir of the methods used by this pipeline.  Getting collisions.
    main:
    workflowName = "PostAdmixture"
    PLINK_FILTER_BY_IDS( bedFile, bimFile, famFile, admixtureIdsByColumn, workflowName )
    MARKER_MISSINGNESS( PLINK_FILTER_BY_IDS.out.bedFile, PLINK_FILTER_BY_IDS.out.bimFile, PLINK_FILTER_BY_IDS.out.famFile, workflowName, PLINK_FILTER_BY_IDS.out.column )
    MINOR_ALLELE_FREQUENCY( MARKER_MISSINGNESS.out.bedFile, MARKER_MISSINGNESS.out.bimFile, MARKER_MISSINGNESS.out.famFile, workflowName, MARKER_MISSINGNESS.out.column )
    MAF_AND_MISSING( MINOR_ALLELE_FREQUENCY.out.bedFile, MINOR_ALLELE_FREQUENCY.out.bimFile, MINOR_ALLELE_FREQUENCY.out.famFile, workflowName, MINOR_ALLELE_FREQUENCY.out.column )
    HARDY_WEINBERG_EQUILIBRIUM(MAF_AND_MISSING.out.bedFile, MAF_AND_MISSING.out.bimFile, MAF_AND_MISSING.out.famFile, workflowName, MAF_AND_MISSING.out.column )
    INDIVIDUALS_MISSINGNESS( HARDY_WEINBERG_EQUILIBRIUM.out.bedFile, HARDY_WEINBERG_EQUILIBRIUM.out.bimFile, HARDY_WEINBERG_EQUILIBRIUM.out.famFile, workflowName, HARDY_WEINBERG_EQUILIBRIUM.out.column )
    INDIVIDUALS_MENDELIAN_ERROR( INDIVIDUALS_MISSINGNESS.out.bedFile, INDIVIDUALS_MISSINGNESS.out.bimFile, INDIVIDUALS_MISSINGNESS.out.famFile, workflowName, INDIVIDUALS_MISSINGNESS.out.column )
    CHECK_SEX( INDIVIDUALS_MENDELIAN_ERROR.out.bedFile, INDIVIDUALS_MENDELIAN_ERROR.out.bimFile, INDIVIDUALS_MENDELIAN_ERROR.out.famFile, workflowName, INDIVIDUALS_MENDELIAN_ERROR.out.column )
        //SEX_INCONSISTENCY() //TODO: Implement this filter
    INDIVIDUALS_AUTOSOMAL_HETEROZYGOSITY( CHECK_SEX.out.bedFile, CHECK_SEX.out.bimFile, CHECK_SEX.out.famFile, workflowName, CHECK_SEX.out.column )
    KING( INDIVIDUALS_AUTOSOMAL_HETEROZYGOSITY.out.bedFile, INDIVIDUALS_AUTOSOMAL_HETEROZYGOSITY.out.bimFile, INDIVIDUALS_AUTOSOMAL_HETEROZYGOSITY.out.famFile, workflowName, INDIVIDUALS_AUTOSOMAL_HETEROZYGOSITY.out.column )
    FILTER_KING( KING.out.kingOut, workflowName )
}

workflow {
    //if( params.convertToHg38 )
    //    CONVERT_REFERENCE_GENOME

    mapFile = params.dataDir + params.dataName + ".map"
    pedFile = params.dataDir + params.dataName + ".ped"
    bedFile = params.dataDir + params.dataName + ".bed"
    bimFile = params.dataDir + params.dataName + ".bim"
    famFile = params.dataDir + params.dataName + ".fam"
    liftover = params.liftoverChainLocation + params.convertFrom //  params.dataDir + params.liftoverChain

    // TODO: Parameterize this and make it optional.
    pedChannel = Channel.fromPath( pedFile )
    mapChannel = Channel.fromPath( mapFile )
    bimChannel = Channel.fromPath( bimFile )
    bedChannel = Channel.fromPath( bedFile )
    famChannel = Channel.fromPath( famFile )

    referenceBed = Channel.fromPath( params.referenceDataDir + params.referenceDataName + ".bed" )
    referenceBim = Channel.fromPath( params.referenceDataDir + params.referenceDataName + ".bim" )
    referenceFam = Channel.fromPath( params.referenceDataDir + params.referenceDataName + ".fam" )

    UPDATE_ANNOTATIONS( pedChannel, mapChannel, bimChannel, bedChannel, famChannel, Channel.fromPath( liftover ) )
    FILTER_MARKERS( UPDATE_ANNOTATIONS.out.bedFile, UPDATE_ANNOTATIONS.out.bimFile, UPDATE_ANNOTATIONS.out.famFile )
    FILTER_INDIVIDUALS( FILTER_MARKERS.out.bedFile, FILTER_MARKERS.out.bimFile, FILTER_MARKERS.out.famFile )
    ADMIXTURE_PIPELINE( FILTER_INDIVIDUALS.out.bedFile, FILTER_INDIVIDUALS.out.bimFile, FILTER_INDIVIDUALS.out.famFile, referenceBed, referenceBim, referenceFam, Channel.fromPath( liftover ) )

    // Have to collect the results of the other process so that the multiple outputs of the admixture pipeline can be executed by multiple tasks.
    POST_ADMIXTURE_PIPELINE( FILTER_INDIVIDUALS.out.bedFile.collect(), FILTER_INDIVIDUALS.out.bimFile.collect(), FILTER_INDIVIDUALS.out.famFile.collect(), ADMIXTURE_PIPELINE.out.idsByColumn.flatten() )


    //POST_ADMIXTURE_ANALYSIS( ADMIXTURE_PIPELINE)
    //CREATE_GRAPHS() //TODO: Add running of runPlotMaker.sh to the workflow.
}
