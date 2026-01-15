#!/usr/bin/env nextflow

//Enable typed processes
nextflow.preview.types = true

// Print help message information
workflow helpMessage {

    if ( params.help ) {
        log.info"""
Usage:  nextflow run main.nf --read_dir [path to fastq files] --sample_sheet [path to sample sheet] \\
		--phagemid_ref [path to reference genome] --matchbox_script [path to matchbox script]

Required Arguments:
--read_dir		: Specify full path of read file(s) location
--sample_sheet          : Specify location of the .csv sample sheet (format: barcode01, sample_x, rat, 1)
--phagemid_ref		: Specify location of the reference genome
--matchbox_script	: Specify matchbox script
--matchbox_parameters	: Specify parameters for matchbox script


Optional Arguments:
--output_dir		: Where the output files will be written to (default: "$projectDir/results")
--barcode_dir           : Whether the input fastq files are located within folders named barcode01 etc (default: false)
-profile		: Specify the profile to run nextflow through
			  Options - [standard, wehi, conda, singularity, local] (default: standard)
""".stripIndent()
    
    exit 1
    }
}