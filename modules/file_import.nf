/*
Import samples using a sample sheet
*/

// Enable typed processes
nextflow.preview.types = true

// ADAPTED FROM https://stackoverflow.com/questions/74039553/nextflow-rename-barcodes-and-concatenate-reads-within-barcodes

process concat_reads {
    tag "${sample_name}"
    label 'process_low'
    publishDir "${params.output_dir}/concat_reads", mode: 'copy', failOnError: true

    input:
    (sample_name, fastq_files, species, report_group): Tuple<String, Path, String, String>

    output:
    tuple(sample_name, file("${sample_name}.${extn}"), species, report_group)

    script:
    if( fastq_files.every { it.name.endsWith('.fastq.gz') } )
        extn = 'fastq.gz'
    else if( fastq_files.every { it.name.endsWith('.fastq') } )
        extn = 'fastq'
    else if( fastq_files.every { it.name.endsWith('.fq.gz') } )
        extn = 'fq.gz'
    else if( fastq_files.every { it.name.endsWith('.fq') } )
        extn = 'fq'
    else
        error "Concatentation of mixed filetypes is unsupported"

    """
    cat ${fastq_files} > "${sample_name}.${extn}"
    """
}

workflow parse_sample_sheet {
    
    take:
    fastq_dir: String
    sample_sheet: Path
    barcode_dir: Boolean
    
    main:
    // update to cover all possible fastq file extensions
    fastq_extns = [ '.fastq', '.fastq.gz' , '.fq', '.fq.gz' ] 

    // deal with the case that fastqs are located in folders named by barcode
    if( barcode_dir ) {
        
        channel.fromPath(sample_sheet)
            .splitCsv(header: true)
            .map{ row ->
                def full_path = fastq_dir + "/" + "${row.barcode}"
                def all_files = file(full_path).listFiles()
                def fastq_files = all_files.findAll { 
                    fn -> fastq_extns.find { fn.name.endsWith( it ) }
                }
                tuple(row.sample_name, fastq_files, row.species, row.report_group)
            }
            .concat_reads
            .set{concatenated_file_tuple}
    }
    // deal with the case that fastqs are all located in the same folder and named barcode01{something}.fq.gz
    else {
        channel.fromPath(sample_sheet)
            .splitCsv(header: true)
            .map{ row ->
                def fastq_files = file(fastq_dir)
                    .listFiles()
                    .findAll {
                        fn -> fn.name.startsWith( "${row.barcode}" ) && fastq_extns.find { 
                            fn.name.endsWith( it ) }
                    }
                tuple(row.sample_name, fastq_files, row.species, row.report_group)
            }
            .concat_reads
            .set{concatenated_file_tuple} 
    }

    emit:
    concatenated_file_tuple
}