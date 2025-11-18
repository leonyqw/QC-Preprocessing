/*
Utilize minimap2 to align Oxford Nanopore DNA read sequences contained in fastq files against a reference genome database.
The aligned sequence is output as a SAM file.
*/

//Enable typed processes
nextflow.preview.types = true

process minimap2 {
	tag "${name}"

	// conda 'bioconda::minimap2'
	// conda (params.enable_conda ? 'bioconda::minimap2=2.30' : null)
	
	// Docker container for conda minimap2 (linux/amd64)
	// container "community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd"

	// Singularity container for conda minimap2 (linux/amd64)
	// container "oras://community.wave.seqera.io/library/minimap2:2.30--3bf3d6cb39a98dae"

    input:
	read_file: Path // Path for DNA sequence fastq files
	reference: Path // Path for reference genome
	name: String // Sample name
	
	output:
	aligned_read: Path = file("${name}_aligned.sam")

	/*
	Run minimap, mapping reads to a reference and outputs a sam file
	-a			Generates CIGAR and outputs alignments in sam format
	-x map-ont	Sets preset for ONT alignment
	-o			Output alignments to sam file
	*/
    script:
    """
	minimap2 \\
	-ax map-ont \\
	${reference} \\
	${read_file} \\
	-o "${name}_aligned.sam"
    """
}