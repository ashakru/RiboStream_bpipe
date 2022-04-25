import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")




star = {
	doc title: "star",
		desc:  """
			Description - sequence alignment with STAR
			Single-end:
				input: 	fastq/*.fastq.gz
						STAR_THREADNUM (number of threads), STAR_REFERENCE (reference genome file): global variables
				output: bam/*.bam
		""",
		constraints: "constraints"

	def star_ref
	def star_gtf

	star_ref=GENOME
	star_gtf=GTF

	folder_check('bam')
	output.dir = "${home_dir}/bam"

	def outputs =  get_filename(input.prefix.prefix) + ".bam"
	def outPrefix = get_filename(input.prefix.prefix)

	if(experiment == "Ribo"){
		produce(outputs){
			exec """mkdir -p  bam/${outPrefix}"""

// Default
			exec """
				STAR --runThreadN ${STAR_THREADNUM} --genomeDir ${star_ref} --readFilesIn ${input} --outFileNamePrefix bam/${outPrefix}/ --outFilterType BySJout --alignSJoverhangMin 8 --outFilterMultimapNmax 9999 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --sjdbGTFfile ${star_gtf}  --alignEndsType EndToEnd
			""", "STAR"

      exec """
				mv ${output.dir}/${outPrefix}/Aligned.toTranscriptome.out.bam bam/${outPrefix}_tx.bam
			""", "localcmd"
			exec """
				mv ${output.dir}/${outPrefix}/Aligned.out.bam ${output}
			""", "localcmd"
		}
	}

	if(experiment == "RNA"){
			produce(outputs){
				exec """mkdir -p  bam/${outPrefix}"""
				exec """echo  ${clip3pNbases}"""

// Default
      exec """
				STAR --runThreadN ${STAR_THREADNUM} --genomeDir ${star_ref} --readFilesIn ${input} --outFileNamePrefix bam/${outPrefix}/ --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --sjdbGTFfile ${star_gtf} --alignEndsType EndToEnd
      """, "STAR"
			exec """
				mv ${output.dir}/${outPrefix}/Aligned.toTranscriptome.out.bam bam/${outPrefix}_tx.bam
			""", "localcmd"
			exec """
				mv ${output.dir}/${outPrefix}/Aligned.out.bam ${output}
			""", "localcmd"
		}
	}
}

alignment = segment{
    if(STAR){
	star
    }
}
