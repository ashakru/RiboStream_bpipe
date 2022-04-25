import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

prealignment_bowtie = {doc title: "bowtie",
	desc:  """
		Description - sequence alignment with BOWTIE
		Single-end:
			input: 	fastq/*.fastq.gz
							BOWTIE2_THREADNUM (number of threads), GENOME: global variables
			output: bam/*.bam

	""",
	constraints: "constraints"



	if(experiment == "Ribo"){

		def bowtie_ref
	       	if (PREALIGNMENT == "rRNA"){
                	bowtie_ref='/home/jak75/genomes/bowtie_rRNA/index'
        	}

        	folder_check('bam')
	       	output.dir = "${home_dir}/fastq"
	
		def outputs =  get_filename(input.prefix.prefix) + "_prealigned.fastq.gz"
		def preoutput =  get_filename(input.prefix.prefix) + "_prealigned.fastq"
		def outPrefix = get_filename(input.prefix.prefix)

		produce(outputs){ 	
			
			exec """ mkdir -p bam/${outPrefix} """
			exec """
bowtie ${bowtie_ref} --best --norc -S --un ${home_dir}/fastq/${preoutput} -q ${input} > bam/${outPrefix}/Prealignment_rRNA.sam 2> bam/${outPrefix}/Prealignment_bowtie.Log		"""
			exec """gzip ${home_dir}/fastq/${preoutput} """
			}
	} else {
		forward input
	}	

}
prealignment = segment{
    if (PREALIGNMENT_TOOL == "BOWTIE"){
		prealignment_bowtie
		}
}
