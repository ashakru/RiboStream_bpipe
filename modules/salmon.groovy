import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")


salmon = {
        doc title: "salmon",
                desc:  """
                        Description - Transcript level counting with Salmon
                """,
                constraints: "constraints"

	if (experiment == "RNA"){

		def salmon_tx

		if (GENOME=="hg38_gencodev29"){
        	        salmon_tx='/home/jak75/genomes/gencode.v29/gencode.v29.transcripts_headers.fa'
        	} else if (GENOME=="hg38_gencodev28"){
                	salmon_tx=''
   	     	}

		folder_check('bam')

    def outPrefix = get_filename(input.prefix)
		def outputs = get_filename(input.prefix) + "_quant.sf"
		output.dir = "${home_dir}/bam"

		produce(outputs){

			exec """salmon quant -t ${salmon_tx} -l SF --fldMean=49 --posBias --gcBias --seqBias --rangeFactorizationBins=4 --useEM -a bam/${outPrefix}_tx.bam -o bam/${outPrefix}""","salmon"
			exec """mv bam/${outPrefix}/quant.sf bam/${outPrefix}_quant.sf"""

			forward input
		}
	} else {
		forward input
	}

}
