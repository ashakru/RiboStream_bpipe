import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

tximport = {
	doc title: "tximport",
                desc:  """
                        Description - Transcript level counting with Tximport
                """,
                constraints: "constraints"

	if (experiment == "RNA"){
		folder_check('RNA_transcriptome')
		def gtf = "/home/jak75/genomes/gencode.v29/gencode.v29.annotation.gtf"
		output.dir = "${home_dir}/RNA_transcriptome"
		def outputs = "RNAseq_tximport_tx.RData"

		produce(outputs){
			exec """Rscript ${prog_dir}/utils/tximport.R ${input} ${gtf}"""
 			forward input
		}

	} else {
	forward input
	}
}
