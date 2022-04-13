import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")


multiQC = {
	doc title: "multiQC",
		desc:  """
			Description - generate a summary of fastQC reports
			Single-end:
				input:  input - config/*.json
				output: output - multiqc_report.html in multiqc/SUBFOLDER/
		""",
		constraints: "constraints"

	folder_check('multiqc')

	currentDate = new Date().format( 'yyyy-MM-dd_HH-mm-ss' )
	output.dir =  "${home_dir}/multiqc/${module}_${currentDate}"

	if (module == "fastqc"){

		def outputs =  "multiqc_report.html"
		produce(outputs){
			exec """
				multiqc fastqc/*fastqc.zip -o ${output.dir}
			""", "multiqc"
			}
		}

		if (module == "trimmed"){
		def outputs =  "multiqc_report.html"
		produce(outputs){
			exec """
				multiqc fastqc/*trimmed_fastqc.zip -o ${output.dir}
			""", "multiqc"
		}
	}
	forward inputs
	}
