//Module to produce fastqc report for each sample
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

fastqdump = {
  doc title: "fastq-dump",
		desc:  """
			Description - downloading FASTQ files from SRA
				input:  SRA accession number
				output: output - fastq/*.fastq.gz
		""",
    constraints: "constraints"

  folder_check('fastq')

  output.dir =  "${home_dir}/fastq"
  def outputs = "${ID}"+".fastq.gz"

	produce("${ID}.fastq.gz"){
		exec """
			fastq-dump --stdout --gzip ${ID} > ${outputs}
		"""

    exec """
			mv ${outputs} ${output.dir}/${ID}.fastq.gz
		"""
		}

  }
