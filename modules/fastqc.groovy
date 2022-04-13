//Module to produce fastqc report for each sample
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

//Define pipeline stage: fastQC
fastQC = {
  doc title: "fastQC",
		desc:  """
			Description - perform quality report with FastQC
			Single-end:
				input:  input - fastq/*.fastq.gz
				output: output1 - fastq/*.fastq.gz
						    output2 - fastqc/*_fastqc.zip
		""",
    constraints: "constraints"

  folder_check('fastqc')
  output.dir =  "${home_dir}/fastqc"

	def outputs =  get_filename(input.prefix.prefix) + "_fastqc.zip"
		produce(outputs){
			exec """
				fastqc -o ${output.dir} --noextract -f fastq  ${input}
			"""
		}
		forward input, output
  }
