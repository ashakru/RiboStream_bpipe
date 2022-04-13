//Module that performs umi_tools extraction and duplication removal
//Implemented tools: Umi tools 0.1

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

umi_extract = {
  	doc title: "umitools_extract",
  		desc:  """
  			Description - Unique Molecular Identifiers extraction
  			input:  fastq/*.fastq.gz
  					    bc-pattern (Type of molecular barcode used in data)
  			output: fastq/*_umi.fastq.gz
                umi.log
  		""",
  		constraints: "constraints"

      if (experiment == "Ribo"){
    		  output.dir = "${home_dir}/fastq"
		      def output = input.prefix.prefix + "_umi.fq.gz"
        	def logfile = input.prefix.prefix + "_umi.log"
    		  produce(output){

          	exec """
             	umi_tools extract --stdin=${input} --bc-pattern=${BCPATTERN} --log=${logfile} --stdout=${output}
    		""",'umi'
    		}
    	} else {
		forward input
	}
}

umi = segment{
  if(UMITOOLS){
    umi_extract
  }
}
