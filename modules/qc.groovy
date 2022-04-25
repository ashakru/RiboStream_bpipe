//Module that performs Ribo-Seq specific quality check with RIbowaltz package

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

qc = {
  	doc title: "Read counting and normalisation",
  		desc:  """
  			Description - Ribo-Seq QC
  			input:  bam/*_unique.bam
  			output: data/periodicity.RData
  		""",
  		constraints: "constraints"

      folder_check("counts")
      folder_check("scripts")

      def outputs = "data/periodicity.RData"
      output.dir = "${home_dir}/counts"

    	produce(outputs){
        exec """
        cp ${prog_dir}/utils/qc.R scripts/
        """

        exec """
        Rscript ${prog_dir}/utils/qc.R
        """
	}
}
