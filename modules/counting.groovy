//Module that performs umi_tools extraction and duplication removal
//Implemented tools: Umi tools 0.1

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

counting = {
  	doc title: "Read counting and normalisation",
  		desc:  """
  			Description - Ribo-Seq and RNA-Seq reads counting
  			input:  bam/*_unique.bam
  			output: counts/Exon_counts.csv,
                counts/FullCDS_counts.csv,
                counts/TrimCDS_counts.csv,
                counts/Counting_stat.csv,
                Counting_Normalising.html
  		""",
  		constraints: "constraints"

      folder_check("counts")
      folder_check("scripts")

      def outputs = "counts/Counting_stat.csv"
      output.dir = "${home_dir}/counts"

    	produce(outputs){
        exec """
        cp ${prog_dir}/utils/counting.R scripts/
        """

        exec """
        Rscript ${prog_dir}/utils/counting.R
        """
	}
}
