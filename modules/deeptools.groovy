//Module that performs umi_tools extraction and duplication removal
//Implemented tools: Umi tools 0.1

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

bamCoverage = {
  	doc title: "deeptools bamCoverage",
  		desc:  """
  			Description - Generates coverage tracks from BAM file
  			input:  bam/*.bam
  			output: bigWigs/*.bw
  		""",
  	constraints: "constraints"

    folder_check("bamCoverage")

    output.dir = "${home_dir}/bamCoverage"
		def outputs = get_filename(input.prefix) + ".bw"

    produce(outputs){
      exec """
      bamCoverage --bam ${input} -o bamCoverage/${outputs} --binSize 1 --normalizeUsing CPM --ignoreForNormalization chrX chrY
    	"""
	}
    forward input
}

computeMatrix = {
  	doc title: "deeptools computeMatrix",
  		desc:  """
  			Description - Performs metagene analysis
  			input:  bam/*.bam
  			output: bigWigs/*.bw
  		""",
  	constraints: "constraints"

    folder_check("computeMatrix")

    output.dir = "${home_dir}/computeMatrix"
		def prefix = get_filename(input.prefix)
    def outputs = get_filename(input.prefix) + "metagene.gz"

    produce(outputs){
      exec """
      computeMatrix
      scale-regions -S bamCoverage/${prefix}.bw
      -R ${GTF}
      -o computeMatrix/${outputs}
      -m ${regionBodyLength}
      -bs ${binSize}
      --missingDataAsZero
      --skipZeros
      --samplesLabel ${prefix}
      -b ${beforeRegionStartLength}
      -a ${afterRegionStartLength}
      --unscaled5prime ${unscaled5prime}
      --unscaled3prime ${unscaled3prime}
      --metagene
      --exonID CDS
      --transcriptID transcript
    	"""
	}
    forward input
}

plotProfile = {
  doc title: "deeptools plotProfile",
    desc:  """
      Description - Performs metagene analysis
      input:  bam/*.bam
      output: bigWigs/*.bw
    """,
  constraints: "constraints"

  folder_check("plotProfile")

  output.dir = "${home_dir}/plotProfile"
  def prefix = get_filename(input.prefix)
  def outputs = get_filename(input.prefix) + "metagene.pdf"

  produce(outputs){
    exec """
    plotProfile
    -m  computeMatrix/${prefix}metagene.gz
    -out plotProfile/${outputs}
    """
}
  forward input
}
