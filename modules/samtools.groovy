//Module that performs umi_tools extraction and duplication removal
//Implemented tools: Umi tools 0.1

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

uniqueBam = {
  	doc title: "samtools view and samtools sort",
  		desc:  """
  			Description - Extract uniquely mapping reads
  			input:  bam/*.bam
  			output: bam/*_unique.bam
  		""",
  		constraints: "constraints"

      output.dir = "${home_dir}/bam"
      def prefix = input.prefix
      def outputs = input.prefix+ "_unique_sorted.bam"

    	produce(outputs){
        exec """samtools view -bh -q 255 ${input} > ${prefix}_unique.bam """
        exec """samtools sort ${prefix}_unique.bam ${prefix}_unique_sorted"""
	}

  forward output
}

indexBam = {
  	doc title: "samtools index",
  		desc:  """
  			Description - Index analysis ready bam file
  			input:  bam/*.bam
  			output: bam/*.bam.bai
  		""",
  		constraints: "constraints"

      output.dir = "${home_dir}/bam"
      def outputs = input + ".bai"

    	produce(outputs){
        exec """samtools index ${input} """
	}
  forward input
}
