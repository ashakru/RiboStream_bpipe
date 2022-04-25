import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

merging = {

		if (module == "fastq"){
				output.dir = "${home_dir}/fastq"

				if (merge == ""){
					succeed "Branch succeeded - file merged with another one."
				} else if (merge == "None"){
					forward input
				} else{
					def outputs =  get_filename(input.prefix.prefix) + "_merged.fastq.gz"
					def filename = "\$(cat ${home_dir}/${ID}_merge.txt)"
					//String toMerge = new File(filename)
					produce(outputs){
						exec """
							cat ${filename} > ${output}
						"""
					}
				}
		}

		if (module == "bam"){
			if (merge == ""){
				succeed "Branch succeeded - file merged with another one."
			} else if (merge == "None"){
				forward input
			} else {
				def outputs =  get_filename(input.prefix) + "_merged.bam"
				output.dir = "${home_dir}/bam"
				def filename = "\$(cat ${home_dir}/${ID}_merge.txt)"
				def merge_var = "${merge}"
				produce(outputs){
					exec """
						samtools merge -f ${output} ${filename}
					""" , "samtools_new"

				}
			}
		}

}
