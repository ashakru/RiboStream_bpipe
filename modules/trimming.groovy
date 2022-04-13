import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")
load (prog_dir + "/utils/utils.groovy")

trimGalore = {
  doc title: "trimGalore",
    desc:  """
      Description - quality based and adapter trimming using trim_galore
        input:  input - fastq/*.fastq.gz
          TRIM_ERROR (max.allowed error-rate), MIN_LENGTH (discard reads of length below this): global variables
          phredScore, adapter, trimQuality - branch variables
          adapterPE, trimQualityPE - branch variables for paired-end sequences
        output: output - fastq/*_trim.fastq.gz
        """, constraints: "constraints"

  output.dir = "${home_dir}/fastq"
  // def outputs = input.prefix.prefix + "_trimmed.fq.gz"
  def outputs = input.prefix.prefix + "_trim.fastq.gz"
  def input_pref=input.prefix.prefix

	produce(outputs){
    //set adapter to the value defined in .json, or leave it blank
		def adapterStr=""
      if (adapter != ""){
        adapterStr = " -a " + adapter.replaceAll(';',' -a ')
      } else {
        def adapterStr2=""
      }

      if (NEXTERA){
        adapterStr2 = " --nextera "
      } else if (ILLUMINA) {
        adapterStr2 = " --illumina "
        } else if (SMALL_RNA){
        adapterStr2 = " --small_rna "
				}

    exec """
      trim_galore --phred33 ${adapterStr2} --fastqc -e ${TRIM_ERROR} --length ${MIN_LENGTH} -q ${trimQuality} ${adapterStr} -o ${output.dir} ${input} &&
      mv ${input_pref}_trimmed.fq.gz ${input_pref}_trim.fastq.gz && mv ${input_pref}_trimmed_fastqc* fastqc/
      """ , "trim"
	}
}

trim = segment{
    if(TRIM_GALORE){
        trimGalore
    }
}
