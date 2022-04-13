import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def home_dir=System.getenv("PWD")
def prog_dir=System.getenv("RIBOSTREAM_HOME")

load (home_dir + "/config.ini")

def folder_check(name){

	// Create a File object representing the folder 'A/B'
	def folder = new File( name )

	// If it doesn't exist, create all folders given in the path, including B
	if( !folder.exists() ) {
	  folder.mkdirs()
	}

}

def file_exist(fileName){
	def File file = new File(fileName);
	return file.exists()
}

def get_filename(filename) {

    def name = filename.split("/")[-1]

    return(name)

}

load_json = {
	doc title: "load_json",
		desc:  "Load the json file and set the branch variables",
		constraints: "input: jsonfile"

	def inputFile = new File(input)
	def inputJSON = new JsonSlurper().parseText(inputFile.text)

	inputJSON = inputJSON[0]

	if(module == "download"){
		branch.jsonFile = inputFile.absolutePath

		branch.ID = inputJSON.ID
		branch.merge = inputJSON.merge
		branch.pairedEnd = inputJSON.pairedEnd
		branch.phredScore = inputJSON.phredScore
		branch.adapter = inputJSON.adapter
		branch.replicate = inputJSON.replicate
		branch.condition = inputJSON.condition
		branch.experiment = inputJSON.experiment
		branch.readLength = inputJSON.readLength

		if(branch.trimQuality == ""){
			if(CUTADAPT == true | TRIM_GALORE == true){
				branch.trimQuality = TRIM_QUALITY
			} else if(TRIMMOMATIC == true){
				branch.trimQuality = REQUIRED_QUALITY
			}
		}

		if(branch.pairedEnd == "false"){
			output("${home_dir}/fastq/" + inputJSON.fileName)
		}else{
			branch.adapterPE = inputJSON.adapterPE
			branch.trimQualityPE = inputJSON.trimQualityPE
			if(branch.trimQualityPE == ""){
				if(CUTADAPT == true | TRIM_GALORE == true){
					branch.trimQualityPE = TRIM_QUALITY
				} else if(TRIMMOMATIC == true){
					branch.trimQualityPE = REQUIRED_QUALITY
				}
			}

			output("${home_dir}/fastq/" + inputJSON.fileName)
			output("${home_dir}/fastq/" + inputJSON.fileNamePE)
		}
	}

	if(module == "preproc"){
		branch.jsonFile = inputFile.absolutePath

		branch.ID = inputJSON.ID
		branch.merge = inputJSON.merge
		branch.pairedEnd = inputJSON.pairedEnd
		branch.phredScore = inputJSON.phredScore
		branch.adapter = inputJSON.adapter
		branch.trimQuality = inputJSON.trimQuality
		branch.replicate = inputJSON.replicate
		branch.condition = inputJSON.condition
		branch.experiment = inputJSON.experiment
		branch.readLength = inputJSON.readLength

		if(branch.trimQuality == ""){
			branch.trimQuality = TRIM_QUALITY
		}else{
			branch.trimQuality = inputJSON.trimQuality
		}

		if(branch.pairedEnd == "false"){
			output("${home_dir}/fastq/" + inputJSON.fileName)
		}else{
			branch.adapterPE = inputJSON.adapterPE
			branch.trimQualityPE = inputJSON.trimQualityPE
			if(branch.trimQualityPE == ""){
				branch.trimQualityPE = TRIM_QUALITY
			}else{
				branch.trimQualityPE = inputJSON.trimQualityPE
			}

			output("${home_dir}/fastq/" + inputJSON.fileName)
			output("${home_dir}/fastq/" + inputJSON.fileNamePE)
		}
	}

	if(module == "trim"){
		branch.jsonFile = inputFile.absolutePath

		branch.ID = inputJSON.ID
		branch.merge = inputJSON.merge
		branch.pairedEnd = inputJSON.pairedEnd
		branch.phredScore = inputJSON.phredScore
		branch.adapter = inputJSON.adapter
		branch.trimQuality = inputJSON.trimQuality
		branch.condition = inputJSON.condition
		branch.replicate = inputJSON.replicate
		branch.readLength = inputJSON.readLength
		branch.experiment = inputJSON.experiment

		if(branch.trimQuality == ""){
			branch.trimQuality = TRIM_QUALITY
		}else{
			branch.trimQuality = inputJSON.trimQuality
		}

		if(branch.pairedEnd == "false"){
			output("${home_dir}/fastq/" + inputJSON.fileName)
		}else{
			branch.adapterPE = inputJSON.adapterPE
			branch.trimQualityPE = inputJSON.trimQualityPE
			if(branch.trimQualityPE == ""){
				branch.trimQualityPE = TRIM_QUALITY
			}else{
				branch.trimQualityPE = inputJSON.trimQualityPE
			}

			output("${home_dir}/fastq/" + inputJSON.fileName)
			output("${home_dir}/fastq/" + inputJSON.fileNamePE)

		}
	}

	if(module == "alignment"){
		branch.jsonFile = inputFile.absolutePath

		branch.ID = inputJSON.ID
		branch.merge = inputJSON.merge
		branch.pairedEnd = inputJSON.pairedEnd
		branch.phredScore = inputJSON.phredScore
		branch.adapter = inputJSON.adapter
		branch.trimQuality = inputJSON.trimQuality
		branch.readLength = inputJSON.readLength
		branch.experiment = inputJSON.experiment

		if(branch.pairedEnd == "false"){
			output("${home_dir}/fastq/" + inputJSON.fileName)

		}else{
			branch.adapterPE = inputJSON.adapterPE
			branch.trimQualityPE = inputJSON.trimQualityPE

			output("${home_dir}/fastq/" + inputJSON.fileName)
			output("${home_dir}/fastq/" + inputJSON.fileNamePE)
		}
	}

	if(module == "merge"){

		branch.jsonFile = inputFile.absolutePath
		branch.ID = inputJSON.ID
		branch.merge = inputJSON.merge
		branch.pairedEnd = inputJSON.pairedEnd
		branch.phredScore = inputJSON.phredScore
		branch.adapter = inputJSON.adapter
		branch.trimQuality = inputJSON.trimQuality
		branch.experiment = inputJSON.experiment
		branch.readLength = inputJSON.readLength

		output("${home_dir}/fastq/" + inputJSON.fileName)
	}
	if(module == "experiment"){
    branch.ID = inputJSON.SampleID
		branch.filename = inputJSON.FileName
    branch.experiment = inputJSON.Experiment
		branch.condition = inputJSON.Condition
		branch.replicate = inputJSON.Replicate
		branch.merged = inputJSON.Merged

		output("${home_dir}/sampleSheet_alignment.csv")

	}
}

save_json = {
	doc title: "save_json",
		desc:  "Save a json file using the branch variables",
		constraints: "save_json.using(module: 'preproc')"

	def outputJSON = [:]

	if(module == "trim"){
		outputJSON.ID = ID
		outputJSON.fileName = get_filename("${input1}")
		outputJSON.merge = merge
		outputJSON.pairedEnd = pairedEnd
		outputJSON.phredScore = phredScore
		outputJSON.adapter = adapter
		outputJSON.trimQuality = trimQuality
		outputJSON.condition = condition
		outputJSON.replicate = replicate
		outputJSON.readLength = readLength
		outputJSON.experiment = experiment

		if(pairedEnd == "true"){
			outputJSON.fileNamePE = get_filename("${input2}")
			outputJSON.adapterPE = adapterPE
			outputJSON.trimQualityPE = trimQualityPE
		}
	}

	if(module == "alignment"){
		outputJSON.ID = ID
		outputJSON.fileName = get_filename("${input1}")
		outputJSON.merge = merge
		outputJSON.pairedEnd = pairedEnd
		outputJSON.phredScore = phredScore
		outputJSON.adapter = adapter
		outputJSON.trimQuality = trimQuality
		outputJSON.condition = condition
		outputJSON.replicate = replicate
		outputJSON.readLength = readLength
		outputJSON.experiment = experiment

		if(pairedEnd == "true"){
			outputJSON.fileNamePE = get_filename("${input2}")
			outputJSON.adapterPE = adapterPE
			outputJSON.trimQualityPE = trimQualityPE
		}
	}

	if(module == "merge"){
		outputJSON.ID = ID
		outputJSON.fileName = get_filename("${input1}")
		outputJSON.merge = merge
		outputJSON.pairedEnd = pairedEnd
		outputJSON.phredScore = phredScore
		outputJSON.adapter = adapter
		outputJSON.trimQuality = trimQuality
		outputJSON.experiment = experiment

		if(pairedEnd == "true"){
			outputJSON.adapterPE = adapterPE
			outputJSON.trimQualityPE = trimQualityPE
		}
	}

	if(module == "bam"){
		outputJSON.ID = ID
		outputJSON.fileName = "${home_dir}/bam/" + get_filename("${input1}")
		outputJSON.merge = merge
		outputJSON.pairedEnd = pairedEnd
		outputJSON.phredScore = phredScore
		outputJSON.adapter = adapter
		outputJSON.trimQuality = trimQuality
		outputJSON.experiment = experiment
	}

	def filePath = "configs/" + ID + "_" +  module + ".json"
	new File(filePath).write(new JsonBuilder([outputJSON]).toPrettyString())
	output(home_dir+"/"+filePath)
}

info = {
	println "INFO ${ID}: ${input}"
	// Print all avaible branch variables

	this.binding.variables.each{ key, value ->
		println "INFO ${ID}: ${key}"
	}

	forward input
}

generateCSV = {
	doc title: "generateCSV",
		desc:  """
			Description - save branch variables into .csv file
			input: all branch variables
			output: sampleSheet_${module}.csv
		""",
		constraints: "constraints"
	produce("sampleSheet_${module}.csv"){
		exec """
			${prog_dir}/utils/run_generateCSV.sh sampleSheet_${module}.csv ${inputs}
		""", "localcmd"
	}
}

generateTXT = {
	doc title: "generateTXT",
		desc:  """
			Description - save branch file names into .txt file (for multiQC)
			input: all branch variables
			output: multiqc_${module}.txt
		""",
		constraints: "constraints"
	produce("multiqc_${module}.txt"){
		exec """
			${prog_dir}/utils/run_generateTXT.sh ${module} ${inputs}
		""", "localcmd"
		forward inputs
	}
}

modifyJSON = {
	doc title: "modifyJSON",
		desc: "Modify .json files",
		constraints: "constraints"

	if(module == "merge" && merge != "" && merge != "None" ){
		def ids = "${merge}".split(";")
		def inputFile, inputJSON, filePath
		new File("${ID}_${module}.txt").withWriter { out ->
			ids.each { id ->
				filePath = "configs/" + id + "_" +  module + ".json"
				inputFile = new File(filePath)
				inputJSON = new JsonSlurper().parseText(inputFile.text)
				inputJSON = inputJSON[0]
				out.println "${home_dir}/fastq/" + inputJSON.fileName
			}
		}
		branch.merge = "${ID}_${module}.txt"
	}

	if(module == "trim" && TRIMMOMATIC=="true" ){
		def ids = "${adapter}".split(";") + "${adapterPE}".split(";")
		if(ids.size()>0){
			new File("${ID}_${module}.txt").withWriter { out ->
				ids.each { id ->
					out.println ">adapter"
					out.println id
				}
			}
			branch.adapter = "ILLUMINACLIP:${ID}_${module}.txt:${SEED_MISMATCHES}:${PALINDROME_CLIP_THRESHOLD}:${SIMPLE_CLIP_THRESHOLD}"
		}else{
			branch.adapter = ""
		}
	}
	forward inputs
}
