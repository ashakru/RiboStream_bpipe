#!/usr/bin/env python

import sys
import os
import json
import csv

WDIR = os.getcwd()
RIBOSTREAM_HOME = os.getenv('RIBOSTREAM_HOME')
sys.path.append( RIBOSTREAM_HOME + '/utils')

#input: fastqc zip file, json file

# print("===============================================")

csvFileName = sys.argv[1]
jsonFiles = sys.argv[2:]

# print(jsonFiles)
# csvTrim = "sampleSheet_trim.csv"
# csvAlign = "sampleSheet_align.csv"

output = []
for jsonFile in jsonFiles:
    with open(jsonFile, 'r') as data_file:    
        jsonData = json.load(data_file)
        output.append(jsonData[0])

# print(output)
with open(csvFileName, 'w') as csvfile:
    header = ['ID', 'fileName', 'fileNamePE', 'experiment', 'condition', 'replicate', 
              'merge', 'pairedEnd', 'phredScore', 'adapter', 'trimQuality',
              'adapterPE', 'trimQualityPE', 'readLength']

    writer = csv.DictWriter(csvfile, fieldnames=header)

    writer.writeheader()
    writer.writerows(output)

# for row in output:
#     row['PhredScore'] = ""
#     row['Adapter'] = ""
#     row['trimQuality'] = ""

# with open(csvAlign, 'w') as csvfile:
#     header = ['ID', 'File', 'Experiment', 'Condition', 'Replicate',
#              'Lane', 'Merge', 'Input', 'PairedEnd', 'macsGroup',
#              'PhredScore', 'Adapter', 'trimQuality']
#     writer = csv.DictWriter(csvfile, fieldnames=header)

#     writer.writeheader()
#     writer.writerows(output)



