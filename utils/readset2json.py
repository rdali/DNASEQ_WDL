
import json
import getopt
import os
import sys


def help():
        print '''
        ---------------------------------------------------------------------------
                                            HELP
        ---------------------------------------------------------------------------

        readset2json.py takes in the GenPipes readset file and outputs a json file
        with the same informafion for the WDL version of GenPipes.\n
        Make sure your readset is valid by using mugqicValidator.py.\n
        -r     --readset     myReadsetFile
        -o     --output      myJsonOutput
        -n     --name        workflowName
        -h     --help



        usage: python readset2json.py -r myReadsetFile -o myJsonOutput -n workflowName

        '''


def reads_to_dict(file_path):
	try:
		# read into rows
		with open(file_path, "r") as fl:
			rows = [line.split("\t") for line in [line.strip() for line in fl.readlines()] if line != "" and line is not None]

		keys = [key.lower() for key in rows[0]]
		rows = rows[1:]

		# all rows into samples
		data = [{key : row[index] if index < len(row) and row[index] != "" else "None" for index, key in enumerate(keys)} for row in rows]

		# build dict with sample : many reads
		wf_samples = {}
		for sample in data:
			if sample['sample'] not in wf_samples:
				wf_samples[sample['sample']]={'sample' : sample['sample'], 'readsets' : [{k:v for k, v in sample.items() if k != 'sample'}]}
			else:
				wf_samples[sample['sample']]['readsets'].append({k:v for k, v in sample.items() if k != 'sample'})

		#print(wf_samples)
		return {name + ".samples" : list(wf_samples.values())}


	except Exception as e:
		#log error
		print("error has occured in reads_to_dict while processing")
		return None


def write_json_to_file(data, output_path):
	try:
		# read into rows
		with open(output_path, "w") as fl:
			fl.write(json.dumps(data))
	except Exception as e:
		#log error
		print("error has occured in reads_to_dict while processing")




if __name__ == "__main__":

	readsetFile = outputFile = name = None

	try:
	    options,remainder = getopt.getopt(sys.argv[1:],'r:o:n:h',["readsetFile=", "outputFile=", "name=", "help"])
	except getopt.GetoptError, e:
	    print ("Error - "+str(e)+". See help ('-h' or '--help')")
	    sys.exit(2)


	for opt,arg in options:
	    if(opt in ["-r","--readset"]):
	        readsetFile=arg
	    elif(opt in ["-o","--output"]):
	        outputFile=arg
	    elif(opt in ["-n","--name"]):
	        name=arg
	    elif(opt in ["-h","--help"]):
	        help()
	        sys.exit()

	if readsetFile is None:
	    print ("Missing readsetFile!")
	    help()
	    sys.exit()

	if outputFile is None:
	    outputFile = os.path.splitext(readsetFile)[0]+ '.json'

	if name is None:
	    name = 'wf'



	data = reads_to_dict(readsetFile)

	if data is not None:
		write_json_to_file(data, outputFile)
