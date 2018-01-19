import sys, os, json, collections, re, msvcrt

def ShowUsage():
	print ('\nUsage: %s [compiled-demographics-file] [migration-type]' % os.path.basename(sys.argv[0]))

def CheckFiles(infilename, outfilename):
	print('Source file:  %s' % infilename)
	print('Destination file: %s' % outfilename)

	if not os.path.exists(infilename):
		print('Source file doesn\'t exist!')
		return False
	
	if os.path.exists(outfilename):
		print('Destination file already exists!  Overwriting...')

		try:
			os.remove(outfilename);
		except:
			print('Error while trying to delete file: %s' % sys.exc_info()[1])
			return False

	return True

def OrderedJsonLoad(in_dict):
	out_dict = collections.OrderedDict([])
	for pair in in_dict:
		out_dict[pair[0]] = pair[1]
	
	return out_dict

def GetMaxValueCountForMigrationType(mig_type):
	return {
		'sea' : 5,
		'local' : 100,
		'regional' : 30,
		'air' : 60
	}.get(mig_type, -1.0)


def main(tool, compiled_demographics, mig_type, outfilename=None):

	infilename = compiled_demographics
	print(infilename)
	if not outfilename:
		outfilename = re.sub('_demographics.compiled.json$', '_%s_migration.bin.json' % mig_type, infilename)
	print(outfilename)
	maxvalcount = GetMaxValueCountForMigrationType(mig_type)
	
	if maxvalcount <= 0:
		print('Invalid migration-type string: %s' % mig_type)
		exit(-1)

	if not CheckFiles(infilename, outfilename):
		exit(-1)

	with open(infilename, 'r') as file:
		demogjson = json.load(file, object_pairs_hook=OrderedJsonLoad) # TODO: should add error-messaging around this loading of the JSON

	migjson = collections.OrderedDict([])
	migjson['Metadata'] = demogjson['Metadata']
	migjson['Metadata']['Tool'] = os.path.basename(tool)
	migjson['Metadata']['DatavalueCount'] = maxvalcount
	
	strmap = demogjson['StringTable']

	demogoffsets = demogjson['NodeOffsets']
	
	migoffsets = ''
	nodecount = 0
	nodecountfound = 0

	for node in demogjson['Nodes']:
		if mig_type != "sea" or strmap['Seaport'] in node[strmap['NodeAttributes']] and node[strmap['NodeAttributes']][strmap['Seaport']] == 1:
			migoffsets += demogoffsets[nodecount*16: nodecount*16 + 8] + '%0.8X' % (nodecountfound * maxvalcount * 12) # 12 -> sizeof(uint32_t) + sizeof(double)
			nodecountfound += 1
		nodecount += 1

	migjson['Metadata']['NodeCount'] = nodecountfound
	migjson['NodeOffsets'] = migoffsets
			
	with open(outfilename, 'w') as file:
		print("MIGRATION HEADER PATH: " + outfilename)
		json.dump(migjson, file, indent=5)


if __name__ == "__main__":
	
	if len(sys.argv) != 3:
		ShowUsage()
		exit(0)
	
	main(sys.argv[0], sys.argv[1:])