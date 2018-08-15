#!/usr/bin/python

import sys, getopt, random
import time
import errno

class Alignment:
	def __init__(self, query_start, query_end, target_start, target_end, target_sequence, alignment_length):
		self.query_start 		= query_start
		self.query_end 			= query_end
		self.target_start 		= target_start
		self.target_end 		= target_end
		self.target_sequence 	= target_sequence
		self.alignment_length	= alignment_length

	def printAlignment(self):
		print "Aligned to " + self.target_sequence + ", length " + str(self.alignment_length)
		print "Query Start: " + str(self.query_start) + ", Query End: " + str(self.query_end)
		print "Target Start: " + str(self.target_start) + ", Target End: " + str(self.target_end)

	def __lt__(self, other):
		return self.alignment_length < other.alignment_length 

class Read:
	def __init__(self, read_id, read_length):
		self.read_id 		= read_id
		self.read_length 	= read_length
		self.alignments 	= list()

start = time.time()

pafFileName = ''
minAlignmentLength = 100
minQuality = 15
minCoverage = 95
reportChimeras = False
fastaFileName = ''
outputFastaName = ''
writeNewFasta = False
printAlignmentStats = True

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:q:c:l:f:o:")
except getopt.GetoptError:
	print "Option not recognised."
	print "chimeraFinder.py -i <paf-file> -q <min alignment quality> -l <min alignment length> -c <min chimera coverage (%)>"
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print "chimeraFinder.py -i <paf-file> -q <min alignment quality> -l <min alignment length> -c <min chimera coverage (%)>"
		sys.exit()
	elif opt in ("-i"):
		pafFileName = arg
	elif opt in ("-q"):
		minQuality = int(arg)
	elif opt in ("-c"):
		minCoverage = int(arg)
	elif opt in ("-l"):
		minAlignmentLength = int(arg)
	elif opt in ("-f"):
		fastaFileName = arg
	elif opt in ("-o"):
		outputFastaName = arg
		writeNewFasta = True


reads = dict()
# try to open and read the file
try:
	with open(pafFileName, 'r') as pafFile:
		# organize the alignments, and filter out poor quality ones
		for line in pafFile:
			fields = line.split()
			alignmentLength = int(fields[3]) - int(fields[2])
			if alignmentLength > minAlignmentLength and fields[11] > minQuality:
				if not fields[0] in reads.keys():
					reads[fields[0]] = Read(fields[0], int(fields[1]))
				reads[fields[0]].alignments.append(Alignment(int(fields[2]), int(fields[3]), int(fields[7]), int(fields[8]), fields[5], alignmentLength))
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not open file " + pafFileName
		sys.exit(2)

estimatedNumberOfReads = len(reads)

if printAlignmentStats:
	averagePercentOfReadsAligned = 0
	for entry in reads.itervalues():
		alignmentList = sorted(entry.alignments, key=lambda alignment: alignment.query_start)
		uniqueAlignedSum = 0
		alignmentSum = 0
		overlapSum = 0
		lastQueryEnd = 0
		averageAlignmentLengthForRead= 0
		for alignment in alignmentList:
			averageAlignmentLengthForRead += alignment.alignment_length
			alignmentSum += alignment.query_end - alignment.query_start
			if lastQueryEnd <= alignment.query_start:
				uniqueAlignedSum += (alignment.query_end - alignment.query_start)
			elif alignment.query_end >= lastQueryEnd :
				uniqueAlignedSum += (alignment.query_end - lastQueryEnd)
				overlapSum += (lastQueryEnd - alignment.query_start)
			else:
				overlapSum += (alignment.query_end - alignment.query_start)
			lastQueryEnd = alignment.query_end
		averageAlignmentLengthForRead = float(averageAlignmentLengthForRead) / len(entry.alignments)
		percentOfReadAligned = float(uniqueAlignedSum) / entry.read_length * 100
		averagePercentOfReadsAligned += percentOfReadAligned
		print "--- Read " + entry.read_id + " stats ---"
		print "Number of alignments: " + str(len(entry.alignments))
		print "Sum of all alignment lengths: " + str(alignmentSum)
		print "Average alignment length: " + str(averageAlignmentLengthForRead)
		print "Alignments overlapped by " + str(overlapSum) + " bp."
		print "Aligned " + str(uniqueAlignedSum) + " / " + str(entry.read_length) + " (" + str(percentOfReadAligned) + "%) bp."
	averagePercentOfReadsAligned = float(averagePercentOfReadsAligned)/len(reads)
	print "----------------------------------------"
	print "Average percent of reads aligned: " + str(averagePercentOfReadsAligned) + "%"


#iterate through the dictionary and remove those entries that only have one alignment
entriesToDelete = list()
for entry in reads.itervalues():
	if len(entry.alignments) <= 1:
		entriesToDelete.append(entry.read_id)
	else:
		#order alignments by size
		entry.alignments.sort()

for read_id in entriesToDelete:
	del reads[read_id]
	
if printAlignmentStats:
	print "Reads with more than one alignment: " + str(len(reads)) + " / " + str(estimatedNumberOfReads) + "."

if not reads:
	print "Could not find Chimeras in file " + pafFileName + "."
	sys.exit()

#look for chimera's
chimeras = dict()
for entry in reads.itervalues():
	minLength = (minCoverage / 100.0) * entry.read_length
	besti = -1
	bestj = -1
	# compare all pairs of alignments
	# pick the pair that has the longest total length
	for i in xrange(len(entry.alignments)):
		for j in xrange(i + 1, len(entry.alignments)):
			if entry.alignments[i].query_end < entry.alignments[j].query_start or entry.alignments[j].query_end < entry.alignments[i].query_start:
				totalAlignmentLength = entry.alignments[i].alignment_length + entry.alignments[j].alignment_length
				if totalAlignmentLength > minLength:
					minLength = totalAlignmentLength
					besti = i
					bestj = j

	if besti != -1:
		chimeras[entry.read_id] = Read(entry.read_id, entry.read_length)
		chimeras[entry.read_id].alignments.append(entry.alignments[besti])
		chimeras[entry.read_id].alignments.append(entry.alignments[bestj])
		if reportChimeras:
			#report chimera:
			print "--------------------------------------------------"
			print "Chimera found in read " + entry.read_id + ":"
			print "-----Alignment 1----- "
			entry.alignments[0].printAlignment()
			print "-----Alignment 2----- "
			entry.alignments[1].printAlignment()
			print "Total length: " + str(entry.alignments[0].alignment_length + entry.alignments[1].alignment_length) + " / " + str(entry.read_length)

numReads = estimatedNumberOfReads

if writeNewFasta:
	try:
	# write the new fasta file with chimeras split into seperate reads
		with open(fastaFileName, 'r') as fastaFile:
			with open(outputFastaName, 'w') as outputFile:
				count = 1
				isChimera = False
				newHeaderLine1 = ''
				newHeaderLine2 = ''
				readID = ''
				for line in fastaFile:
					if count % 2 == 1:
						fields = line.split()
						readID = fields[0][1:]
						# Check to see if this read is in the list of chimeras 
						if readID in chimeras.keys():
							isChimera = True
							newHeaderLine1 = ">" + readID + "_chimera_1"
							newHeaderLine2 = ">" + readID + "_chimera_2"
							for i in xrange(1, len(fields)):
								newHeaderLine1 += " " + fields[i]
								newHeaderLine2 += " " + fields[i]
							newHeaderLine1 += "\n"
							newHeaderLine2 += "\n"
						else:
							isChimera = False
							outputFile.write(line)
					else:
						if isChimera:
							# find the splitting point
							# Probably should do something more intelligent than this - maybe keep the highest quality alignment intact and add the junk onto the other one?
							splittingPoint = 0
							alignment1 = chimeras[readID].alignments[0]
							alignment2 = chimeras[readID].alignments[1]
							if alignment1.query_end < alignment2.query_start:
								splittingPoint = int((alignment1.query_end + alignment2.query_start) / 2) + 1
							else:
								splittingPoint = int((alignment2.query_end + alignment1.query_start) / 2) + 1
							outputFile.write(newHeaderLine1)
							outputFile.write(line[:splittingPoint] + "\n")
							outputFile.write(newHeaderLine2)
							outputFile.write(line[splittingPoint:])
						else:
							outputFile.write(line)	
					count += 1	
				numReads = count - 1			
	except (OSError, IOError) as e: 
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print "Could not open file " + fastaFile
			sys.exit(2)

numChimeras = len(chimeras)
percent = (float(numChimeras) / numReads) * 100
print "Found " + str(numChimeras) + " chimera's out of " + str(numReads) + " reads (" + str(percent) + "%), in reads:"
for read in chimeras.itervalues():
	print read.read_id

end = time.time()
print("chimeraFinder time: " + str(end - start))
