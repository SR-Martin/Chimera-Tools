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
		self.read_id 			= read_id
		self.read_length 		= read_length
		self.alignments 		= list()
		self.alignment_coverage = 0
		self.alignment_overlap	= 0
		self.alignment_sum		= 0

	def printAlignmentReportForRead(self):
		print "--------------------------------------------------"
		print "Chimera suspected in read " + self.read_id + ", consisting of " + str(len(self.alignments)) + " parts:"
		for i in xrange(0, len(self.alignments)):
			print "-----Alignment " + str(i+1) + "----- "
			self.alignments[i].printAlignment()
		print "Total length: " + str(self.alignment_sum) + " / " + str(self.read_length) + " (" + str(self.alignment_coverage * 100) + ")."		

start = time.time()

pafFileName = ''
minAlignmentLength = 100
minQuality = 15
minCoverage = 95
reportChimeras = False
fastaFileName = ''
outputFastaName = ''
writeNewFasta = False
printAlignmentStats = False
polyChimera = 0

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:q:c:l:f:o:")
except getopt.GetoptError:
	print "Option not recognised."
	print "chimeraFinder.py -i <paf-file> -q <min alignment quality> -l <min alignment length> -c <min chimera coverage (%) -f <fasta containing original reads> -o <filename for new fasta with chimeras split into separate reads>"
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print "chimeraFinder.py -i <paf-file> -q <min alignment quality> -l <min alignment length> -c <min chimera coverage (%)  -f <fasta containing original reads> -o <filename for new fasta with chimeras split into separate reads>"
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

averagePercentOfReadsAligned = 0
for entry in reads.itervalues():
	entry.alignments.sort(key=lambda alignment: alignment.query_start)
	alignmentList = entry.alignments
	uniqueAlignedSum = 0
	lastQueryEnd = 0
	averageAlignmentLengthForRead = 0
	for alignment in alignmentList:
		averageAlignmentLengthForRead += alignment.alignment_length
		entry.alignment_sum += alignment.query_end - alignment.query_start
		if lastQueryEnd <= alignment.query_start:
			uniqueAlignedSum += (alignment.query_end - alignment.query_start)
		elif alignment.query_end >= lastQueryEnd :
			uniqueAlignedSum += (alignment.query_end - lastQueryEnd)
			entry.alignment_overlap += (lastQueryEnd - alignment.query_start)
		else:
			entry.alignment_overlap += (alignment.query_end - alignment.query_start)
		lastQueryEnd = alignment.query_end
	averageAlignmentLengthForRead = float(averageAlignmentLengthForRead) / len(entry.alignments)
	entry.alignment_coverage = float(uniqueAlignedSum) / entry.read_length
	percentOfReadAligned = entry.alignment_coverage * 100
	averagePercentOfReadsAligned += percentOfReadAligned

	if printAlignmentStats:
		print "--- Read " + entry.read_id + " stats ---"
		print "Number of alignments: " + str(len(entry.alignments))
		print "Sum of all alignment lengths: " + str(entry.alignment_sum)
		print "Average alignment length: " + str(averageAlignmentLengthForRead)
		print "Alignments overlapped by " + str(entry.alignment_overlap) + " bp."
		print "Aligned " + str(uniqueAlignedSum) + " / " + str(entry.read_length) + " (" + str(percentOfReadAligned) + "%) bp."

averagePercentOfReadsAligned = float(averagePercentOfReadsAligned)/len(reads)
print "\nAverage percent of reads aligned: " + str(averagePercentOfReadsAligned) + "%"

#look for chimerae
chimeras = dict()
for entry in reads.itervalues():
	if len(entry.alignments) > 1 and entry.alignment_overlap == 0 and entry.alignment_coverage * 100 > minCoverage:
		chimeras[entry.read_id] = Read(entry.read_id, entry.read_length)
		chimeras[entry.read_id].alignments = entry.alignments
		if reportChimeras:
			entry.printAlignmentReportForRead()
		if len(entry.alignments) > 2:
			polyChimera += 1

#free this memory
del reads

if not chimeras:
	print "Could not find Chimera's in file " + pafFileName + "."
	sys.exit()

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
				headerFields = ""
				for line in fastaFile:
					if count % 2 == 1:
						headerFields = line.split()
						readID = headerFields[0][1:]
						# Check to see if this read is in the list of chimeras 
						if readID in chimeras.keys():
							isChimera = True
						else:
							isChimera = False
							outputFile.write(line)
					else:
						if isChimera:
							read = chimeras[readID]
							lastSplittingPoint = 0
							for i in xrange(0, len(read.alignments)):
								header = ">" + readID + "_chimera_" + str(i+1)
								for j in xrange(1, len(headerFields)):
									header += " " + headerFields[j]
								header +="\n"
								outputFile.write(header)
								if i != len(read.alignments) - 1:
									assert(read.alignments[i].query_end <= read.alignments[i+1].query_start)
									splittingPoint = ((read.alignments[i].query_end + read.alignments[i+1].query_start) / 2) + 1
									assert(splittingPoint - lastSplittingPoint >= 100)
									outputFile.write(line[lastSplittingPoint:splittingPoint] + "\n")
									lastSplittingPoint = splittingPoint
								else:
									outputFile.write(line[lastSplittingPoint:])
						else:
							outputFile.write(line)	
					count += 1	
				numReads = count/2		
	except (OSError, IOError) as e: 
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print "Could not open file " + fastaFile
			sys.exit(2)

numChimeras = len(chimeras)
percent = (float(numChimeras) / numReads) * 100
print "Found " + str(numChimeras) + " chimera's out of " + str(numReads) + " reads (" + str(percent) + "%)."
print "Found " + str(polyChimera) + " poly chimera's. Amazing."

end = time.time()
print("chimeraFinder time: " + str(end - start))
