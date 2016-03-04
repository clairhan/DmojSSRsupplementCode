'''
Created on Feb 8, 2016

@author: PiTav
'''

import cPickle
import glob

fileList = glob.glob('*Or98a_2.fasta')
fastaWriter = open('Or98a_2.fasta','w')


translationTable = cPickle.load(open('/tigress/ANDOLFATTO/Andrew/clustering2/data/translationTable.p'))
#translationTable = cPickle.load(open('../../../Clustering/data/translationTable.p'))
degenerateTable = {'R':('A','G'),'Y':('C','T'),'S':('G','C'),'W':('A','T'),'K':('G','T'),'M':('A','C')}

def degenExpansion(seq):
	maskedSequence = ''.join([('NNN' if sum([1 if seq[x:x+3][i] in degenerateTable else 0 for i in range(0,3)]) > 1 else seq[x:x+3]) for x in range(0,len(seq),3)])
	seqOne = []
	seqTwo = []
	for i in range(0,len(seq)):
		if seq[i] in degenerateTable:
			seqOne.append(degenerateTable[seq[i]][0])
			seqTwo.append(degenerateTable[seq[i]][1])
		else:
			seqOne.append(seq[i])
			seqTwo.append(seq[i])
	seqOne = ''.join(seqOne)
	seqTwo = ''.join(seqTwo)
	return (seqOne,seqTwo)

def translate(seq):
	#return ''.join(['-' if (seq[x:x+3] not in translationTable) else translationTable[seq[x:x+3]] for x in range(0,len(seq),3)])
	return ''.join(['-' if ('N' in seq[x:x+3]) else translationTable[seq[x:x+3]] for x in range(0,len(seq),3)])

for currFile in fileList:
	fastaWriter.write('>' + currFile + '_1\n')
	currReader = open(currFile)
	currReader.readline()
	sequence = currReader.readline()
	sequence = sequence.rstrip('\n')
	tempSeq = degenExpansion(sequence)
	seqOne = tempSeq[0]
	seqTwo = tempSeq[1]

	fastaWriter.write(translate(seqOne) + '\n')
	fastaWriter.write('>' + currFile + '_2\n')
	fastaWriter.write(translate(seqTwo) + '\n')
