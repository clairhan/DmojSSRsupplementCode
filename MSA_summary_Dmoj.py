import numpy
import collections
import glob
import string

fileName = 'Or98a_2.fasta'
fastaReader = open(fileName)
alignmentWriter = open('Or98a_2_MSA.fa','w')
summaryAlignmentWriter = open('Or98a_2_MSA_summary.txt','w')

def condenseProteinSequences(seqHolder):
    # Takes a number of sequences and returns the consensus (being the most common), and omitting "-"
    # If two amino acids are equally common, the first one alphabetically is chosen
    consensusSequence = []
    assert len(seqHolder) > 1
    for ind in range(0,len(seqHolder[0])):
        aaList = numpy.array([seqHolder[x][ind] for x in range(0,len(seqHolder)) if "-" != seqHolder[x][ind]])
        counter = collections.Counter(iter(aaList)).most_common()
        if len(counter) == 0:
            consensusSequence.append('-')
        elif len(counter) > 1 and counter[0][1] == counter[1][1]:
            consensusSequence.append('-')
        else:
            consensusSequence.append(counter[0][0])
    return ''.join(consensusSequence)

sequenceHolder = dict()

currLine = fastaReader.readline()

while currLine != "":
    if currLine[0] == '>':
        tempKey = currLine[1:].rstrip('\n')
        tempValue = fastaReader.readline().rstrip('\n')
        sequenceHolder[tempKey] = tempValue
    currLine = fastaReader.readline()

sequenceKeys = ['G0-11-MAPQ20_allsites_' + fileName + '_1','G0-11-MAPQ20_allsites_' + fileName + '_2',
                'G1-11-MAPQ20_allsites_' + fileName + '_1','G1-11-MAPQ20_allsites_' + fileName + '_2',
                'G0-12-MAPQ20_allsites_' + fileName + '_1','G0-12-MAPQ20_allsites_' + fileName + '_2',
                'G1-12-MAPQ20_allsites_' + fileName + '_1','G1-12-MAPQ20_allsites_' + fileName + '_2',
                'G0-13-MAPQ20_allsites_' + fileName + '_1','G0-13-MAPQ20_allsites_' + fileName + '_2',
                'G1-13-MAPQ20_allsites_' + fileName + '_1','G1-13-MAPQ20_allsites_' + fileName + '_2',
                'G0-14-MAPQ20_allsites_' + fileName + '_1','G0-14-MAPQ20_allsites_' + fileName + '_2',
                'G1-14-MAPQ20_allsites_' + fileName + '_1','G1-14-MAPQ20_allsites_' + fileName + '_2',
                'G0-16-MAPQ20_allsites_' + fileName + '_1','G0-16-MAPQ20_allsites_' + fileName + '_2',
                'G1-16-MAPQ20_allsites_' + fileName + '_1','G1-16-MAPQ20_allsites_' + fileName + '_2',
                'G0-18-MAPQ20_allsites_' + fileName + '_1','G0-18-MAPQ20_allsites_' + fileName + '_2',
                'G1-18-MAPQ20_allsites_' + fileName + '_1','G1-18-MAPQ20_allsites_' + fileName + '_2']


A997bHaplotype1 = ['G0-11-MAPQ20_allsites_' + fileName + '_1','G1-11-MAPQ20_allsites_' + fileName + '_1',
                   'G0-12-MAPQ20_allsites_' + fileName + '_1','G1-12-MAPQ20_allsites_' + fileName + '_1']
A997bHaplotype2 = ['G0-11-MAPQ20_allsites_' + fileName + '_2','G1-11-MAPQ20_allsites_' + fileName + '_2',
                   'G0-12-MAPQ20_allsites_' + fileName + '_2','G1-12-MAPQ20_allsites_' + fileName + '_2']

OPNMhaplotype1 = ['G0-13-MAPQ20_allsites_' + fileName + '_1','G1-13-MAPQ20_allsites_' + fileName + '_1',
                  'G0-14-MAPQ20_allsites_' + fileName + '_1','G1-14-MAPQ20_allsites_' + fileName + '_1']
OPNMhaplotype2 = ['G0-13-MAPQ20_allsites_' + fileName + '_2','G1-13-MAPQ20_allsites_' + fileName + '_2',
                  'G0-14-MAPQ20_allsites_' + fileName + '_2','G1-14-MAPQ20_allsites_' + fileName + '_2']

SQ59aHaplotype1 = ['G0-16-MAPQ20_allsites_' + fileName + '_1','G1-16-MAPQ20_allsites_' + fileName + '_1',
                   'G0-18-MAPQ20_allsites_' + fileName + '_1','G1-18-MAPQ20_allsites_' + fileName + '_1']
SQ59aHaplotype2 = ['G0-16-MAPQ20_allsites_' + fileName + '_2','G1-16-MAPQ20_allsites_' + fileName + '_2',
                   'G0-18-MAPQ20_allsites_' + fileName + '_2','G1-18-MAPQ20_allsites_' + fileName + '_2']

sequenceHolder['A997b_Hap1'] = condenseProteinSequences([sequenceHolder[x] for x in A997bHaplotype1])
sequenceHolder['A997b_Hap2'] = condenseProteinSequences([sequenceHolder[x] for x in A997bHaplotype2])
sequenceHolder['OPNM_Hap1 '] = condenseProteinSequences([sequenceHolder[x] for x in OPNMhaplotype1])
sequenceHolder['OPNM_Hap2 '] = condenseProteinSequences([sequenceHolder[x] for x in OPNMhaplotype2])
sequenceHolder['SQ59a_Hap1'] = condenseProteinSequences([sequenceHolder[x] for x in SQ59aHaplotype1])
sequenceHolder['SQ59a_Hap2'] = condenseProteinSequences([sequenceHolder[x] for x in SQ59aHaplotype2])

for x in sequenceKeys:
    del sequenceHolder[x]

sequenceKeys = ['A997b_Hap1','A997b_Hap2','OPNM_Hap1 ','OPNM_Hap2 ','SQ59a_Hap1','SQ59a_Hap2']

consensusList = []
positionList = []
consensusHolder = dict()

for currKey in sequenceKeys:
    consensusHolder[currKey] = []

for ind in range(0,len(sequenceHolder[sequenceKeys[0]])):
    aaList = numpy.array([sequenceHolder[x][ind] for x in sequenceKeys])
    counter = collections.Counter(iter(aaList)).most_common()
    if len(counter) == len(sequenceKeys):
        consensusList.append('.')
        positionList.append(ind + 1)
    else:
        consensusList.append(counter[0][0])
        positionList.append(ind + 1)
    for currOrganism in sequenceKeys:
        if sequenceHolder[currOrganism][ind] != counter[0][0]:
            consensusHolder[currOrganism].append(sequenceHolder[currOrganism][ind])
        else:
            consensusHolder[currOrganism].append('.')

for currOrganism in consensusHolder:
    consensusHolder[currOrganism] = ''.join(consensusHolder[currOrganism])

consensusSeq = ''.join(consensusList)

for x in range(0,len(consensusList),80):
    alignmentWriter.write(''.join([' ' for y in range(0,14)]) + consensusSeq[x:x+80] + '\n')
    for currOrganism in sequenceKeys:
        alignmentWriter.write(currOrganism + ''.join([' ' for y in range(0,4)]) + consensusHolder[currOrganism][x:x+80] + '\n')
    alignmentWriter.write('\n')

alignmentWriter.close()

positionHolder = []
summaryHolder = dict()
summaryHolder['Reference'] = []
for currOrganism in sequenceKeys:
    summaryHolder[currOrganism] = []

for ind in range(0,len(consensusSeq)):
    aaList = numpy.array([consensusHolder[x][ind] for x in sequenceKeys])
    if not all(numpy.logical_or(aaList == '.',aaList == '-')):
        positionHolder.append(ind)
        for currOrganism in sequenceKeys:
            summaryHolder[currOrganism].append(consensusHolder[currOrganism][ind])
        summaryHolder['Reference'].append(consensusSeq[ind])

for x in summaryHolder:
    summaryHolder[x] = ''.join(summaryHolder[x])

positionHolder = numpy.array(positionHolder) + 1

for w in range(0,len(summaryHolder['Reference']),80):
    summaryAlignmentWriter.write(''.join([' ' for y in range(0,14)]) + ''.join([str(x)[0] if len(str(x)) == 3 else ' ' for x in positionHolder[w:w+80]]) + '\n')
    summaryAlignmentWriter.write(''.join([' ' for y in range(0,14)]) + ''.join([str(x)[1] if len(str(x)) == 3 else (str(x)[0] if len(str(x)) == 2 else ' ') for x in positionHolder[w:w+80]]) + '\n')
    summaryAlignmentWriter.write(''.join([' ' for y in range(0,14)]) + ''.join([str(x)[2] if len(str(x)) == 3 else (str(x)[1] if len(str(x)) == 2 else str(x)[0]) for x in positionHolder[w:w+80]]) + '\n')
    summaryAlignmentWriter.write('\n')
    summaryAlignmentWriter.write(''.join([' ' for y in range(0,14)]) + summaryHolder['Reference'][w:w+80] + '\n')
    for currOrganism in sequenceKeys:
        summaryAlignmentWriter.write(currOrganism + ''.join([' ' for y in range(0,4)]) + summaryHolder[currOrganism][w:w+80] + '\n')
    summaryAlignmentWriter.write('\n')

summaryAlignmentWriter.close()




