def fastqParser(fileLoc, revComp=True):
    """Reads in a fastq formatted file and returns a list of list of reads"""

    fastq = open(fileLoc, 'r').readlines()
    reads =[]

    for i in range(0,len(fastq)-1):
        if (fastq[i][0] == '@') & (fastq[i+1][0] != '@'): # sometimes quality score string will start with @ symbol. This removes them.
            reads.append(fastq[i+1].rstrip())
    if revComp:
        reads = revCompList(reads)  # If necessary to get the reverse complement of the reads
    return reads
def revComp(read):
    """Returns the reverse complement of a single read"""

    comp=[]

    for j in read:
        if j == 'A': comp.append('T')
        elif j == 'T': comp.append('A')
        elif j == 'C': comp.append('G')
        elif j == 'G': comp.append('C')
        else: comp.append('N')
    comp.reverse()
    comp=''.join(comp)
    return(comp)
def revCompList(readList):
    """Takes a list of reads and and returns a new list of reverse complements"""

    newList = []
    for i in range(0,len(readList)):
        newList.append(revComp(readList[i]))
    return newList
def hitFinder(target, r1, r2):
    hits = []
    for i in range(0, len(r2)):
        if target in r2[i]:
            hits.append(r1[i])
    return hits
def readCounter(readList, ranMerLen):
    """Takes a list of string type reads. Builds a new list of lists:
    [(str)The full read, (int)number of times the read is found, (int)unique number of times the read is found]
     ranMerLen has to be at least 1, meaning at least one nt needs to get trimmer, or else it breaks"""

    counts = []

    for i in range(0, len(readList)):
        flag = False
        if len(readList[i]) <= ranMerLen: continue  # if for some strange reason the read is shorter than the random mer
        for j in range(0, len(counts)):
            if readList[i][0:-ranMerLen] == counts[j][0][0:-ranMerLen]:  # looks for matches
                counts[j][1] += 1
                counts[j][2].append(readList[i][-ranMerLen:])
                flag = True
        if not flag:
            counts.append([readList[i],1,[readList[i][-ranMerLen:]]])  # Also keeps list of unique random mers
        if i%1000 == 0: print (round(i/len(readList)*100),"%")  # Progress bar
    for i in range(0,len(counts)):
        #counts[i][2] = len(set(counts[i][2]))  #Exact matches: Converts ranmer list to a set and returns len to get unique reads
        counts[i][2] = uniqueFinder(counts[i][2])

    return counts
def score(seq1, seq2):
    """Takes two sequences and returns an int for how many bp are different between them"""

    score = len(seq1)

    for i in range(0, len(seq1)):
        if len(seq1)!= len(seq2):
            print('the fuck?')
            continue
        if seq1[i] == seq2[i]: score -= 1
    return score
def uniqueFinder(ranMerList, maxMispair=1):
    uniqueList = []
    ranMerList = list(set(ranMerList))
    for i in range(0,len(ranMerList)):
        flag = False
        if not uniqueList:
            uniqueList.append(ranMerList[i])
        else:
            for j in range(0,len(uniqueList)):
                if score(ranMerList[i], uniqueList[j]) <= maxMispair:
                    flag = True
            if not flag:
                uniqueList.append(ranMerList[i])
    return len(uniqueList)
def writeCSV(counts,ranMerLen, outLoc, name):
    """Writes counts to a CSV file with a header. ranMerLen has to be at least 1"""

    outFile = open(outLoc+name+".csv", 'w')

    outFile.write('Sequence,Total Reads,Unique Reads\n')  # header
    counts = sorted(counts, key=lambda x:x[1], reverse=True)  # sorts by total reads term
    for i in range(0, len(counts)):
        outFile.write(counts[i][0][:-ranMerLen]+','+str(counts[i][1])+','+str(counts[i][2])+'\n')



#######################################################################

r1Loc1 = "C:\\Users\\Tim\\Desktop\\RawData\\041516 Gene Specific and RIPs-29649369\\Rea3-35047762\\Data\\Intensities\\BaseCalls\\Rea3_S5_L001_R1_001.fastq"
r2Loc2 = "C:\\Users\\Tim\\Desktop\\RawData\\041516 Gene Specific and RIPs-29649369\\Rea3-35047762\\Data\\Intensities\\BaseCalls\\Rea3_S5_L001_R2_001.fastq"
outLoc = "C:\\Users\\Tim\Desktop\\ProcessedData\\705\\"
manifestLoc = "C:\\Users\\Tim\\Desktop\\RawData\\16-04-15_703.csv"
r1 = fastqParser(r1Loc1)
print('r1 read in success!')
r2 = fastqParser(r2Loc2, revComp=False)
print('r2 read in success!')
assert len(r1) == len(r2)


manifest = open(manifestLoc, 'r')
for line in manifest:
    line=line.split(",")
    name = line[0]
    barcode = line[2].rstrip() + line[3].rstrip()
    ranMerLen = int(line[4]) + 2
    print(name, barcode, ranMerLen)
    hits = hitFinder(barcode, r1, r2)
    print(hits[1],"\n",hits[2],'\n',hits[3])
    counts = readCounter(hits, ranMerLen)
    writeCSV(counts, ranMerLen, outLoc,name)

'''
hits = hitFinder('AAAGGGACTGAACAAGCTGCTGGGC', r1, r2)
print(len(hits))
counts = readCounter(hits, 12)
writeCSV(counts, 12, outLoc)
'''
