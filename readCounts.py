import sys
from multiprocessing import Pool

def fastqParser(fileLoc, revComp=True):
    """Reads in a fastq formatted file and returns a list of list of reads"""

    fastq = open(fileLoc, 'r').readlines()
    reads =[]

    for i in range(0,len(fastq)):
        if fastq[i][0] == '@':
            if fastq[i+1].rstrip() == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN': pass  # removes placeholders
            elif fastq[i+1][0] != '@':  # sometimes quality score string will start with @ symbol. This removes them.
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

def score(seq1, seq2):
    """Takes two sequences and returns an int for how many bp are different between them"""

    score = len(seq1)

    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]: score -= 1
    return score

def readCounter(readList, ranMerLen=13):
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
        counts[i][2] = len(set(counts[i][2]))  #Exact matches: Converts ranmer list to a set and returns len to get unique reads
        #readList[i][2] = uniqueFinder(readList[i][2]) #accepts mismatches instead

    return counts

def uniqueFinder(ranMerList, maxMispair=1):
    uniqueList = []

    for i in range(0,len(ranMerList)):
        flag = false
        if not uniqueList:
            uniqueList.append(ranMerList[i])
        else:
            for j in range(0, len(uniqueList)):
                if score(ranMerList[i], uniqueList[j]) <= maxMispair:
                    flag = True
                    break
        if not flag:
            uniqueList.append(ranMerList)
    return len(uniqueList)

def writeCSV(counts,ranMerLen, outLoc):
    """Writes counts to a CSV file with a header. ranMerLen has to be at least 1"""

    outFile = open(outLoc, 'w')

    outFile.write('Sequence,Total Reads,Unique Reads\n')  # header
    counts = sorted(counts, key=lambda x:x[1], reverse=True)  # sorts by total reads term
    for i in range(0, len(counts)):
        outFile.write(counts[i][0][:-ranMerLen]+','+str(counts[i][1])+','+str(counts[i][2])+'\n')

def combiner(pooledCounts):
    for i in range(0, len(pooledCounts)):
        if i == 0:
            master = pooledCounts[i][:]
            continue
        for j in range (0, len(pooledCounts[i])):
            seq = pooledCounts[i][j][0]
            for k in range(0, len(master)):
                if pooledCounts[i][j][0] == master[k][0]:
                    master[k][1] += pooledCounts[i][j][1]
                    master[k][2] += pooledCounts[i][j][2]
                else:
                    master.append(pooledCounts[i][j])
    return master


############################################################################

inLoc = "C:\\Users\\Tim\\Desktop\\siTOE-F3_S1_L001_R1_001.fastq"
outLoc = "C:\\Users\\Tim\\Desktop\\threadtest.csv"
ranMerLen = 13

if __name__=='__main__':
    #reads = fastqParser(inLoc)
    reads = ['GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'TCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCTA',\
             'TCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA',\
             'GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'GCTACGCCTGTCTGACCGTCGCTTAGTACCACGCGA', 'GCTACGCCTGTCTGGGCGTCGCTTAGTACCACGCGA',\
             'GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'TCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA', 'GCTACGCCTGTCTGAGCGTCGCTTAGTACCACGCGA']
    slice = int(len(reads)/4)
    print (slice)
    with Pool(processes = 8) as pool:
        pooledCounts = pool.map(readCounter, (reads[:slice],reads[slice:slice*2],reads[slice*2:slice*3], reads[slice*3:]))
    counts = combiner(pooledCounts)
    print(counts)
    #writeCSV(counts, ranMerLen, outLoc)