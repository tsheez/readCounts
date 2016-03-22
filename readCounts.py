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

def readCounter(readList, ranMerLen=12):
    """Takes a list of string type reads. Builds a new list of lists:
    [(str)The full read, (int)number of times the read is found, (int)unique number of times the read is found]
     ranMerLen has to be at least 1, meaning at least one nt needs to get trimmer, or else it breaks"""

    counts = []

    for i in range(0, len(readList)):
        flag = False
        if len(readList[i]) <= ranMerLen: continue  # if for some strange reason the read is shorter than the random mer
        for j in range(0, len(counts)):
            if readList[i][0:-ranMerLen] == counts[j][0]:  # looks for matches
                counts[j][1] += 1
                counts[j][2].append(readList[i][-ranMerLen:])
                flag = True
                break
        if not flag:
            counts.append([readList[i][:-ranMerLen],1,[readList[i][-ranMerLen:]]])  # Also keeps list of unique random mers
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

def writeCSV(counts, outLoc):
    """Writes counts to a CSV file with a header. ranMerLen has to be at least 1"""

    outFile = open(outLoc, 'w')

    outFile.write('Sequence,Total Reads,Unique Reads\n')  # header
    counts = sorted(counts, key=lambda x:x[1], reverse=True)  # sorts by total reads term
    for i in range(0, len(counts)):
        outFile.write(counts[i][0]+','+str(counts[i][1])+','+str(counts[i][2])+'\n')

def splitter(reads):
    AA,AC,AG,AT = [],[],[],[]
    CA,CC,CG,CT = [],[],[],[]
    GA,GC,GG,GT = [],[],[],[]
    TA,TC,TG,TT = [],[],[],[]
    for i in range(0, len(reads)):
        if reads[i][:2] == 'AA': AA.append(reads[i])
        elif reads[i][:2] == 'AC': AC.append(reads[i])
        elif reads[i][:2] == 'AG': AG.append(reads[i])
        elif reads[i][:2] == 'AT': AT.append(reads[i])
        elif reads[i][:2] == 'CA': CA.append(reads[i])
        elif reads[i][:2] == 'CC': CC.append(reads[i])
        elif reads[i][:2] == 'CG': CG.append(reads[i])
        elif reads[i][:2] == 'CT': CT.append(reads[i])
        elif reads[i][:2] == 'GA': GA.append(reads[i])
        elif reads[i][:2] == 'GC': GC.append(reads[i])
        elif reads[i][:2] == 'GG': GG.append(reads[i])
        elif reads[i][:2] == 'GT': GT.append(reads[i])
        elif reads[i][:2] == 'TA': TA.append(reads[i])
        elif reads[i][:2] == 'TC': TC.append(reads[i])
        elif reads[i][:2] == 'TG': TG.append(reads[i])
        elif reads[i][:2] == 'TT': TT.append(reads[i])
        else: print('fudge-ums')
    return[AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT]

def splitter2(splits):
    out = []
    for i in range(0,len(splits)):
        tempA=[]
        tempC=[]
        tempG=[]
        tempT=[]
        for j in range(0,len(splits[i])):
            if splits[i][j][2] == 'A': tempA.append(splits[i][j])
            elif splits[i][j][2] == 'C': tempC.append(splits[i][j])
            elif splits[i][j][2] == 'G': tempG.append(splits[i][j])
            elif splits[i][j][2] == 'T': tempT.append(splits[i][j])
            else: print('Fudge-ums')
        out.append(tempA)
        out.append(tempC)
        out.append(tempG)
        out.append(tempT)
    return out
def splitter3(splits):
    out = []
    for i in range(0,len(splits)):
        tempA=[]
        tempC=[]
        tempG=[]
        tempT=[]
        for j in range(0,len(splits[i])):
            if splits[i][j][3] == 'A': tempA.append(splits[i][j])
            elif splits[i][j][3] == 'C': tempC.append(splits[i][j])
            elif splits[i][j][3] == 'G': tempG.append(splits[i][j])
            elif splits[i][j][3] == 'T': tempT.append(splits[i][j])
            else: print('Fudge-ums')
        out.append(tempA)
        out.append(tempC)
        out.append(tempG)
        out.append(tempT)
    return out

def betterSplitter(counts, pos):
    out = []

    if pos == 0: counts=[counts,]
    for i in range(len(counts)):
        tempA=[]
        tempC=[]
        tempG=[]
        tempT=[]
        for j in range(len(counts[i])):
            if counts[i][j][pos] == 'A': tempA.append(counts[i][j])
            elif counts[i][j][pos] == 'C': tempC.append(counts[i][j])
            elif counts[i][j][pos] == 'G': tempG.append(counts[i][j])
            elif counts[i][j][pos] == 'T': tempT.append(counts[i][j])
            else: print('Fudge-ums', counts[i][j])
        out.append(tempA)
        out.append(tempC)
        out.append(tempG)
        out.append(tempT)

    return out





############################################################################

inLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\2016-03-21-MiSeqRaw\\Tails032116-29290347\\siLuc-34379961\\Data\\Intensities\\BaseCalls\\siLuc_S2_L001_R1_001.fastq"
outLoc = "C:\\Users\\Tim\\Desktop\\test.csv"

if __name__=='__main__':
    reads = fastqParser(inLoc)
    print('Read in successful')


    split = betterSplitter(betterSplitter(betterSplitter(betterSplitter(reads, 0),1),2),3)
    print(len(split))

    with Pool(processes=12) as pool:
        pooledCounts = pool.map(readCounter, split)
    counts =[]
    for i in pooledCounts:
        counts+=i
    writeCSV(counts, outLoc)





