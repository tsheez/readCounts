import subprocess, os
from readCounts import countReads
from tailSeqAnalyzer import tailSeqAnalyzer, seqParser2, csvWriter, splitter
from functools import partial
from multiprocessing import Pool

def queryMaker(readList, outLoc="queryTemp.txt"):
    f= open(outLoc, 'w')
    for i in range(len(readList)):
        f.write(">"+str(i)+"\n")
        f.write(readList[i][0]+"\n")
    f.close()
def subListMaker(nameList, RNASeqList ):
    subList=[]
    for name in nameList:
        for i in range(len(RNASeqList)):
            if name in RNASeqList[i][0]:
                subList.append([RNASeqList[i][0], RNASeqList[i][1]])
                break
    return subList
def subTailSeeker(RNASeqList, countList):
    tails =[]

    for i in range(len(countList)):
        tailtemp=[]
        tailtemp2=[]
        if not countList[i][3]:
            tails.append([countList[i][0],countList[i][2], "No local db blast match", "n/a", "n/a", "n/a" ])
        else:
            for name in countList[i][3]:
                sub = subListMaker([name], RNASeqList)
                tailtemp.append(tailSeqAnalyzer(sub, countList[i]))
            if not tailtemp:
                tails.append([countList[i][0],countList[i][2], "Tail Analyzer Failed", "n/a", "n/a", "n/a" ])
                continue

            for tail in tailtemp:
                if type(tail[3])==int: tailtemp2.append(tail)
            if tailtemp2:
                tailtemp2=sorted(tailtemp2, key=lambda x:abs(x[3]))
                tails.append(tailtemp2[0])
            else: tails.append([countList[i][0],countList[i][2], "Tail Analyzer Failed", "n/a", "n/a", "n/a" ])

        if i%1000 == 0: print (round(i/len(countList)*100),"%") #Progress bar
    return tails

def main(inLoc, outLoc, dbLoc, ranMerLen=13, blastLoc="blastn.exe", processors = 6):
    f=open("blastTemp.txt", 'w')
    f.close()
    f=open("queryTemp.txt", 'w')
    f.close()
    counts = countReads(inLoc, ranMerLen)
    print("Read Counting Successful")
    queryMaker(counts)
    print("Blasting...")
    subprocess.call([blastLoc, '-db', dbLoc , '-query', "queryTemp.txt", '-out', 'blastTemp.txt', '-outfmt', '6'])
    print("Blast successful!")

    f = open("blastTemp.txt", 'r')
    hits = f.readlines()
    f.close()
    for x in counts:
        x.append([])
    for line in hits:
        line = line.split('\t')
        i = int(line[0])
        counts[i][3].append(line[1])

    RNASeqList = seqParser2(dbLoc)

    print("Analyzing tails...")
    func = partial(subTailSeeker, RNASeqList)
    split = splitter(counts, processors)

    with Pool(processes=processors) as pool:
        tailPool = pool.map(func, split )
    tails = []
    for i in tailPool:
        tails+=i

    print("Writing out to CSV")
    csvWriter(tails, outLoc)
    os.remove("queryTemp.txt")
    os.remove("blastTemp.txt")

if __name__== "__main__":
    filepath="C:\\Users\\Lab Admin\\Desktop\\Tim\\"
    main(filepath+"siTOE-F3_S1_L001_R1_001.fastq", filepath+"test.csv", "superset.fa")


