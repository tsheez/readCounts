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
        if not countList[i][3]:
            tails.append([countList[i][0],countList[i][2], "No local db blast match", "n/a", "n/a", "n/a" ])
        else:
            sublist = subListMaker(countList[i][3], RNASeqList)
            tails.append(tailSeqAnalyzer(sublist, countList[i]))
        if i%10000 == 0: print (round(i/len(countList)*100),"%") #Progress bar
    return tails

def main(inLoc, outLoc, dbLoc, ranMerLen=12, blastLoc="blastn.exe", processors = 12):
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
    main("C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\2016-03-09-MiSeqRaw\\C-33927852\\Data\\Intensities\\BaseCalls\\C_S3_L001_R1_001.fastq",\
         "C:\\Users\\Tim\\Desktop\\test.csv", "all_small_RNA_trim.fa")

    '''
    inLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\2016-03-09-MiSeqRaw\\C-33927852\\Data\\Intensities\\BaseCalls\\C_S3_L001_R1_001.fastq"
    outLoc = "C:\\Users\\Tim\\Desktop\\test.csv"
    dbLoc = "all_small_RNA_trim.fa"

    counts = countReads(inLoc, ranMer=12)

    print("count successful", len(counts), "reads")

    queryMaker(counts)

    print("query created")

    subprocess.call(["blastn.exe", '-db', dbLoc , '-query', "temp.txt", '-out', 'temp2.txt', '-outfmt', '6'])

    print("blast finished")

    f = open("temp2.txt", 'r')
    hits = f.readlines()
    f.close()
    for x in counts:
        x.append([])
    for line in hits:
        line = line.split('\t')
        i = int(line[0])
        counts[i][3].append(line[1])

    print("counts modified")

    RNASeqList = seqParser2("all_small_RNA_trim.fa")

    print("seqlist created")

    tails = []
    for i in range(len(counts)):
        try:
            if not counts[i][3]:
                tails.append([counts[i][0],counts[i][2], "No local db blast match", "n/a", "n/a", "n/a" ])
            else:
                sublist = subListMaker(counts[i][3], RNASeqList)
                tails.append(tailSeqAnalyzer(sublist, counts[i]))
        except IndexError:
            print(counts[i])
        if i%1000 == 0: print (round(i/len(counts)*100),"%")


    csvWriter(tails, outLoc)
    '''

