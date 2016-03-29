import subprocess
from readCounts import countReads
from tailSeqAnalyzer import tailSeqAnalyzer, seqParser2, csvWriter

def queryMaker(readList):
    f= open("temp.txt", 'w')
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




if __name__== "__main__":
    inLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\2016-03-09-MiSeqRaw\\C-33927852\\Data\\Intensities\\BaseCalls\\C_S3_L001_R1_001.fastq"
    outLoc = "C:\\Users\\Tim\\Desktop\\test.csv"

    counts = countReads(inLoc)

    print("count successful")

    queryMaker(counts)

    print("query created")

    subprocess.call(["blastn.exe", '-db', "all_small_RNA_trim.fa", '-query', "temp.txt", '-out', 'temp2.txt', '-outfmt', '6'])

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
        if not counts[i][3]:
            tails.append([counts[i][0],counts[i][2], "No local db blast match", "n/a", "n/a", "n/a" ])
        else:
            sublist = subListMaker(counts[i][3], RNASeqList)
            tails.append(tailSeqAnalyzer(sublist, counts[i]))
        if i%1000 == 0: print (round(i/len(counts)*100),"%")


    csvWriter(tails, outLoc)


