import subprocess
from readCounts import countReads

def queryMaker(readList):
    f= open("temp.txt", 'w')
    for i in range(len(readList)):
        f.write(">"+str(i)+"\n")
        f.write(readList[i][0]+"\n")
    f.close()




def blaster(query, database, out):

    subprocess.call(["blastn.exe", '-db', database, '-query', query, '-out', 'temp2.txt', '-outfmt', '6'])


if __name__== "__main__":
    inLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\2016-03-09-MiSeqRaw\\C-33927852\\Data\\Intensities\\BaseCalls\\C_S3_L001_R1_001.fastq"
    counts = countReads(inLoc)

    queryMaker(counts)

    subprocess.call(["blastn.exe", '-db', "all_small_RNA_trim.fa", '-query', "temp.txt", '-out', 'temp2.txt', '-outfmt', '6'])
    f = open("temp2.txt", 'r')
    hits = f.readlines()
    f.close()
    for x in counts:
        x.append([])

    for line in hits:
        line = line.split('\t')
        i = int(line[0])
        counts[i][3].append(line[1])
