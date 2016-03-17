import csv

def geneSeparator(key, fileLoc):
    out = []

    f = open(fileLoc, 'r').readlines()

    for i in range(0, len(f)):
        if i == 0: continue
        line = f[i].rstrip()
        line = line.split(',')
        if key in line[2]:
            out.append(line)

    return out

def csvWriter(sepGenes, outLoc):
    f = open(outLoc, 'w')
    f.write("Sequence, Unique Reads, Gene, 3' end, Tail Length, TailSeq\n")
    for i in range(0,len(sepGenes)):
        for j in range(0,(len(sepGenes[i]))):
            f.write(str(sepGenes[i][j]).replace(',', '|'))
            f.write(',')
        f.write('\n')
    f.close()


###########################################
fileLoc = "C:\\Users\\Lab Admin\\Desktop\\TimTemp\\test.csv"
outLoc = "C:\\Users\\Lab Admin\\Desktop\\TimTemp\\sepGenes.csv"
key = "RNU6"

sepGenes = geneSeparator(key, fileLoc)
csvWriter(sepGenes, outLoc)


