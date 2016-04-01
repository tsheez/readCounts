def lncRNAParser(seqLoc):
    f=open(seqLoc,'r')
    RNAseqs = f.readlines()
    f.close()
    RNAseqlist = []
    seqFiltered = []
    seqFiltered2 = []
    for i in range(0, len(RNAseqs)):
        if RNAseqs[i][0] == ">":
            RNAseqlist.append([RNAseqs[i].rstrip(),RNAseqs[i+1].rstrip()])

    for seq in RNAseqlist:
        if len(seq[1])<=400:
            seqFiltered.append(seq)
    for seq in seqFiltered:
        if len(seq[1])>=100:
            seqFiltered2.append(seq)


    return seqFiltered2





########################
if __name__=="__main__":
    f = open('temp.txt', 'r').readlines()
    out = open("temp2.txt", 'w')
    for line in f:
        line =line.rstrip()
        out.write(line)
    out.close()








