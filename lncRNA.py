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
    inLoc1="C:\\Users\\Tim\\Dropbox\Data\\Resources\\\FASTA_Subsets\\all_small_RNA.fa"
    inLoc2="C:\\Users\\Tim\\Dropbox\Data\\Resources\\\FASTA_Subsets\\gencode.v24.lncRNA_transcripts.fa"
    outLoc = "C:\\Users\\Tim\\Dropbox\Data\\Resources\\\FASTA_Subsets\\superset.fa"


    in1 = open(inLoc1, 'r')
    in1 = in1.readlines()
    print(len(in1))
    out = open(outLoc, 'w')
    for line in in1:
        if line =="--\n":
            continue
        elif line[0]==">":
            line=line.rstrip()
            front = line.find("|-|-|")+ len("|-|-|")
            back = line.find("|", front)
            out.write(">"+line[front:back]+"|"+line[line[:-1].rfind("|")+1:line.rfind("|")]+"\n")
        else:
            out.write(line)

    in2 = open(inLoc2, 'r').readlines()
    print(len(in2))
    for line in in2:
        if line[0]==">":
            x = line.rfind("OTTHUMT")
            y = line.find("|", x)+1
            z = line.find("|", y)
            out.write(">"+line[y:z]+"|lncRNA\n")
        else:
            out.write(line)

    out = open(outLoc, 'r').readlines()
    print(len(out))
    for x in range(20):
        print(out[x])











