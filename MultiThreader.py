def seqParser2(seqLoc):
    f=open(seqLoc,'r')
    RNAseqs = f.readlines()
    f.close()
    RNAseqlist = []
    for i in range(0, len(RNAseqs)):
        if RNAseqs[i][0] == ">":
            RNAseqlist.append([RNAseqs[i].rstrip(),RNAseqs[i+1].rstrip()])
    return RNAseqlist

seqlist = seqParser2("C:\\Users\\Tim\\Dropbox\\Data\\Resources\\FASTA_Subsets\\all_small_RNA.fa")