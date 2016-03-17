def seqParser(seqLoc):

    f = open(seqLoc, 'r')
    RNAseqs = f.read()
    f.close()

    RNAseqlist = []
    x = 1
    for line in RNAseqs:
        a = RNAseqs.find('>', x-1)
        b = RNAseqs.find('\n', a)
        c = RNAseqs.find('\n', b)
        d = RNAseqs.find('>', c)
        if d == -1:
            d = len(RNAseqs)
        while RNAseqs[d-1] == '\n':
            d -= 1

        RNAseqlist.append([RNAseqs[a+1:b],RNAseqs[c+1:d].replace("\n","")])

        x = RNAseqs.find('>', c)





    print(RNAseqlist[1])
#output = nested list [[name, seq],...]











#####test code
seqLoc = "C:\\Users\\Lab Admin\\Desktop\\TimTemp\\Human small RNA sequences.txt"
seqParser(seqLoc)
