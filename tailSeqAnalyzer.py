from multiprocessing import Pool
from functools import partial


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

        RNAseqlist.append([RNAseqs[a+1:b].replace(",",""),RNAseqs[c+1:d].replace("\n","")])

        x = RNAseqs.find('>', c)





    return RNAseqlist
#output = nested list [[name, seq],...]
def seqParser2(seqLoc):
    f=open(seqLoc,'r')
    RNAseqs = f.readlines()
    f.close()
    RNAseqlist = []
    for i in range(0, len(RNAseqs)):
        if RNAseqs[i][0] == ">":
            RNAseqlist.append([RNAseqs[i].rstrip(),RNAseqs[i+1].rstrip()])
    return RNAseqlist
def countParser(countLoc):
    f=open(countLoc, 'r')
    counts = f.readlines()
    f.close()

    countList=[]
    for i in range(0, len(counts)):
        if i == 0: continue
        counts[i] = counts[i].split(',')
        counts[i][1] = int(counts[i][1].rstrip())
        counts[i][2] = int(counts[i][2].rstrip())
        countList.append(counts[i])


    return countList
#output = nested list [[read, total counts, unique counts],...]
def tailSeqAnalyzer(RNAseqlist, read):
    #Sequence, UniqueReads, Gene, 3'end, tailLength, TailSeq
    #read = [sequence, total reads, unique read]

    seq = read[0]
    uniqueReads = read[2]

    #Identifies the best matching gene(s) for the current read
    linenumber = -1
    topgenescore = 0
    matchgeneseq = []
    matchgenelist = []
    matchlinenumber = []
    matchposition = []
    matchlength = []
    for geneline in RNAseqlist:

        linenumber += 1
        geneseq = geneline[1]
        x = 0
        genescore = 0

        while x <= len(geneseq) - 10:
            y = 0
            m = 0
            while y < (len(seq)) and (x + y) < len(geneseq) and m == 0:
                if seq[y] == geneseq[x+y]:
                    y += 1
                else:
                    m += 1

            if y > genescore:
                genescore = y
                geneposition = x
            x += 1

        if genescore == topgenescore:
            matchgeneseq.append(geneseq)
            matchgenelist.append(geneline[0])
            matchlinenumber.append(linenumber)
            matchposition.append(geneposition)
            matchlength.append(genescore)
        if genescore > topgenescore:
            matchgeneseq = [geneseq]
            matchgenelist = [geneline[0]]
            matchlinenumber = [linenumber]
            matchposition = [geneposition]
            matchlength = [genescore]
            topgenescore = genescore



    #Identifies the 3'end, tail length and tail sequence

    g = -1
    toss = 0
    errmsg = ''

    if matchlength[0] < 10: #Requires 10 or more matches
        toss = 1
        errmsg = 'ERR(Short)'

    else:

        for geneseq in matchgeneseq:
            g += 1

            if matchlength[g] < len(seq) and toss == 0:

                x = 1
                match = 0
                while (x + matchlength[g]) < len(seq):

                    if (matchposition[g] + matchlength[g] + x) < len(geneseq):
                        if seq[matchlength[g] + x] == geneseq[matchposition[g] +
                                    matchlength[g] + x]:
                            match += 1
                    x += 1

                    if x > 3 and float(match)/(x-1) >= 0.75: #Tosses the read if it matches the gene 75% at any point 3n after the mismatch
                        toss = 1
                        errmsg = 'ERR(Mismatch)'

                if x == 3 and match == 2: #Tosses if last 2n of seq match gene after mismatch
                        toss = 1
                        errmsg = 'ERR(Mismatch)'

                if toss == 0 and x > 1: #Looks for deletion or insert if the tail is longer than 1n
                    y = 1
                    match = 0
                    while (y + matchlength[g] + 1) < len(seq):

                        if (matchposition[g] + matchlength[g] + y) < len(geneseq):
                            if seq[matchlength[g] + y + 1] == geneseq[matchposition[g]
                                    + matchlength[g] + y]:
                                match += 1
                        y += 1

                        if y > 3 and float(match)/(y-1) >= 0.75: #Tosses the read if it matches the gene 75% at any point 3n after the mismatch
                            toss = 1
                            errmsg = 'ERR(Insert)'


                    y = 1
                    match = 0
                    while (y + matchlength[g]) < len(seq) and (matchposition[g] +
                                matchlength[g] + y + 1) < len(geneseq):
                        if matchposition[g] + matchlength[g] + y + 1 < len(geneseq):
                            if seq[matchlength[g] + y] == geneseq[matchposition[g] +
                                        matchlength[g] + y + 1]:
                                match += 1
                        y += 1

                        if y > 3 and float(match)/(y-1) >= 0.75: #Tosses the read if it matches the gene 75% at any point 3n after the mismatch
                            toss = 1
                            errmsg = 'ERR(Deletion)'


                if toss == 0 and x > matchlength[0]: #Tosses if tail is longer than sequenced part of gene
                    toss = 1
                    errmsg = 'ERR(Tail longer than body)'

            if matchlength[g] == len(seq):
                x = 0



    if toss == 0:
        end = matchposition[0] + matchlength[0] - len(geneseq)
        taillength = x
        tailseq = seq[matchlength[0]:]

    else:
        matchgenelist = []
        end = errmsg
        taillength = ''
        tailseq = ''

    return [seq, uniqueReads, matchgenelist, end, taillength, tailseq]
def csvWriter(tails, outLoc):
    tails = sorted(tails, key=lambda x:x[1], reverse=True)
    f = open(outLoc, 'w')
    f.write("Sequence, Unique Reads, Gene, 3' end, Tail Length, TailSeq\n")
    for i in range(0,len(tails)):
        for j in range(0,(len(tails[i]))):
            f.write(str(tails[i][j]).replace(',', '|'))
            f.write(',')
        f.write('\n')
    f.close()
def splitter(counts, numSplits):
    split =[]
    for i in range(0,numSplits):
        split.append([])
    for i in range(0,len(counts)):
        split[i%numSplits].append(counts[i])
    return split
def tailSeeker(RNAseqlist, readList):
    tails =[]
    for i in range(len(readList)):
        tails.append(tailSeqAnalyzer(RNAseqlist, readList[i]))
        if i%10 == 0: print (round(i/len(readList)*100),"%")  # Progress bar
    return tails


#####test code



seqLoc = "C:\\Users\\Tim\\Dropbox\\Data\\Resources\\FASTA_Subsets\\gencode.v24.lncRNA_transcripts.fa"
countLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\2016-03-14-MiSeqAnalysis\\siTOE_no5S_Counts.csv"
outLoc = "C:\\Users\\Tim\\Desktop\\test.csv"
processors = 8

'''
seqLoc = sys.argv[1]
countLoc = sys.argv[2]
outLoc = sys.argv[3]
processors = int(sys.argv[4])
'''


if __name__=='__main__':
    RNAseqlist = seqParser2(seqLoc)
    print(len(RNAseqlist))
    counts = countParser(countLoc)
    print('read in successful')


    func = partial(tailSeeker, RNAseqlist)
    split = splitter(counts, processors)

    with Pool(processes=processors) as pool:
        tailPool = pool.map(func, split )

    tails = []
    for i in tailPool:
        tails+=i

    csvWriter(tails, outLoc)
