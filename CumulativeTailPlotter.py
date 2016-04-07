import csv, os

def threePosMaker(inLoc):
    limit = -10 #sets a limit to the accepted truncation of RNAs to include in the analysis
    inLoc = inLoc



    fcsv = open(inLoc)
    Masterlist = []
    for line in csv.reader(fcsv):
        Masterlist.append(line)
    fcsv.close()

    #Sums up total number of mappable reads from the Masterlist
    totalreads = 0
    totalfulllengthreads = 0
    totaltailedreads = 0

    for line in Masterlist[1:]:
        if line[2][0] == '[' and int(line[3]) >= limit:
            totalreads += int(line[1])
            if int(line[3]) > -1:
                totalfulllengthreads += int(line[1])
            if int(line[4]) > 0:
                totaltailedreads += int(line[1])



    #Generates the Positionlist and Taillists
    Endpositionlist = [['Position','%Reads','%Tailed (of total)','%Tail>=2','%Tail>=3','%Tail>=4','%Tail>=5',
        '%Tail>=6','%Tail>=7','%Tail>=8','%Tail>=9','%Tail>=10','Cumulative 3-end','Cumulative tailposition', 				'%Tailed (for each position)','Cumulative total length']]

    l = -200
    while l < 51:
        Endpositionlist.append([l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        l += 1

    Taillengthlist = [['Tail length','%Encoded tails','%Unencoded tails']]
    Tailcomposition = [['Tail Position','A','T','C','G']]

    l = 0
    while l < 26:
        Taillengthlist.append([l, 0, 0])
        Tailcomposition.append([l+1, 0, 0, 0, 0])
        l += 1


    for line in Masterlist[1:]:
        if line[2][0] == '[' and int(line[3]) >= limit: #Requires an identified gene

            Endpositionlist[int(line[3])+201][1] += 100*float(line[1])/totalreads

            x = int(line[3])
            while x < 51: #Generates the cumulative 3'end list
                Endpositionlist[x + 201][12] += 100*float(line[1])/totalreads
                x += 1

            x = int(line[3]) + int(line[4])
            while x < 51: #Generates the cumulative total length list
                Endpositionlist[x + 201][15] += 100*float(line[1])/totalreads
                x += 1

            if int(line[4]) < 26:
                Taillengthlist[int(line[4])+1][2] += 100*float(line[1])/totalreads

            if int(line[4]) > 0:
                x = 0
                while x < int(line[4]) and x < 10: #Calculates tail lengths
                    Endpositionlist[int(line[3])+201][x+2] += 100*float(line[1])/totalreads
                    x += 1

                x = int(line[3])
                while x < 51: #Generates the cumulative tailposition list
                    Endpositionlist[x + 201][13] += 100*float(line[1])/totaltailedreads
                    x += 1

                x = 0
                for n in line[5]: #Calculates average tail composition
                    x += 1
                    if x < 27:
                        if n == 'A':
                            Tailcomposition[x][1] += 100*float(line[1])/totalreads
                        if n == 'T':
                            Tailcomposition[x][2] += 100*float(line[1])/totalreads
                        if n == 'C':
                            Tailcomposition[x][3] += 100*float(line[1])/totalreads
                        if n == 'G':
                            Tailcomposition[x][4] += 100*float(line[1])/totalreads

            if int(line[3]) > -1:
                if int(line[3]) < 26:
                    Taillengthlist[int(line[3])+1][1] += 100*float(line[1])/totalfulllengthreads
                else:
                    Taillengthlist[26][1] += 100*float(line[1])/totalfulllengthreads


    l = 0
    for line in Endpositionlist[1:]:
        l += 1
        if float(line[1]) >= 0.1:
            Endpositionlist[l][14] += 100*float(line[2])/line[1]

    '''
    savefile = inLoc.replace("tails", "Positiondata")
    f = open(savefile, 'w', newline='')
    csv.writer(f).writerows(Endpositionlist)
    f.close()

    savefile = inLoc.replace("tails", "Taillengthdata")
    f = open(savefile, 'w', newline='')
    csv.writer(f).writerows(Taillengthlist)
    f.close()

    savefile =inLoc.replace("tails","Tailcomposition")
    f = open(savefile, 'w', newline='')
    csv.writer(f).writerows(Tailcomposition)
    f.close()

    print (inLoc.split('/')[-1] + ": n= " + str(totalreads))
    '''

    return Endpositionlist
def tailLenMaker(inLoc):
    limit = -10 #sets a limit to the accepted truncation of RNAs to include in the analysis
    inLoc = inLoc



    fcsv = open(inLoc)
    Masterlist = []
    for line in csv.reader(fcsv):
        Masterlist.append(line)
    fcsv.close()

    #Sums up total number of mappable reads from the Masterlist
    totalreads = 0
    totalfulllengthreads = 0
    totaltailedreads = 0

    for line in Masterlist[1:]:
        if line[2][0] == '[' and int(line[3]) >= limit:
            totalreads += int(line[1])
            if int(line[3]) > -1:
                totalfulllengthreads += int(line[1])
            if int(line[4]) > 0:
                totaltailedreads += int(line[1])



    #Generates the Positionlist and Taillists
    Endpositionlist = [['Position','%Reads','%Tailed (of total)','%Tail>=2','%Tail>=3','%Tail>=4','%Tail>=5',
        '%Tail>=6','%Tail>=7','%Tail>=8','%Tail>=9','%Tail>=10','Cumulative 3-end','Cumulative tailposition', 				'%Tailed (for each position)','Cumulative total length']]

    l = -200
    while l < 51:
        Endpositionlist.append([l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        l += 1

    Taillengthlist = [['Tail length','%Encoded tails','%Unencoded tails']]
    Tailcomposition = [['Tail Position','A','T','C','G']]

    l = 0
    while l < 26:
        Taillengthlist.append([l, 0, 0])
        Tailcomposition.append([l+1, 0, 0, 0, 0])
        l += 1


    for line in Masterlist[1:]:
        if line[2][0] == '[' and int(line[3]) >= limit: #Requires an identified gene

            Endpositionlist[int(line[3])+201][1] += 100*float(line[1])/totalreads

            x = int(line[3])
            while x < 51: #Generates the cumulative 3'end list
                Endpositionlist[x + 201][12] += 100*float(line[1])/totalreads
                x += 1

            x = int(line[3]) + int(line[4])
            while x < 51: #Generates the cumulative total length list
                Endpositionlist[x + 201][15] += 100*float(line[1])/totalreads
                x += 1

            if int(line[4]) < 26:
                Taillengthlist[int(line[4])+1][2] += 100*float(line[1])/totalreads

            if int(line[4]) > 0:
                x = 0
                while x < int(line[4]) and x < 10: #Calculates tail lengths
                    Endpositionlist[int(line[3])+201][x+2] += 100*float(line[1])/totalreads
                    x += 1

                x = int(line[3])
                while x < 51: #Generates the cumulative tailposition list
                    Endpositionlist[x + 201][13] += 100*float(line[1])/totaltailedreads
                    x += 1

                x = 0
                for n in line[5]: #Calculates average tail composition
                    x += 1
                    if x < 27:
                        if n == 'A':
                            Tailcomposition[x][1] += 100*float(line[1])/totalreads
                        if n == 'T':
                            Tailcomposition[x][2] += 100*float(line[1])/totalreads
                        if n == 'C':
                            Tailcomposition[x][3] += 100*float(line[1])/totalreads
                        if n == 'G':
                            Tailcomposition[x][4] += 100*float(line[1])/totalreads

            if int(line[3]) > -1:
                if int(line[3]) < 26:
                    Taillengthlist[int(line[3])+1][1] += 100*float(line[1])/totalfulllengthreads
                else:
                    Taillengthlist[26][1] += 100*float(line[1])/totalfulllengthreads


    l = 0
    for line in Endpositionlist[1:]:
        l += 1
        if float(line[1]) >= 0.1:
            Endpositionlist[l][14] += 100*float(line[2])/line[1]

    '''
    savefile = inLoc.replace("tails", "Positiondata")
    f = open(savefile, 'w', newline='')
    csv.writer(f).writerows(Endpositionlist)
    f.close()

    savefile = inLoc.replace("tails", "Taillengthdata")
    f = open(savefile, 'w', newline='')
    csv.writer(f).writerows(Taillengthlist)
    f.close()

    savefile =inLoc.replace("tails","Tailcomposition")
    f = open(savefile, 'w', newline='')
    csv.writer(f).writerows(Tailcomposition)
    f.close()

    print (inLoc.split('/')[-1] + ": n= " + str(totalreads))
    '''

    return Taillengthlist
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

def geneSepWriter(sepGenes, outLoc):
    f = open(outLoc, 'w')
    f.write("Sequence, Unique Reads, Gene, 3' end, Tail Length, TailSeq\n")
    for i in range(0,len(sepGenes)):
        for j in range(0,(len(sepGenes[i]))):
            f.write(str(sepGenes[i][j]).replace(',', '|'))
            f.write(',')
        f.write('\n')
    f.close()

def cumulativeTailPlotter(inLoc, geneList):
    master3Pos=[]
    masterTails = []
    for gene in geneList:
        geneSeparator(gene, inLoc)
        geneSepWriter(geneSeparator(gene, inLoc), "SepGeneTemp.txt")

        temps = threePosMaker("SepGeneTemp.txt")
        threepos = []
        for temp in temps:
            threepos.append([temp[0], temp[12]])
        if not master3Pos:
            column = []
            for item in threepos:
                column.append(item[0])
            master3Pos.append(column)
        threepos[0][1] = gene
        column = []
        for item in threepos:
            column.append(item[1])
        master3Pos.append(column)

        temps = tailLenMaker("SepGeneTemp.txt")
        tailLen = []
        for temp in temps:
            tailLen.append([temp[0], temp[2]])
        if not masterTails:
            column=[]
            for item in tailLen:
                column.append(item[0])
            masterTails.append(column)
        tailLen[0][1] = gene
        column = []
        for item in tailLen:
            column.append(item[1])
        masterTails.append(column)



        os.remove("SepGeneTemp.txt")



##############################
inLoc = 'C:\\Users\\Tim\\Dropbox\\Data\\TLS005\\TLS006 Analysis\\siLuc_tails.csv'
geneList = ['RNU6','RNU4']

cumulativeTailPlotter(inLoc, geneList)