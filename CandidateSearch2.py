from scipy import stats
import numpy as np

def tailParser(inLoc):
    f = open(inLoc, 'r')
    tails = f.readlines()
    f.close()

    tailList = []

    for i in range(len(tails)):
        if i==0: continue
        line = tails[i].split(',')
        tailList.append(line)

    return tailList
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches
def findAll(string, sub):
    return(list(find_all(string, sub)))
def getNames(tailList):
    allNames = []
    key = '|-|-|'
    for i in range(len(tailList)):
        lineName = tailList[i][2]
        locs = findAll(lineName, key)
        for j in locs:
            front = lineName.find('|',j+len(key)) + 1
            back = lineName.find('|', front)
            allNames.append(lineName[front:back])
    return list(set(allNames))
def filter(tails):
    filtTails = []
    filtTails2 = []
    for x in tails:
        if x[2] != '[]': filtTails.append(x)
    for x in filtTails:
        if int(x[3])>-10: filtTails2.append(x)
    return filtTails2
def tailStats(tail1, tail2, gene):
    threeLoc1 = []
    threeLoc2 = []
    tailLen1 = []
    tailLen2 = []
    for tail in tail1:
        if gene in tail[2]:
            repeater(int(tail[3]), threeLoc1, int(tail[1]))
            repeater(int(tail[4]), tailLen1, int(tail[1]))
            #threeLoc1.append(int(tail[3]))
            #tailLen1.append(int(tail[4]))
    for tail in tail2:
        if gene in tail[2]:
            repeater(int(tail[3]), threeLoc2, int(tail[1]))
            repeater(int(tail[4]), tailLen2, int(tail[1]))
            #threeLoc2.append(int(tail[3]))
            #tailLen2.append(int(tail[4]))
    if not threeLoc1 or not threeLoc2:
        pLoc = "nan"
        pTail = "nan"
    else:
        pLoc = stats.ttest_ind(threeLoc1, threeLoc2)[1]
        pTail = stats.ttest_ind(tailLen1, tailLen2)[1]
    return gene, len(threeLoc1), np.average(threeLoc1), np.average(tailLen1), len(threeLoc2), np.average(threeLoc2), np.average(tailLen2), pLoc, pTail
def CSVWriter(stats, outLoc):
    f = open(outLoc, 'w')
    f.write('Gene'+','+'1-Reads'+','+'1-3\' Loc Avg'+','+'1-Tail Length Avg'+','+'2-Reads'+','+'2-3\' Loc Avg'+','+'2-Tail Length Avg'+','+'P-Value 3\' Loc'+','+'P-Value Tail Length\n')
    for line in stats:
        for item in line:
            f.write(str(item)+',')
        f.write('\n')
def repeater(item, list, reps):
    for i in range(reps):
        list.append(item)
    return



#######################################

tailLoc1 = "C:\\Users\\Tim\\Desktop\\siLuc_no5S_tails.csv"
tailLoc2 = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\siTOE_No5S_tails.csv"
outLoc = "C:\\Users\\Tim\\Desktop\\temp.csv"

tails1 = tailParser(tailLoc1)
tails2 = tailParser(tailLoc2)
tails1 = filter(tails1)
tails2 = filter(tails2)
genes1 = getNames(tails1)
genes2 = getNames(tails2)

compiledGenes = set(genes1+genes2)

data = []
for gene in compiledGenes:
    data.append(tailStats(tails1, tails2, gene))

print(data)
CSVWriter(data, outLoc)
