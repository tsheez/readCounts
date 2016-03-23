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
    return set(allNames)
def readPercent(tails, key):
    count = 0
    totReads = 0
    for i in range(len(tails)):
        if i == 0: continue
        totReads += int(tails[i][1])
    for i in range(len(tails)):
        if key in tails[i][2]:
            count += int(tails[i][1])
    return count/totReads
def filter(tails):
    filtTails = []
    filtTails2 = []
    for x in tails:
        if x[2] != '[]': filtTails.append(x)
    for x in filtTails:
        if int(x[3])>-10: filtTails2.append(x)
    return filtTails2
def analyzer(tails, name):
    count = 0
    endTot = 0


######################
inLoc = "C:\\Users\\Tim\\Desktop\\siTOE_No5S_tails.csv"

tails = tailParser(inLoc)

typeList = ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA',\
            'rRNA', 'scRNA', 'snRNA', 'snoRNA',\
            'ribozyme', 'sRNA', 'scaRNA', '[]']
tails = filter(tails)
names = getNames(tails)
for x in typeList:
    print(x,'\t', readPercent(tails, x))

for name in names:
    count = 0
    endTot = 0
    tailLen = 0
    for tail in tails:
        if name in tail[2]:
            numReads = int(tail[1])
            endTot += (int(tail[3])*numReads)
            tailLen += (int(tail[4])*numReads)
            count += numReads
    print(name, count, endTot/count, tailLen/count)
