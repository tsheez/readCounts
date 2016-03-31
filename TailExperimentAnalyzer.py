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
def getNames(tails):
    names = []
    for tail in tails:
        if tail[2][0]=='[':
            names.append(tail[2][tail[2].find('[\'>')+3:tail[2].find('|')])
    return list(set(names))
def tailFilter(tails):
    tail1=[]
    tail2=[]
    for tail in tails:
        if 'n/a' not in tail[3]:
            tail1.append(tail)
    for tail in tail1:
        if int(tail[3])>=-10:
            tail2.append(tail)
    return tail2
def typeCounter(tails):
    typeCounts=[]
    typeList = ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA',\
            'rRNA', 'scRNA', 'snRNA', 'snoRNA',\
            'ribozyme', 'sRNA', 'scaRNA', 'lncRNA', \
            'No local db blast match', 'Tail Analyzer Failed']
    for item in typeList:
        count = 0
        for tail in tails:
            if item in tail[2]:
                count += int(tail[1])
        typeCounts.append((item, count))
    return typeCounts







##################################
filepath = "C:\\Users\\Lab Admin\\Desktop\\Tim\\"
inLoc = filepath + "test.csv"

tails = tailParser(inLoc)
print(len(tails))
print(typeCounter(tails))
tails = tailFilter(tails)
print(len(tails))
print(typeCounter(tails))
print(tails[5])