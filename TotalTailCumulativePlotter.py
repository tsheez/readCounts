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
def geneSeparator(gene, tails):
    out=[]
    for tail in tails:
        if gene in tail[2]:
            out.append(tail)
    return out
def totalTailGrabber(tails):
    totalTails = []
    for tail in tails:
        totalTails.append((tail[2], int(tail[1]),int(tail[3])+int(tail[4])))
    return totalTails
def cumulativePlotter(totalTails):
    out=[]
    totalTails = sorted(totalTails, key = lambda x:x[1])
    len = 0
    for tail in totalTails:
        len+=tail[1]
    for i in range(-10, 20):
        count = 0
        for tail in totalTails:
            if tail[2]<= i: count+=tail[1]
        out.append((i,count/len))
    return out



def main(inLoc, gene):
    tails = tailFilter(tailParser(inLoc))
    tails = geneSeparator(gene, tails)
    tails = totalTailGrabber(tails)
    plot = cumulativePlotter(tails)
    for line in plot:
        print(line)

####################################
inLoc = "C:\\Users\\Tim\\Dropbox\\ProcessedData\\Seq1\\siTOE_tails.csv"
gene = "U1"
main(inLoc, gene)

