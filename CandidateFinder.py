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
        #pLoc = stats.ttest_ind(threeLoc1, threeLoc2)[1]
        #pTail = stats.ttest_ind(tailLen1, tailLen2)[1]
        pLoc = stats.ks_2samp(threeLoc1, threeLoc2)[1]
        pTail = stats.ks_2samp(tailLen1, tailLen2)[1]
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

def main (inLoc1, inLoc2, outLoc):

    tails1 = tailFilter(tailParser(inLoc1))
    tails2 = tailFilter(tailParser(inLoc2))
    compiledGenes = set(getNames(tails1) + getNames(tails2))

    data = []
    for gene in compiledGenes:
        data.append(tailStats(tails1, tails2, gene))

    CSVWriter(data, outLoc)


############################
if __name__=="__main__":
    inLoc1 = "C:\\Users\\Lab Admin\\Desktop\\AnalysisStrategy2\\siLuc_total_tails2.csv"
    inLoc2 = "C:\\Users\\Lab Admin\\Desktop\\AnalysisStrategy2\\DE_total_tails2.csv"
    outLoc = "C:\\Users\\Lab Admin\\Desktop\\AnalysisStrategy2\\Candidates\\siLuc_DE_KS.csv"

    main(inLoc1, inLoc2, outLoc)


