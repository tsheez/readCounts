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
    typeCounts = sorted(typeCounts, key=lambda x:x[1], reverse=True)
    return typeCounts
def transcriptCounter(tails):
    names = getNames(tails)
    transcripts=[]
    for name in names:
        count = 0
        for tail in tails:
            if name in tail[2]:
                count+=int(tail[1])
        transcripts.append((name, count))
    transcripts=sorted(transcripts, key=lambda x:x[1], reverse=True)
    return transcripts

def main(inLoc, outLoc):
    tails = tailParser(inLoc)
    out = open(outLoc, 'w')

    out.write(inLoc+"\n\nUnfiltered\nType,Reads\n")
    types = typeCounter(tails)
    for item in types:
        out.write(item[0]+','+str(item[1])+'\n')

    tails = tailFilter(tails)
    out.write("\n\nfiltered\nType,Reads\n")
    types = typeCounter(tails)
    for item in types:
        out.write(item[0]+','+str(item[1])+'\n')

    transcripts = transcriptCounter(tails)
    out.write("\n\nTranscript,Reads\n")
    for item in transcripts:
        out.write(item[0]+','+str(item[1])+'\n')
    out.close()









##################################

inLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\AnalysisStrategy2\\WT_total_tails2.csv"
outLoc = "C:\\Users\\Tim\\Dropbox\\Data\\TLS004\\AnalysisStrategy2\\WT_total_tails2_analysis.csv"

main(inLoc, outLoc)
