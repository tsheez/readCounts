def analysisParser(inLoc):
    transcriptCount=[]

    file = open(inLoc, 'r')
    file = file.readlines()
    flag = 0
    for line in file:
        if flag:
            line = line.rstrip().split(",")
            transcriptCount.append([line[0],int(line[1])])
        if "Transcript" in line:
            flag = 1
    return transcriptCount
def scriptPercentage(parsed):
    len = 0
    for x in parsed:
        len+=x[1]
    for x in parsed:
        x.append(x[1]/len)
    return parsed











#####################################
inLoc = "C:\\Users\\Tim\\Dropbox\\ProcessedData\\Seq1\\DE_tails_analysis.csv"

parsed = analysisParser(inLoc)
for x in scriptPercentage(parsed):
    print(x)