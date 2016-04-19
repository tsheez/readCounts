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




def main(inLoc1, inLoc2):
    out = []
    Tot = scriptPercentage(analysisParser(inLoc1))
    RIP = scriptPercentage(analysisParser(inLoc2))

    for x in RIP:
        flag = 0
        for y in Tot:
            if x[0] == y[0]:
                x.append(y[2])
                flag = 1
                break
        if not flag:
            x.append('nan')
    for x in RIP:
        if x[3] == 'nan': pass
        else:
            out.append((x[0],(x[2]-x[3])/x[3]))
    for x in out:
        print(x)









#####################################
inLoc1 = "C:\\Users\\Tim\\Dropbox\\ProcessedData\\Seq1\\WT_tails_analysis.csv"
inLoc2 = "C:\\Users\\Tim\\Dropbox\\ProcessedData\\041516_RIPS\\WT_RIP_Tails_Analysis.csv"


main(inLoc1, inLoc2)
