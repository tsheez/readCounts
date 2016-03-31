def tailParser(inLoc):
    f = open(inLoc, 'r')
    tails = f.readlines()
    f.close()

    tailList = []

    for i in range(len(tails)):
        if i==0: continue
        line = tails[i].split(',')
        tailList.append(line)






##################################
filepath = "C:\\Users\\Lab Admin\\Desktop\\Tim\\"
inLoc = filepath + "test.csv"

tails = tailParser(inLoc)