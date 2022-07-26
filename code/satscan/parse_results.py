import numpy as np
import re



#open and read input file
inTxt = open(r'C:\Users\ahohl\Google Drive\Research Projects\greenChicago\data\satscan\results_v2.txt', 'r')
inArr = inTxt.readlines()

#open output file
outFile = open(r'C:\Users\ahohl\Google Drive\Research Projects\greenChicago\data\satscan\results_v2_arr.txt', 'w')


#create array that holds the first line position of each cluster
firstLine = []
i = 0
while i < len(inArr):
    try:
        if inArr[i].split('.')[1] == 'Location IDs included':
            firstLine.append(i)
    except:
        pass
    i += 1

#append position of last line
firstLine.append(len(inArr))

#itearte through inArr using line positions from firstLine
i = 0
#for each break point
while i < len(firstLine)-1:
    fr = firstLine[i]
    to = firstLine[i+1]

    j = fr
    # print(inArr[fr])

    # x-coord, y-coord, radius, Number of cases, Mean inside, Mean outside, Variance, StDev, LLR, p-value
    clu = [0,0,0,0,0,0,0,0,0,0,0,0,0]


    clu[0] = inArr[fr].split('.')[0] #cluster ID

    #for each line between break points
    while j < to:
        #print(inArr[j])
        if inArr[j].split(':')[0].strip() == 'Coordinates / radius..':
            #print(inArr[j])
            clu[0] = re.sub('[()]', '', inArr[j].split(': ')[1].strip().split('/')[0].split(',')[0])    #x-coord
            clu[1] = re.sub('[()]', '', inArr[j].split(': ')[1].strip().split('/')[0].split(',')[1])    #y-coord
            clu[2] = inArr[j].split(': ')[1].strip().split('/')[1].strip()                              # radius
            clu[3] = inArr[j + 1].split(': ')[1].strip()                                                #Number of cases
            clu[4] = inArr[j + 2].split(': ')[1].strip()                                                #Mean inside
            clu[5] = inArr[j + 3].split(': ')[1].strip()                                                #Mean outside
            clu[6] = inArr[j + 4].split(': ')[1].strip()                                                #Variance
            clu[7] = inArr[j + 5].split(': ')[1].strip()                                                #StDev
            clu[8] = inArr[j + 6].split(': ')[1].strip()                                                #LLR
            clu[9] = inArr[j + 7].split(': ')[1].strip()                                                #p-value

            print(clu[1])
            outFile.write(clu[0] + ',' +
                          clu[1] + ',' +
                          clu[2] + ',' +
                          clu[3] + ',' +
                          clu[4] + ',' +
                          clu[5] + ',' +
                          clu[6] + ',' +
                          clu[7] + ',' +
                          clu[8] + ',' +
                          clu[9] + '\n')

        j += 1
    #print('break')

    i += 1


outFile.close()
