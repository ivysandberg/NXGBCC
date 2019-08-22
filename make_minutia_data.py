'''
Create fake minutia and fingerprint data to test 

'''

import numpy as np

# make some example data to sample workflow

num_prints = 1000
min_minutia = 25
max_minutia = 125
xlim = (0,100)
ylim = (0,100)

def makeMinutia(xlim, ylim):
    x = np.random.uniform(xlim[0], xlim[1])
    y = np.random.uniform(ylim[0], ylim[1])
    rad = np.random.uniform(0, 2*math.pi)
    output = (x,y,rad)
    return output

def makePrint(min_minutia, max_minutia, xlim, ylim):
    output = []
    num_minutia = int(np.random.uniform(min_minutia, max_minutia))
    for i in range(num_minutia):
        output.append(makeMinutia(xlim, ylim))
    return output

def makeNumPrints(nPrints, min_minutia, max_minutia, xlim, ylim):
    output = []
    for i in range(nPrints):
        output.append(makePrint(min_minutia, max_minutia, xlim, ylim))
    return output

prints = makeNumPrints(3, 3, 5, (0,10),(0,10))

# print (prints)


# print ("minutia list: ", prints[0])
