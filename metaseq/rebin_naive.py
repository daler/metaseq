from numpy import *
from math import floor

#          0   1   2   3   4   5   6   7   8   9
y = array([0,  0,  1,  2,  2,  1,  1,  0,  0,  0])
# bins=5     0       1       2       3       4
# bins=4     0       1       2       3
# bins=20  01  23  45  67  89  01   23 45  67  89
oldbins = arange(len(y))
newbins = arange(5)
expected = array([0, 3, 3, 1, 0])



def f(y, bins):
    yi = zeros( (bins,))
    spacing = len(y) / bins
    if spacing > 0:
        for i in range(bins):
            #print
            i1 = i * spacing
            i2 = i1 + spacing
            s = sum(y[i1:i2])
            #print i, y[i1:i2], '->', s,
            yi[i] = s
        if i2 < len(y):
            s = sum(y[i2:])
            yi[-1] += s
            #print '+=',  s, '->', yi[-1]
    else:
        # we need to expand...
        # Convention will be to use the lowest bin.
        padding = int(floor(1. / (len(y) / float(bins))))
        c = 0
        #print
        for i in range(len(y)):
            yi[c] = y[i]
            c += 1
            #print i * padding, y[i], '->', y[i]
            for j in range(padding-1):
                if c < bins:
                    yi[c] = 0
                    c += 1
                    #print i * padding + j, '0 -> 0'
        while c < bins:
            yi[c] = 0
            c += 1
            #print i * padding + j, '0 -> 0'

    return array(yi)

def p(x):
    return str(int(x))

yi = f(y, 5)
assert all(yi == array([0, 3, 3, 1, 0]))
print
yi = f(y, 4)
assert all(yi == array([0, 3, 3, 1]))

# how to handle "expansions"?

print
yi = f(y, 20)
print ' '.join(map(p,yi))
assert all(yi == array([0,0,0,0,1,0,2,0,2,0,1,0,1,0,0,0,0,0,0,0]))

print
yi = f(y, 22)
print ' '.join(map(p,yi))
assert all(yi == array([0,0,0,0,1,0,2,0,2,0,1,0,1,0,0,0,0,0,0,0,0,0]))
