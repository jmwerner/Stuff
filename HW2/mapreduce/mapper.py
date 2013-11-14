#!/usr/bin/env python

import math
import sys


for line in sys.stdin:
    #Remove leading and trailing whitespace
    line = line.strip()
    #Split the input line
    split = line.split()

    x_lo = math.floor(float(split[0])*10)/10
    x_hi = x_lo + .1
    y_lo = math.floor(float(split[1])*10)/10
    y_hi = y_lo + .1

    print '%.1f,%.1f,%.1f,%.1f\t%i' % (x_lo, x_hi, y_lo, y_hi, 1)

exit()