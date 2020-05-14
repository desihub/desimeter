"""
Utility functions
"""

import numpy as np

def parse_fibers(fiber_string) :
    """
    Short func that parses a string containing a comma separated list of
    integers, which can include ":" or ".." or "-" labeled ranges

    Args:
        fiber_string (str) : list of integers or integer ranges

    Returns (array 1-D):
        1D numpy array listing all of the integers given in the list,
        including enumerations of ranges given.

    Note: this follows python-style ranges, i,e, 1:5 or 1..5 returns 1, 2, 3, 4
    """
    if fiber_string is None :
        return np.array([])
    else:
        fiber_string = str(fiber_string)

    if len(fiber_string.strip(' \t'))==0:
        return np.array([])

    fibers=[]


    for sub in fiber_string.split(',') :
        sub = sub.replace(' ','')
        if sub.isdigit() :
            fibers.append(int(sub))
            continue

        match = False
        for symbol in [':','..','-']:
            if not match and symbol in sub:
                tmp = sub.split(symbol)
                if (len(tmp) == 2) and tmp[0].isdigit() and tmp[1].isdigit() :
                    match = True
                    for f in range(int(tmp[0]),int(tmp[1])) :
                        fibers.append(f)

        if not match:
            print("parsing error. Didn't understand {}".format(sub))
            import sys
            sys.exit(1)

    return np.array(fibers)
