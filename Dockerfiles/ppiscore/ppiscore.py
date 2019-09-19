#!/usr/bin/env python3
import numpy as np
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python3 %s <MEGADOCK output filename>\n" % sys.argv[0])
        sys.exit(1)
    with open(sys.argv[1], "r") as f:
        next(f) # skip comment lines
        next(f) # skip comment lines
        next(f) # skip comment lines
        next(f) # skip comment lines
        data = np.array(list(map(lambda s: s.split()[6], f)), dtype=np.float32)
        std = np.std(data)
        mean = np.mean(data)
        maximum = np.max(data)
        score = (maximum - mean) / std
        pair_columns = sys.argv[1].rsplit('/', maxsplit=1)[1][:-4].replace('-', '\t').replace('_', '\t')
        
        # print("%s\t%f" % (pair_columns, score))
        print("%s, E = %f" % (sys.argv[1], score))