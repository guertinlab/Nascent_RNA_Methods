#! /usr/bin/python3
import sys
import getopt
import pandas as pd
import os
import itertools
import random
import os
import glob
from decimal import *

def function(infile, scale, outfilename):
    df_in = pd.read_csv(infile, sep="\t", header=None)
    df_in[3] = df_in[3].multiply(float(scale))
    df_in[3] = df_in[3].abs()
    df_in.to_csv(outfilename, sep="\t", index=False, header=False)
    #outfile = open(outfilename, 'w')
    #for line in open(infile):
    #        #print line
    #    splitline=line.split()
    #    signal=abs(Decimal(splitline[3])*Decimal(scale))
    #    if ( signal > 0 ):
    #        outfile.write("%s\t%s\t%s\t%s\n"%(splitline[0], splitline[1], splitline[2], signal))
    #outfile.close()
    
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "i:s:o:h", ["infile=", "scale=", "outfilename=","help"])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    infile = False
    scale = False
    outfilename = False
    for opt, arg in opts:
        if opt in ('-i', '--infile'):
            infile = arg
        elif opt in ('-s', '--scale'):
            scale = arg
        elif opt in ('-o', '--outfilename'):
            outfilename = arg
        elif opt in ('-h', '--help'):
            print('\nnormalize_bedGraph.py -i HEK293T_TIR1_Cl4_rep1_minus_body_0-mer.bg -s 1.2 -o HEK293T_TIR1_Cl4_rep1_minus_body_0-mer.scaled.bg')
            sys.exit()
    if infile and scale and outfilename:
        print(infile)
        function(infile, scale, outfilename)
if __name__ == "__main__":
    main(sys.argv[1:])
