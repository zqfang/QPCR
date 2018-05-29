'''
This script used for read output from ABI qPCR machine.
Extract Data for further computing
'''

from __future__ import print_function
import argparse
import os
import sys
import pandas as pd

# parse command line args
parser = argparse.ArgumentParser(description="Extract qpcr data from ABi machine")

parser.add_argument("-m","--machine",choices=['viia7','q7','7900'],default="viia7",dest="machine",
                     help="Choose platform which data was generated from: viia7 or 7900. default: viia7.")
parser.add_argument("-d","--data", action="store", dest="data", required=True,
                    help="the excel file you want to analysis ")
parser.add_argument("-s","--sheetName", action="store",default="Results", dest="sheet",
                     help="the sheet name of your excel file you want to analysis ")
parser.add_argument("--header", action="store",type=int,dest="head", default=35,
                     help="header row you want to start with. Default: 35")
parser.add_argument("--tail",action="store",type=int,dest="tail",default=5,
                     help="the tail rows of your excel file you want to drop, default: 5")
# parser.add_argument("-r","--referenceControl", action="store", default="GAPDH", dest="rc",
#                     help="the reference gene name of your sample, default: GAPDH")
# parser.add_argument("-c","--experimentalControl",action="store",dest="ec",
#                     help="the control group name which your want to compare, e.g. hESC")
parser.add_argument("-o","--outFileNamePrefix",action="store",default="foo",dest="out",
                    help="the output file name")
parser.add_argument("--version",action="version",version="%(prog)s 1.0")
args = parser.parse_args()

print("ExeclFile        =", args.data)
print("Platform         =", args.machine)
print("SheetName        =", args.sheet)
print("headerRow        =", args.head)
print("tailRow          =", args.tail)
# print("ReferenceControl =", args.rc)
# print("ExperimentControl=", args.ec)
print("outFileName      =", args.out)

col_viia = ['Sample Name','Target Name','CT','Ct Mean','Ct SD']
col_7900 = ['Sample Name','Detector Name','Ct','Ct Mean','Ct StdEV']
# checking flies and parameters.
if not os.path.exists(args.data) :
   print("InputFile doesn't exist, please check your file path!")
   sys.exit(1)

print("Input File Checking passed !")


# Read data into pandas DataFrame Object
if args.machine in ['viia7','q7']:
    if args.machine == 'q7': args.head = 44  
    data = pd.read_excel(args.data, sheet_name=args.sheet, comment='#',
                         header= args.head, skipfooter=args.tail)
    dat = data[col_viia]
elif args.machine == '7900':
    args.head, args.tail = 10, 7
    data = pd.read_table(args.data, comment='#',
                         header= args.head, skipfooter=args.tail)
    dat = data[col_7900]
else:
    print("-m args error, please refine your args")
    sys.exit(1)

# write to a new excel file
print("Writing out put files.")
dat.to_csv("%s_%s.txt"%(args.out, args.machine),sep="\t",index=False,)
#
print("Done!")
