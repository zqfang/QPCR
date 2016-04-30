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

parser.add_argument("-m","--machine",choices=['viia7','7900'],default="viia7",dest="machine",\
                     help="Choose platform which data was generated from: viia7 or 7900. default: viia7.")
parser.add_argument("-d","--data", action="store", dest="data",help="the excel file you want to analysis ")
parser.add_argument("-s","--sheetName", action="store",default="Results", dest="sheet", \
                     help="the sheet name of your excel file you want to analysis ")
parser.add_argument("--header", action="store",type=int,dest="head", default=35,\
                     help="header row you want to start with. Default: 35")
parser.add_argument("--tail",action="store",type=int,dest="tail",default=5,\
                     help="the tail rows of your excel file you want to drop, default: 5")
parser.add_argument("-r","--referenceControl", action="store", default="GAPDH", dest="rc", \
                     help="the reference gene name of your sample, default: GAPDH")
parser.add_argument("-c","--experimentalControl",action="store",dest="ec",\
                     help="the control group name which your want to compare, e.g. hESC")
parser.add_argument("-o","--outFileNamePrefix",action="store",default="foo",dest="out",\
                     help="the output file name")
parser.add_argument("--version",action="version",version="%(prog)s 1.0")
args = parser.parse_args()

print("ExeclFile        =", args.data)
print("Platform         =", args.machine)
print("SheetName        =", args.sheet)
print("headerRow        =", args.head)
print("tailRow          =", args.tail)
print("ReferenceControl =", args.rc)
print("ExperimentControl=", args.ec)
print("outFileName      =", args.out)

col_viia = ['Sample Name','Target Name','CT','Ct Mean','Ct SD']
col_7900 = ['Sample Name','Detector Name','Ct','Ct Mean','Ct StdEV']
# checking flies and parameters.
if not os.path.exists(args.data) :
  print("InputFile doesn't exist, please check your file path!")
  sys.exit(1)

print("Input File Checking passed !")


# Read data into pandas DataFrame Object
if args.machine == 'viia7':
    data0 = pd.read_excel(args.data, sheetname=args.sheet, header= args.head)
    
elif args.machine == '7900':
    data0 = pd.read_table(args.data, header= args.head)
else:
    print("-m args error, plesea refine your args")
    sys.exit(1)

data = data0.drop(data0.tail(args.tail).index) # remove last 5 row

if args.machine == 'viia7':
   dat = data[col_viia] # extract data wich we are interest
else:
   dat = data[col_7900]

# write to a new excel file
print("Writing out put files.")

dat.to_csv(args.out+"_interest_data.csv",index=False, na_values = ['NA'])
data.to_csv(args.out+".csv",index=False, na_values = ['NA'])

#
print("Done!")
