'''
This script used for read excel output from ABI ViiA7 qPCR machine.
Extract Data for further computing
'''
from __future__ import print_function
import argparse


# parse command line args
parser = argparse.ArgumentParser(description="Extract qpcr data from ABi machine")

parser.add_argument("-m","--machine",choices=['viia7','7900'],default="viia7",dest="machine",\
                     help="Choice platform which data was generated from.")
parser.add_argument("-f","--file", action="store", dest="file",help="the excel file you want to analysis ")
parser.add_argument("-s","--sheetName", action="store",default="Results", dest="sheet", \
                     help="the sheet name of your excel file you want to analysis ")
parser.add_argument("--header", action="store",type=int,dest="head", default=35,\
                     help="header row you want to start with")
parser.add_argument("--tail",action="store",type=int,dest="tail",default=5,\
                     help="the tail rows of your excel file you want to drop")
parser.add_argument("-r","--referenceControl", action="store", default="GAPDH", dest="rc", \
                     help="the reference gene name of your sample, e.g. GAPDH")
parser.add_argument("-c","--experimentalControl",action="store",dest="ec",\
                     help="the control group name which your want to compare, e.g. hESC")
parser.add_argument("-o","--outFileNamePrefix",action="store",default="foo",dest="out",\
                     help="the output file name")
parser.add_argument("--version",action="version",version="%(prog)s 1.0")
args = parser.parse_args()

print("ExeclFile        =", args.file)
print("Machine          =", args.machine)
print("SheetName        =", args.sheet)
print("headerRow        =", args.head)
print("tailRow          =", args.tail)
print("ReferenceControl =", args.rc)
print("ExperimentControl=", args.ec)
print("outFileName      =", args.out)

columns = ['Sample Name','Target Name','CT','Ct Mean','Ct SD']

# Read data into pandas DataFrame Object
data0 = pd.read_excel(args.file,sheetname=args.sheet,header= args.head)
data = data0.drop(data0.tail(tailRowDrop).index) # remove last 5 row
dat = data[columns] # extract data wich we are interest

# write to a new excel file
dat.to_excel(args.out+"_interest_data.xls",'Results',index_col= None, na_values = ['NA'])
data.to_excel(args.out+".xls",'Results',index_col= None, na_values = ['NA'])
