'''
This script used for calculateing DeltaCt, DeltaDeltaCt, FoldChange for QPCR 
'''
from __future__ import print_function
import argparse
import sys
import os
import pandas as pd

########################### parse command line args############################################
parser = argparse.ArgumentParser(description="Calculating DeltaCt, DeDeltaCt, FoldChange for QPCR results.")

parser.add_argument("-d","--data", action="store", dest="data", required=True, 
                    help="the file you want to analysis. ")
parser.add_argument("-s","--sheetName", action="store",default="Results", dest="sheet", required=True,
                     help="the sheet name of your excel file you want to analysis ")
parser.add_argument("-i","--internalControl", action="store", dest="ic", 
                    required=True, help="the internal control gene name of your sample, e.g. GAPDH")
parser.add_argument("-e","--experimentalControl", action="store", dest="ec",
                    required=True, help="the control group name which your want to compare, e.g. hESC")
parser.add_argument("-o","--outFileNamePrefix", action="store", default="foo", dest="out",\
                     help="the output file name")
parser.add_argument("-m", "--mode", action="store", dest="mode", type=str,
                    choices=("bioRep", "techRep" ), default="techRep", 
                    help="calculation mode. Choose from {'bioRep', 'techRep'}."+\
                          "bioRep: using all data to caluclate mean DeltaCT, techRep: only use first entry of replicates. Default: techRep.")
parser.add_argument("--header", action="store", type=int, dest="head", default=0,\
                     help="header row you want to start with")
parser.add_argument("--tail", action="store",type=int, dest="tail", default=0,\
                     help="the tail rows of your excel file you want to drop")

parser.add_argument("--version",action="version",version="%(prog)s 1.0")
args = parser.parse_args()

print("InputFile        =", args.data)
print("SheetName        =", args.sheet)
print("headerRow        =", args.head)
print("tailRow          =", args.tail)
print("ReferenceControl =", args.ic)
print("ExperimentControl=", args.ec)
print("outFileName      =", args.out)

####################################################################


ref_ctrl = args.ic
exp_ctrl = args.ec


# checking flies and parameters.
if not os.path.exists(args.data) :
  print("InputFile doesn't exist, please check your file path!")
  sys.exit(1)
# read data into a dataFrame
data = pd.read_excel(args.data, args.sheet)


# rename column name
if 'Target Name' not in data:
    if 'Detector Name' in data:
        data.rename(columns={'Ct': 'CT','Detector Name':'Target Name','Ct StdEV':'Ct SD'}, inplace=True)
    else:
        print("column name error! Plesase rename your column name!")
        sys.exit(1)


# calculate Ct mean values for each replicates
if args.mode == 'bioRep':
    data2 = data.groupby(['Sample Name','Target Name']).mean().set_index(['Sample Name','Target Name'])
else:

# instead of using groupby, you can remove duplicates using drop_duplicates
    data2 = data.drop_duplicates(['Sample Name','Target Name']).set_index(['Sample Name','Target Name'])


sample = data['Sample Name'].unique()
print("Your Samples are: ")
print(sample)

  
#calculate Delta_Ct value.
DelCt=pd.DataFrame()
for i in range(len(sample)):
    deltaCt = data2.loc[sample[i]]- data2.loc[(sample[i], ref_ctrl)]
    deltaCt['Sample Name'] = sample[i]
    DelCt = DelCt.append(deltaCt)
       

#reshape your dataFrame,export to a csv file
DelCt.index.name = 'Target Name'
DelCt2 = DelCt.reset_index()
DelCt3 = DelCt2.set_index(['Sample Name','Target Name'])

#rename column names
DelCt3.rename(columns={'Ct Mean': 'DeltaCt'}, inplace=True)


#calculate Delta_Delta_Ct
DDelCt =pd.DataFrame()
for i in range(len(sample)):
    DDeltaCt = DelCt3.loc[sample[i]]-DelCt3.loc[exp_ctrl]
    DDeltaCt['Sample Name'] = sample[i]
    DDelCt = DDelCt.append(DDeltaCt)
    

#reshape your dataFrame
DDelCt.index.name = 'Target Name'
DDelCt2 = DDelCt.reset_index()
DDelCt3 = DDelCt2.set_index(['Sample Name','Target Name'])
DDelCt3.rename(columns={'DeltaCt': 'DDeltaCt'}, inplace=True)
DDelCt4=DDelCt3.dropna(how='all')


# calculate FoldChange, and export to a csv file
foldChange0 = pow(2,-DDelCt4)
foldChange0.rename(columns={'DDeltaCt':'FoldChange'}, inplace=True)
foldChange = foldChange0.drop(exp_ctrl,level=0)
foldChange1 = foldChange0.drop(ref_ctrl,level=1)

  
#reshape your final results
#extract columns you needed,remain as DataFrame object using double bracket.
#for merging columns esayliy

#merge data,export to a csv file.
m1=pd.merge(data2[['Ct Mean']], DelCt3[['DeltaCt']], left_index=True, right_index=True)
m2=pd.merge(m1, DDelCt4[['DDeltaCt']], left_index=True, right_index=True)
final_report=pd.merge(m2, foldChange0[['FoldChange']], left_index=True, right_index=True)
final_report.index.names=['Sample Name','Target Name']
final_report.to_csv(args.out+'_final_results.csv')

print('The first 5 rows in your Final results: ')
print(final_report.head(5))

print("Program ran successfully")
print("Good Job! Cheers!!!" )

