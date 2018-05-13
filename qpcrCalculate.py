'''
This script used for calculating DeltaCt, DeltaDeltaCt, FoldChange for QPCR
'''
from __future__ import print_function
import argparse
import sys
import os
import numpy as np
import pandas as pd

########################### parse command line args############################################
def parse_cli():
    parser = argparse.ArgumentParser(description="Calculate Delta Ct, DDelta Ct, Fold Changes for QPCR results.")

    parser.add_argument("-d","--data", action="store", dest="data", required=True,
                        help="the file you want to analysis. ")
    parser.add_argument("-s","--sheetName", action="store",default=0, dest="sheet",
                         help="str, int. the sheet name of your excel file you want to analysis."+\
                               "Strings are used for sheet names, Integers are used in zero-indexed sheet positions.")
    parser.add_argument("-i","--internalControl", action="store", dest="ic",
                        required=True, help="the internal control gene name of your sample, e.g. GAPDH")
    parser.add_argument("-e","--experimentalControl", action="store", dest="ec",
                        required=True, help="the control group name which your want to compare, e.g. hESC")
    parser.add_argument("-o","--outFileNamePrefix", action="store", default="foo", dest="out",\
                         help="the output file name")
    parser.add_argument("-m", "--mode", action="store", dest="mode", type=str,
                        choices=("bioRep", "techRep","dropOutlier" ), default="bioRep",
                        help="calculation mode. Choose from {'bioRep', 'techRep','dropOutlier'}."+\
                              "bioRep: using all data to calculate mean DeltaCT.\n"+\
                              "techRep: only use first entry of replicates. Default: bioRep.")
    parser.add_argument("--header", action="store", type=int, dest="head", default=0,
                         help="Row (0-indexed) to use for the column labels of the parsed DataFrame")
    parser.add_argument("--tail", action="store",type=int, dest="tail", default=0,
                         help="the tail rows of your excel file you want to skip (0-indexed)")

    parser.add_argument("--version",action="version",version="%(prog)s 1.0")
    args = parser.parse_args()

    return args

def parse_input(args):
    # checking flies and parameters.
    if not os.path.exists(args.data) :
        print("InputFile doesn't exist, please check your file path!")
        sys.exit(1)
    # read data into a dataFrame
    suffix = args.data.split(".")[-1]
    if suffix in ['xls','xlsx']:
        data = pd.read_excel(io=args.data, sheetname=args.sheet,
                         header=args.head, skip_footer=args.tail)
    elif suffix == 'csv':
        data = pd.read_csv(args.data, header=args.head, skip_footer=args.tail)
    elif suffix == 'txt':
        data = pd.read_table(args.data, header=args.head, skip_footer=args.tail)
    else:
        print("Unsupported file format input. Please use xls,xlsx. csv, txt format!")
        sys.exit(1)

    # rename column name
    if 'Target Name' not in data:
        if 'Detector Name' in data:
            # ABI 7900
            data.rename(columns={'Ct': 'CT','Detector Name':'Target Name','Ct StdEV':'Ct SD'}, inplace=True)
        else:
            print("User defined Data Input.")
            print("You input data columns are:")
            print(data.columns.tolist())

            print("""
                  Column name error! Please rename your column name exactly like this:
    
                  | Sample Name | Target Name | CT | Ct Mean | Ct SD |
    
                  Note: these columns are optional:
                        | Ct SD |
                  """)
            sys.exit(1)


    args.df = data
    return args

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return data[s<m]


def mininal_min(arr):
    arr_std = np.std(arr)
    if arr_std < 0.5:
        mmin = np.mean(arr)
    else:
        temp = reject_outliers(arr)
        mmin = temp.mean()
    return mmin


def calculate(args):
    """ main program"""
    ref_ctrl = args.ic
    exp_ctrl = args.ec
    data = args.df
    # calculate Ct mean values for each replicates
    if args.mode == 'bioRep':
        data2 = data.groupby(['Sample Name','Target Name']).mean()
    elif args.mode == 'techRep':
        # instead of using groupby, you can remove duplicates using drop_duplicates
        # tech replicates need to drop outliers
        data2 = data.drop_duplicates(['Sample Name','Target Name']).set_index(['Sample Name','Target Name'])
    elif args.mode == 'dropOutlier':
        # find to drop outliers in each data point, the calculate CT mean
        data2 = data.groupby(['Sample Name', 'Target Name'])['CT'].apply(mininal_min)
        data2 = pd.DataFrame(data2)
        data2.rename(columns={'CT': 'Ct Mean'}, inplace=True)
    else:
        print("No supported method for further calculation")
        sys.exit(1)


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
    DelCt3 = DelCt.reset_index().set_index(['Sample Name','Target Name'])
    #rename column names
    DelCt3.rename(columns={'Ct Mean': 'Delta Ct'}, inplace=True)


    #calculate Delta_Delta_Ct
    DDelCt =pd.DataFrame()
    for i in range(len(sample)):
        DDeltaCt = DelCt3.loc[sample[i]]-DelCt3.loc[exp_ctrl]
        DDeltaCt['Sample Name'] = sample[i]
        DDelCt = DDelCt.append(DDeltaCt)


    #reshape your dataFrame
    DDelCt.index.name = 'Target Name'
    DDelCt3 = DDelCt.reset_index().set_index(['Sample Name','Target Name'])
    DDelCt3.rename(columns={'Delta Ct': 'DDelta Ct'}, inplace=True)
    DDelCt4=DDelCt3.dropna(how='all')


    # calculate FoldChange, and export to a csv file
    foldChange0 = pow(2,-DDelCt4)
    foldChange0.rename(columns={'DDelta Ct':'Fold Changes'}, inplace=True)
    foldChange = foldChange0.drop(exp_ctrl,level=0)
    foldChange1 = foldChange0.drop(ref_ctrl,level=1)


    #reshape your final results
    #extract columns you needed,remain as DataFrame object using double bracket.
    #for merging columns easily

    #merge data,export to a csv file.
    m1=pd.merge(data2[['Ct Mean']], DelCt3[['Delta Ct']], left_index=True, right_index=True)
    m2=pd.merge(m1, DDelCt4[['DDelta Ct']], left_index=True, right_index=True)
    final_report=pd.merge(m2, foldChange0[['Fold Changes']], left_index=True, right_index=True)
    final_report.index.names=['Sample Name','Target Name']
    final_report.to_csv(args.out+'_final_results.csv')
    print('The first 5 rows in your Final results: ')
    print(final_report.head(5))


if __name__ == "__main__":
    args = parse_cli()
    print("InputFile        =", args.data)
    print("SheetName        =", args.sheet)
    print("headerRow        =", args.head)
    print("tailRow          =", args.tail)
    print("ReferenceControl =", args.ic)
    print("ExperimentControl=", args.ec)
    print("outFileName      =", args.out)

    # run program
    calculate(parse_input(args))
    print("Program ran successfully")
    print("Good Job! Cheers!!!")
