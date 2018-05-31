'''
This script used for calculate delta Ct, delta delta Ct, fold changes and p-values.
version: 0.1
date: 2018.05.28
author: Zhuoqing Fang
email: fzq518@gmail.com
'''

from __future__ import print_function
import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations


def parse_cli():
    parser = argparse.ArgumentParser(description="Calculate Delta Ct, DDelta Ct, Fold Changes, P-values for QPCR results.")

    parser.add_argument("-d", "--data", action="store", dest="data", required=True,
                        help="the file(s) you want to analysis. For multi-file input, separate each file by comma. ")
    parser.add_argument("-s", "--sheetName", action="store", default=0, dest="sheet",
                        help="str, int. the sheet name of your excel file you want to analysis." + \
                             "Strings are used for sheet names, Integers are used in zero-indexed sheet positions.")
    parser.add_argument("-i", "--internalControl", action="store", dest="ic", default=None,
                        help="the internal control gene name of your sample, e.g. GAPDH")
    parser.add_argument("-e", "--experimentalControl", action="store", dest="ec",
                        required=True, help="the control group name which your want to compare, e.g. hESC")
    parser.add_argument("-o", "--outFileNamePrefix", action="store", default="foo", dest="out",
                        help="the output file name")
    parser.add_argument("-m", "--mode", action="store", dest="mode", type=str,
                        choices=("bioRep", "techRep", "dropOut","stat"), default="stat",
                        help="calculation mode. Choose from {'bioRep', 'techRep','dropOut'.'stat'}." + \
                             "'bioRep': using all data to calculate mean DeltaCT.." + \
                             "'techRep': only use first entry of replicates." + \
                             "'dropOut': if sd < 0.5, reject outlier and recalculate mean CT."+\
                             "'stat': statistical testing for each group vs experimental control. Default: dropOut.")
    parser.add_argument("--std", action="store", type=float, dest="std", default=0.5,
                        help="Minimal standard deviation to filtering CT values. Only affects CT means. Default: 0.5")
    parser.add_argument("--header", action="store", type=int, dest="head", default=0,
                        help="Row (0-indexed) to use for the column labels of the parsed DataFrame")
    parser.add_argument("--tail", action="store", type=int, dest="tail", default=0,
                        help="the tail rows of your excel file you want to skip (0-indexed)")

    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    args = parser.parse_args()

    return args

def read_input(args):
    # checking flies and parameters.
    if not os.path.exists(args.path):
        print("InputFile: %s doesn't exist, please check your file path!"%args.path)
        sys.exit(1)
    # read data into a dataFrame
    suffix = args.path.split(".")[-1]
    if suffix in ['xls', 'xlsx']:
        data = pd.read_excel(io=args.path, sheet_name=args.sheet, comment='#',
                             header=args.head, skipfooter=args.tail)
    elif suffix == 'csv':
        data = pd.read_csv(args.path, comment='#', header=args.head, skipfooter=args.tail)
    elif suffix == 'txt':
        data = pd.read_table(args.path, comment='#', header=args.head, skipfooter=args.tail)
    else:
        print("Unsupported file format input. Please use xls, xlsx. csv, txt format!")
        sys.exit(1)

    # rename column name
    if 'Target Name' not in data:
        if 'Detector Name' in data:
            # ABI 7900
            data.rename(columns={'Ct': 'CT', 'Detector Name': 'Target Name', 'Ct StdEV': 'Ct SD'}, inplace=True)
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
    if args.ec not in data['Sample Name'].unique(): 
        print("Experimental control Name: %s Not found in %s!"%(args.ec,args.path))
        sys.exit(1)
    if args.mode != 'stat':    
        if args.ic not in data['Target Name'].unique(): 
            print("intternal control Name:%s Not found in %s!"%(args.ic,args.path) )
            sys.exit(1)
    # convert "undetermined to NA"
    data.iloc[:, 2:] = data.iloc[:, 2:].apply(pd.to_numeric, errors='coerce', axis=1)

    return data


def parse_input(args):
    """parse multi-file input"""
    if args.std <=0:
        print("standard deviation could not be negative!!!")
        sys.exit(1)

    path = args.data.split(",")
    data =[]
    for p in path:
        args.path = p
        data.append(read_input(args))
    args.groups = data
    return args

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.nanmedian(data))
    mdev = np.nanmedian(d)
    s = d/(mdev if mdev else 1.)
    return data[s<m]


def min_mean(arr, std=0.5):
    arr_std = np.nanstd(arr,ddof=1)
    if arr_std < std:
        mmean = np.nanmean(arr)
    else:
        temp = reject_outliers(arr)
        mmean = temp.mean()
    return mmean

def min_mean2(arr, std=0.5):
    na = np.isnan(arr)
    na_num,arr_len = na.sum(),len(arr)
    if na_num == arr_len: return np.NaN
    arr_std = np.nanstd(arr)
    if arr_std < std:
        mmean = np.nanmean(arr)
    else:
        mmean = np.array(list(combinations(arr[~na], arr_len-na_num))).mean(axis=0).min()
    return mmean


def stars(p):
   if p < 0.0001:
       return "****"
   elif (p < 0.001):
       return "***"
   elif (p < 0.01):
       return "**"
   elif (p < 0.05):
       return "*"
   else:
       return "-"

def ttest(arr1, arr2, axis=1):
    """arr1, arr2 could be two dimension. refer to scipy.stats.ttest_ind
       for pandas dataframe input, note the axis

    """
    # z, pm = stats.mannwhitneyu(arr1, arr2)
    # pm = pm * 2 # two tailed
    t, p = stats.ttest_ind(arr1,arr2, axis=axis) # p, already two tailed
    # s = stars(p)
    s = [stars(pv) for pv in p]
    return t,p,s

def reshape(df):


    df['CTs'] = 'CT Rep' + df.groupby(['Sample Name', 'Target Name'], as_index=False).cumcount().astype(str)
    df_pivot = df.pivot_table(index=['Sample Name','Target Name'],columns='CTs',values='CT')
    df_stat = df.groupby(['Sample Name', 'Target Name'])['CT'].agg(['mean','std'])
    df_stat.rename(columns={'mean': 'Ct Mean (old)','std': 'Ct SD'}, inplace=True)
    df2 = df_pivot.merge(df_stat, left_index=True, right_index=True)
    return df2


def run(args):
    """run program"""

    if args.mode != 'stat':
        if args.ic is None:
            print('"-i", "--internalControl" is required')
            sys.exit(1)
    
    # run mode
    outname = args.out
    idxs = len(args.groups)
    if args.mode == 'bioRep':
        for idx, data in enumerate(args.groups):
            if idxs >1: args.out = outname+str(idx)
            args.ct = reshape(data)
            data2 = data.groupby(['Sample Name', 'Target Name'])['CT'].mean()
            data2 = pd.DataFrame(data2)
            data2.rename(columns={'CT': 'Ct Mean'}, inplace=True)
            calc_fc(args, data2)
    elif args.mode == 'techRep':
        # instead of using groupby, you can remove duplicates using drop_duplicates
        # tech replicates need to drop outliers by hand
        # this means you already have a 'CT mean' column exists
        for idx, data in enumerate(args.groups):
            if idxs >1: args.out = outname+str(idx)
            args.ct = reshape(data)
            data2 = data.drop_duplicates(['Sample Name', 'Target Name'])
            data2 = data2.set_index(['Sample Name', 'Target Name'])
            calc_fc(args, data2)
    elif args.mode == 'dropOut':
        # find and drop outliers automatically, the calculate CT mean...
        for idx, data in enumerate(args.groups):
            if idxs >1: args.out = outname+str(idx)
            args.ct = reshape(data)
            data2 = data.groupby(['Sample Name', 'Target Name'])['CT'].apply(min_mean, std=args.std)
            data2 = pd.DataFrame(data2)
            data2.rename(columns={'CT': 'Ct Mean'}, inplace=True)
            calc_fc(args, data2)
    elif args.mode == 'stat':
        args.df = pd.concat(args.groups)
        if 'Delta Ct' in args.df:
            calc_stats(args)
        else:
            args.mode = 'dropOut'
            run(args)
    else:
        print("No supported method for further calculation")
        sys.exit(1)

def calc_fc(args, df):
    """calculate Delta Ct, Fold Changes for one experiment with tech replicates"""
    sample =  df.index.get_level_values(0).unique().tolist()
    print("Your Samples are: ")
    print(sample)

    # calculate Delta_Ct value.
    DelCt = pd.DataFrame()
    for i in range(len(sample)):
        deltaCt = df.loc[sample[i]] - df.loc[(sample[i], args.ic)]
        deltaCt['Sample Name'] = sample[i]
        DelCt = DelCt.append(deltaCt)

    # reshape your dataFrame,export to a csv file
    DelCt.index.name = 'Target Name'
    DelCt3 = DelCt.reset_index().set_index(['Sample Name', 'Target Name'])
    # rename column names
    DelCt3.rename(columns={'Ct Mean': 'Delta Ct'}, inplace=True)

    # calculate Delta_Delta_Ct
    DDelCt = pd.DataFrame()
    for i in range(len(sample)):
        DDeltaCt = DelCt3.loc[sample[i]] - DelCt3.loc[args.ec]
        DDeltaCt['Sample Name'] = sample[i]
        DDelCt = DDelCt.append(DDeltaCt)

    # reshape your dataFrame
    DDelCt.index.name = 'Target Name'
    DDelCt3 = DDelCt.reset_index().set_index(['Sample Name', 'Target Name'])
    DDelCt3.rename(columns={'Delta Ct': 'DDelta Ct'}, inplace=True)
    DDelCt4 = DDelCt3.dropna(how='all')

    # calculate FoldChange, and export to a csv file
    foldChange0 = np.power(2, -DDelCt4)
    foldChange0.rename(columns={'DDelta Ct': 'Fold Changes'}, inplace=True)
    # reshape your final results
    final_report = pd.concat([args.ct, df[['Ct Mean']], DelCt3[['Delta Ct']],
                              DDelCt4[['DDelta Ct']], foldChange0[['Fold Changes']]], axis=1, )
    final_report.to_csv(args.out + '.csv')
    print("\n\n########################################")
    print('The first 8 rows in your Final results:\n')
    print(final_report.head(8))


def calc_stats(args):
    """calculate Delta Ct, Fold Changes for n independent experiments with biological replicates"""

    # ref_ctrl = args.ic
    # exp_ctrl = args.ec
    data = args.df
    data['Replicates'] = 'Rep' + data.groupby(['Sample Name', 'Target Name']).cumcount().astype(str)
    df_ctrl = data[data['Sample Name'] == args.ec]
    DCt_ctrl = df_ctrl.groupby(['Target Name'])['Delta Ct'].mean()

    sample = data['Sample Name'].unique().tolist()
    print("Your Samples are: ")
    print(sample)

   # calculate Delta_Delta_Ct
    df = data.set_index(['Sample Name', 'Target Name'])['Delta Ct']
    DDelCt = pd.DataFrame()
    for i in range(len(sample)):
        DDeltaCt = df.loc[sample[i]] - DCt_ctrl
        DDeltaCt = DDeltaCt.reset_index()
        DDeltaCt['Sample Name'] = sample[i]
        DDelCt = DDelCt.append(DDeltaCt)

    # reshape your dataFrame
    DDelCt['Replicates'] = 'Rep' + DDelCt.groupby(['Sample Name', 'Target Name']).cumcount().astype(str)
    DDelCt = DDelCt.set_index(['Sample Name', 'Target Name','Replicates'])
    DDelCt.rename(columns={'Delta Ct': 'DDelta Ct'}, inplace=True)
    DDelCt = DDelCt.dropna(how='all')
    # calculate FoldChange, and export to a csv file
    foldChange0 = np.power(2, -DDelCt)
    foldChange0.rename(columns={'DDelta Ct': 'Fold Changes'}, inplace=True)

    # reshape your final results
    # extract columns you needed,remain as DataFrame object using double bracket.
    # for merging columns easily
    df2 = data.set_index(['Sample Name', 'Target Name','Replicates'])[['Ct Mean', 'Delta Ct']]
    final = pd.concat([df2, DDelCt, foldChange0], axis=1,)
    #  statistical testing
    # final = final.reset_index()
    # final_pivot = final.loc[:,['Sample Name', 'Target Name', 'Replicates', 'Fold Changes']]
    # final_pivot = final_pivot.pivot_table(index=['Sample Name','Target Name'],
    #                                      columns='Replicates',values='Fold Changes')
    final_pivot = final.unstack(level=2)
    t_stat = pd.DataFrame()
    test0 = final_pivot.loc[args.ec, ['Fold Changes']]
    for i in range(len(sample)):
        # if sample[i] == exp_ctrl: continue
        test1 = final_pivot.loc[sample[i], ['Fold Changes']]
        t, p, s = ttest(test1, test0, axis=1)
        test1['Sample Name'] = sample[i]
        test1['T-statistic'] = t
        test1['P-value'] = p
        test1['Stars'] = s
        t_stat = t_stat.append(test1)
    # column names
    t_stat = t_stat.reset_index().set_index(['Sample Name', 'Target Name']).reset_index()
    # write output
    writer = pd.ExcelWriter('%s_stats.xls'%args.out)
    final.reset_index().to_excel(writer, 'Full')
    t_stat.to_excel(writer,'Stats')
    writer.save()
    print("\n\n#######################################")
    print('The first 8 rows in your Final results:\n')
    print(t_stat.head(8))


if __name__ == "__main__":
    args = parse_cli()
    print("###########################################")
    print("Input parameters:")
    print("InputFile        =", args.data)
    print("SheetName        =", args.sheet)
    print("headerRow        =", args.head)
    print("tailRow          =", args.tail)
    print("ReferenceControl =", args.ic)
    print("ExperimentControl=", args.ec)
    print("outFileName      =", args.out)
    print("\n###########################################")

    # run program
    run(parse_input(args))
    print("\nProgram ran successfully")
    print("Good Job! Cheers!!!")
