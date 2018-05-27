from __future__ import print_function
import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats



def parse_cli():
    parser = argparse.ArgumentParser(description="Calculate Delta Ct, DDelta Ct, Fold Changes for QPCR results.")

    parser.add_argument("-d", "--data", action="store", dest="data", required=True,
                        help="the file you want to analysis. ")
    parser.add_argument("-s", "--sheetName", action="store", default=0, dest="sheet",
                        help="str, int. the sheet name of your excel file you want to analysis." + \
                             "Strings are used for sheet names, Integers are used in zero-indexed sheet positions.")
    parser.add_argument("-i", "--internalControl", action="store", dest="ic",
                        required=True, help="the internal control gene name of your sample, e.g. GAPDH")
    parser.add_argument("-e", "--experimentalControl", action="store", dest="ec",
                        required=True, help="the control group name which your want to compare, e.g. hESC")
    parser.add_argument("-o", "--outFileNamePrefix", action="store", default="foo", dest="out",
                        help="the output file name")
    parser.add_argument("-m", "--mode", action="store", dest="mode", type=str,
                        choices=("bioRep", "techRep", "dropOut","stat"), default="bioRep",
                        help="calculation mode. Choose from {'bioRep', 'techRep','dropOut'}." + \
                             "'bioRep': using all data to calculate mean DeltaCT." + \
                             "'techRep': only use first entry of replicates." + \
                             "'dropOut': if sd < 0.5, reject outlier and recalculate mean CT."+\
                             "'stat': statistical testing for each group and target names. Default: bioRep.")
    parser.add_argument("--header", action="store", type=int, dest="head", default=0,
                        help="Row (0-indexed) to use for the column labels of the parsed DataFrame")
    parser.add_argument("--tail", action="store", type=int, dest="tail", default=0,
                        help="the tail rows of your excel file you want to skip (0-indexed)")

    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    args = parser.parse_args()

    return args


def parse_input(args):
    # checking flies and parameters.
    if not os.path.exists(args.data):
        print("InputFile doesn't exist, please check your file path!")
        sys.exit(1)
    # read data into a dataFrame
    suffix = args.data.split(".")[-1]
    if suffix in ['xls', 'xlsx']:
        data = pd.read_excel(io=args.data, sheet_name=args.sheet,
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
    # convert "undetermined to NA"
    data.iloc[:, 2:] = data.iloc[:, 2:].apply(pd.to_numeric, errors='coerce', axis=1)

    args.df = data
    return args

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

def calculate(args):

    ref_ctrl = args.ic
    exp_ctrl = args.ec
    data = args.df
    # calculate Ct mean values for each replicates
    # if args.mode == 'bioRep':
    #     data2 = data.groupby(['Sample Name', 'Target Name']).mean()
    # elif args.mode == 'techRep':
    #     # instead of using groupby, you can remove duplicates using drop_duplicates
    #     # tech replicates need to drop outliers
    #     # this means you already have a 'CT mean' column exists
    #     data2 = data.drop_duplicates(['Sample Name', 'Target Name']).set_index(['Sample Name', 'Target Name'])
    # elif args.mode == 'dropOut':
    #     # find to drop outliers in each data point, the calculate CT mean
    #     data2 = data.groupby(['Sample Name', 'Target Name'])['CT'].apply(min_mean2)
    #     data2 = pd.DataFrame(data2)
    #     data2.rename(columns={'CT': 'Ct Mean'}, inplace=True)
    # elif args.mode == 'stat':
    #     # get ctrl Delta CT mean value
    #     df_ctrl = data[data['Sample Name'] == exp_ctrl]
    #     DCt_ctrl = df_ctrl.groupby(['Sample Name', 'Target Name']).mean()
    #
    # else:
    #     print("No supported method for further calculation")
    #     sys.exit(1)

    data['Replicates'] = 'Rep' + data.groupby(['Sample Name', 'Target Name']).cumcount().astype(str)
    df_ctrl = data[data['Sample Name'] == exp_ctrl]
    DCt_ctrl = df_ctrl.groupby(['Target Name'])['Delta Ct'].mean()

    sample = data['Sample Name'].unique()
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
    test0 = final_pivot.loc[exp_ctrl, ['Fold Changes']]
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
    t_stat = t_stat.reset_index().set_index(['Sample Name', 'Target Name'])
    # write output
    writer = pd.ExcelWriter('%s_final_stats.xls'%args.out)
    final.to_excel(writer, 'Full')
    t_stat.to_excel(writer,'Stats')
    writer.save()
    print('The first 8 rows in your Final results: ')
    print(t_stat.head(8))


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
