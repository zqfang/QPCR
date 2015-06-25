'''
This script used for read excel output from ABI ViiA7 qPCR machine.
Extract Data for further computing
'''

import pandas as pd
import numpy as np
import os
import getopt
import sys
####################################################################
def	show_help():
	print "USAGE:\tpython readExel.py <--input_file=input_file> [OPTIONS]"
	print "OPTIONS"
	print "-i or --input "
	sys.exit()


def read_cmd_line ():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["input_file=", "ic=", "interctrl=", \
                                                       "internal_control=", "ec=", "expCtrl=", \
                                                       "experiment_control=", "help"])
    except getopt.GetoptError:
        print "illegal params!"
        sys.exit()
	global input_file, internal_control, experiment_control
    for op, value in opts:
		if op == '--input_file':
			input_file=value
		if op == '--ic' or op == '--internal_control' or re.search('^--interCtrl', op, re.I) :
			internal_control=value
		if op == '--ec' or op == '--experiment_control' or re.search('^--expCtrl', op, re.I):
			experiment_control=value
		elif op == '--help' or op =='-h':
			show_help()
##############################################################################



# Params to specify 
input_file="/path/to/your/.xls/file"
output_file_name = "foo.xlsx"
sheetName="sheetname"
header=35  #header line begin
columns = ['Sample Name','Target Name','CT','Ct Mean','Ct SD']

# Read data into pandas DataFrame Object
data0 = pd.read_excel(file,sheetname=sheetName,header=header)
data = data0.drop(data0.tail(5).index) # remove last 5 row
dat = data[columns] # extract data wich we are interest

# write to a new excel file
dat.to_excel(output_file_name,'Results',index_col= None, na_values = ['NA'])
