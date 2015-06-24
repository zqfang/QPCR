'''
This script used for read excel output from ABI ViiA7 qPCR machine.
Extract Data for further computing
'''
import numpy as np
import pandas as pd

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
