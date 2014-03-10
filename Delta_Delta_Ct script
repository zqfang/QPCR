#    FileName:BioNinja
#    Desc: Calculating Delta_Ct and Delta_Delta_Ct values from 
#    Quantitative real time polymerase chain reaction(qRT-PCR)
#    Author:Fang Zhuoqing, Wang Sishuo
#    Email:fangzhuoqing@sibs.ac.cn
#    HomePage:
#    Version: 1.2.1
#    LastChange:
#    History:2014/03/10



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import getopt
import sys
import re
from Tkinter import *


####################################################################
def	show_help():
	print "USAGE:\tpython pandas.py <--input_file=input_file> [OPTIONS]"
	print "OPTIONS"
	print "--ic or --internal_control or --expCtrl	default: GAPDH"
	print "--ec or --experiment_control	--expCtrl	default:	\'NT ES\'"
	print "-h or --help	looking for help information"
	sys.exit()

def GUI ():
	def GUI_choose_input_file ():
		import tkFileDialog
		print "starting GUI"
		tk_cwd = os.getcwd()
		input_file = tkFileDialog.askopenfilename(initialdir = tk_cwd)
		print input_file
		return (input_file)
	
	def GUI_give_Ctrl_value ():
		def show_entry_fields():
			global e1_value, e2_value
			e1_value = e1.get()
			e2_value = e2.get()
			print("internal control: %s\nexperiment control: %s" % (e1.get(), e2.get()))
		master = Tk()
		Label(master, text="internal control").grid(row=0)
		Label(master, text="experiment control").grid(row=1)
		e1 = Entry(master)
		e2 = Entry(master)
		e1.insert(10, 'GAPDH')
		e2.insert(10, "NT ES")
		e1.grid(row=0, column=1)
		e2.grid(row=1, column=1)
		Button(master, text='Quit', command=master.quit).grid(row=3, column=0, sticky=W, pady=4)
		Button(master, text='Show', command=show_entry_fields).grid(row=3, column=1, sticky=W, pady=4)
		mainloop()
		print e1_value, e2_value
		global internal_control, experiment_control
		[internal_control, experiment_control] = [e1_value, e2_value]

	global input_file
	(input_file) = GUI_choose_input_file()
	GUI_give_Ctrl_value()
	return ()


def read_cmd_line ():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["input_file=", "ic=", "interctrl=", "internal_control=", "ec=", "expCtrl=", "experiment_control=", "help"])
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

def read_param():
	if len(sys.argv) == 1:
		GUI()
	else:
		read_cmd_line()

####################################################################
read_param()

if (not 'input_file' in dir()) or (not os.path.exists(input_file)):
	print "The input_file is not given or does not exist."
	show_help()
if not 'internal_control' in dir():
	internal_control = 'GAPDH'
if not 'experiment_control' in dir():
	experiment_control = 'NT ES'

print "input_file is " + input_file;
print "internal_control is set to " + internal_control;
print "experiment_control is set to " + experiment_control;

####################################################################
# read data into a dataFrame
data = pd.read_csv(input_file)
print("The first 10 row in your original data is: ")
print(data.head(10))


# calculate Ct mean values for each replicates

# data2 = data.groupby(['SampleName','DetectorName']).mean()

# instead of using groupby, you can remove duplicates using drop_duplicates
data0 = data.drop_duplicates(['Sample Name','Detector Name']) 
data2 = data0.set_index(['Sample Name','Detector Name'])

print("The first 10 row in pre-processed(duplicates and NaN are filtered out) data is: ")
print(data2.head(10))


#Assign your sample name and detector name like this
#sample = ['NT ES','NT MESEN','NT NE2','NT NE4','NT NE6','NT TROPH','RNAi ES','RNAi MESEN','RNAi NE2','RNAi NE4','RNAi NE6','RNAi TROPH'] 

# or using set() to build an unordered collection of unique name of sample. However, set object is not itreatable, so we can convert set to list

sample0 = set(data['Sample Name'])
sample = list(sample0)

print("Your Samples are: ")
print(sample)
  
#calculate Delta_Ct value, for better interpration,  I use 'GAPDH' as a internal control as demo.

DelCt=pd.DataFrame()
for i in range(len(sample)):
    deltaCt =data2.loc[sample[i]]- data2.loc[(sample[i],internal_control)]
    deltaCt['Sample'] = sample[i]
    DelCt = DelCt.append(deltaCt)
       

#reshape your dataFrame,export to a csv file
DelCt.index.name = 'Detector'
DelCt2 = DelCt.reset_index()
DelCt3 = DelCt2.set_index(['Sample','Detector'])

# If you want to rename column names, using this code 
#DelCt3.rename(columns={'Ct': 'Delta Ct', 'Ct Mean': 'Delta CtMean'}, inplace=True)
DelCt3.rename(columns={'Ct Mean': 'DeltaCt'}, inplace=True)

DelCt3.to_csv('Delta_Ct.csv')
print("The first 10 row in your Delta_Ct value is: ")
print(DelCt3.head(10))

#calculate Delta_Delta_Ct, for demo, I use 'NT ES' as our experiment control group.

DDelCt =pd.DataFrame()
for i in range(len(sample)):
    DDeltaCt = DelCt3.loc[sample[i]]-DelCt3.loc[experiment_control]
    DDeltaCt['Sample'] = sample[i]
    DDelCt = DDelCt.append(DDeltaCt)
    



#reshape your dataFrame,export to a csv file
DDelCt.index.name = 'Detector'
DDelCt2 = DDelCt.reset_index()
DDelCt3 = DDelCt2.set_index(['Sample','Detector'])
DDelCt3.rename(columns={'DeltaCt': 'DDeltaCt'}, inplace=True)
DDelCt4=DDelCt3.dropna(how='all')
DDelCt4.to_csv('Delta_Delta_Ct.csv')

print("The first 10 row in your Delta_Delta_Ct values is: ")
print(DDelCt3.head(10))

# calculate FoldChange, and export to a csv file

foldChange0 = pow(2,-DDelCt4)
foldChange0.rename(columns={'DDeltaCt':'FoldChange'},inplace=True)
# foldChange = foldChange0.drop(experiment_control,level=0)
foldChange1 = foldChange0.drop(internal_control,level=1)
foldChange1.to_csv('FoldChange.csv')

print("The first 10 row in your FoldChange values is: ")
print(foldChange1.head(10))

    
#reshape your final results
#extract columns you needed

d1=data2[['Ct Mean']]
d2=DelCt3[['DeltaCt']]
d3=DDelCt4[['DDeltaCt']]
d4=foldChange0[['FoldChange']]

#merge data,export to a csv file.

m1=pd.merge(d1,d2,left_index=True,right_index=True)
m2=pd.merge(m1,d3,left_index=True,right_index=True)
final_report=pd.merge(m2,d4,left_index=True,right_index=True)
final_report.index.names=['Sample','Detector']
final_report.to_csv('FinalResults.csv')

print('The first 10 rows in your Final results: ')
print(final_report.head(10))
