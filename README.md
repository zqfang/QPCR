# qPCR 

### Process You qRT-PCR Results Easy  
====


     
This python script aims at calculating the **Delta_Ct, Delta_Delta_Ct and FoldChange** value which
produced by Quantitative real time polymerase chain reaction.
    
**Pandas**, a python data analysis tool, is required for using this script.
    
The following file formats are supported: **xls, xlsx, csv, txt**. 
    
In addition, your should specify **reference control** name and **expriment control** name for your own data sets.
    
For ABi 7900 users, data column names must be 'Sample Name','Detector Name','Ct','Ct Mean'. But 'Ct StdEV' is optional.

For ABi ViiA 7 users, you can use qpcrRead.py to extract data computing results directly.

##Useage

Python2.7 or Python3 users

#####Before use this module, see help
python qpcrRead.py -h 
python Delta_Delta_Ct.py -h

#####Extract Data from Original Data output
e.g. python qpcrRead.py -f foo.xls --header 35 --tail 5  -o 20150625_NPC_Knockdown

#####Calculate Detal_Ct, Delta_Delta_Ct, Fold_Change and generate output file
e.g. python Delta_Delta_Ct.py -f foo.csv -r GAPDH -c hESC -o 20150625_NPC_Knockdown

###Generate Plots: line, bar plot

##To Do List

1. process data automatically without pre-filter outliner mannully.
2. Generate Plots using Matplotlib automatically 
3. Generate Matplotlib Plotting Scripts for customized modification

