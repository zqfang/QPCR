qPCR
====


     
This python script aims at calculating the Delta_Ct, Delta_Delta_Ct and FoldChange value which
produced by Quantitative real time polymerase chain reaction.
    
Pandas, a python data analysis tool, is required for using this script. I use Ipython as my interpreter.
    
  
    
The file path and file format mannually need to be specified.
The following file formats are supported: csv,txt et. al.
    
In addition, your should specify internal control name and expriment control name for your own data sets.
    
Last but not least, your data set struture needs to be the same as my trainning data set.That's to say, 
column names must be 'Sample Name','Detector Name','Ct','Ct Mean'. But 'Ct StdEV' is optional.
