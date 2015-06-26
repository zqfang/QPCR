




####################################################################
# read data into a dataFrame
data = pd.read_csv(input_file)
print("The first 10 row in your original data is: ")
print(data.head(5))


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
print(DelCt3.head(5))

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
print(DDelCt3.head(5))

# calculate FoldChange, and export to a csv file

foldChange0 = pow(2,-DDelCt4)
foldChange0.rename(columns={'DDeltaCt':'FoldChange'},inplace=True)
# foldChange = foldChange0.drop(experiment_control,level=0)
foldChange1 = foldChange0.drop(internal_control,level=1)
foldChange1.to_csv('FoldChange.csv')

print("The first 10 row in your FoldChange values is: ")
print(foldChange1.head(5))

    
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
print(final_report.head(5))
