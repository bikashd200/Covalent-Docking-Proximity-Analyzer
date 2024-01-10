import pandas as pd
import sys
from rdkit import Chem

from rdkit import Chem
from rdkit.Chem import AllChem
#'C=CC(=O)N'

fname = "input_txt/"+sys.argv[1]
with open(fname) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

sdf = Chem.SDMolSupplier(lines[4])
writer = Chem.SDWriter(lines[5])
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        #print(smi)
        #print('ID1: ', ID1)
        ID = mol.GetProp('_Name')
        #print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol.SetProp("ID",(str(ID)+'_P'+str(count)))
                if count > 3:
                    continue
        else:
                count = 1
                mol.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol)

writer.close()


df = pd.read_csv(lines[3], sep=',', index_col= 1)
df = df.drop("TC", axis=1)
df = df.drop("Unnamed: 0", axis=1)
df.replace(to_replace =True, value = 1, inplace = True)
df.replace(to_replace =False, value = 0, inplace = True)
df.columns=df.columns.str.replace("(",'')
df.columns=df.columns.str.replace(")",'')
df.columns=df.columns.str.replace("',",':')
df.columns=df.columns.str.replace("'",'')
print(df.columns)
#df1 = df[df['CYS87.B: CustomHydrophobic']>0] #distance 4.0
df1 = df[df['CYS87.B: Hydrophobic']>0] #distance 4.5
print(df.shape)
print(df1.shape)
my_sdf = Chem.SDMolSupplier(lines[5], removeHs = False)
writer = Chem.SDWriter(lines[6])
for mol in my_sdf:
    if mol.GetProp('_Name') in df1.index.values.tolist():
        writer.write(mol)
writer.close()
