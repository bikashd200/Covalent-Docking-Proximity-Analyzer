from rdkit import Chem
from rdkit.Chem import AllChem
#'C=CC(=O)N'
sdf = Chem.SDMolSupplier('ube2n_Enamine_wh_40K.sdf')
writer = Chem.SDWriter('ube2n_Enamine_wh_40K_rename.sdf')
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
        else:
                count = 1
                mol.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol)

writer.close()
