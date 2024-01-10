from rdkit import Chem
from rdkit.Chem import AllChem

import sys
import os

fname = "input_txt/"+sys.argv[1]
with open(fname) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

#filename_protein = lines[0]

#sdf = Chem.SDMolSupplier('covalentLib_150eachWarhead_3D_1.sdf')
#writer = Chem.SDWriter('subtract.sdf')

sdf = Chem.SDMolSupplier(lines[4])
writer = Chem.SDWriter(lines[1])

# Nitrofurylsulphonyl
core = Chem.MolFromSmiles('[O-][N+](=O)c1ccc(S)o1')
sub_core = Chem.MolFromSmiles('S')
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        ID = mol.GetProp('_Name')
        print('ID: ', ID) 
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #print(c1h1.GetProp('_Name'))
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        #print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)




# Acrylamides
#core = Chem.MolFromSmiles('[#7]-[#6](=[#8])-[#6]=[#6]')
core = Chem.MolFromSmiles('NC(=O)C=C')
sub_core = Chem.MolFromSmiles('C=C')
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        ID = mol.GetProp('_Name')
        print('ID: ', ID)
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #print(c1h1.GetProp('_Name'))
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        #print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)







# Alkynes
core = Chem.MolFromSmiles('C#C')
sub_core = Chem.MolFromSmiles('C#C')
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)
#Allylester
core = Chem.MolFromSmiles('C=CC(O)=O')
sub_core = Chem.MolFromSmiles('C=C')
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

# Azaridine

core = Chem.MolFromSmiles('C1N(C)C1')
sub_core = Chem.MolFromSmiles('CC')
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

#Disulfide
core = Chem.MolFromSmiles('SS')
sub_core = Chem.MolFromSmiles('SS')
ID1 = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)


# Ketoamide
core = Chem.MolFromSmiles('O=C(C=O)N')
sub_core = Chem.MolFromSmiles('O=C(C=O)')
ID1 = []
ID_list = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        ID_list.append([ID])
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

# Ketohalogen, chloroacetamide

core = Chem.MolFromSmarts('[#7]-[#6](=[#8])-[#6]-[Cl,Br,I,F]')
sub_core = Chem.MolFromSmarts('[#6]-[Cl,Br,I,F]')
#core = Chem.MolFromSmiles('C(C(=O)N)X')
#sub_core = Chem.MolFromSmiles('CX')

ID1 = []
ID_list = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        ID_list.append([ID])
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

# Sulfonylallyl, VinylSulfones

#core = Chem.MolFromSmarts('[#7]-[#6](=[#8])-[#6]-[Cl,Br,I,F]')
#sub_core = Chem.MolFromSmarts('[#6]-[Cl,Br,I,F]')
core = Chem.MolFromSmiles('S(C=C)(=O)=O')
sub_core = Chem.MolFromSmiles('C=C')

ID1 = []
ID_list = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        ID_list.append([ID])
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

# Thiol

core = Chem.MolFromSmarts('[S;H]')
#sub_core = Chem.MolFromSmarts('[S;X1]')
#core = Chem.MolFromSmiles('SH')
sub_core = Chem.MolFromSmiles('S')

ID1 = []
ID_list = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        ID_list.append([ID])
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

# Cyanoacrylates

#core = Chem.MolFromSmarts('[S;H]')
#sub_core = Chem.MolFromSmarts('[S;X1]')
core = Chem.MolFromSmiles('NC(=O)C(=C)C#N')
sub_core = Chem.MolFromSmiles('C=C')

ID1 = []
ID_list = []
count = 1
for mol in sdf:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        tmp = Chem.ReplaceSidechains(mol, core)
        if tmp is None:
            continue
        #print(smi)
        du = Chem.MolFromSmiles('*')
        c1h=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        tmp = Chem.ReplaceSidechains(c1h, sub_core)
        if tmp is None:
            continue
        du = Chem.MolFromSmiles('*')
        mol1=AllChem.ReplaceSubstructs(tmp,du,Chem.MolFromSmiles('[H]'),True)[0]
        #mol1 = Chem.RemoveHs(c1h1)
        #print('ID1: ', ID1)
        ID = mol1.GetProp('_Name')
        print('ID: ', ID)
        ID_list.append([ID])
        if ID == ID1:
                count+=1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        else:
                count = 1
                mol1.SetProp("_Name",(str(ID)+'_P'+str(count)))
                mol1.SetProp("ID",(str(ID)+'_P'+str(count)))
        ID1 = ID
        writer.write(mol1)

writer.close()
