import MDAnalysis as mda
import prolif as plf
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import PandasTools
import sys
import os
import glob


fname = "input_txt/"+sys.argv[1]
with open(fname) as file:
	lines = file.readlines()
	lines = [line.rstrip() for line in lines]

# load protein
#prot = mda.Universe("/home/bdebnath_championsoncology_com/SiFt/rec.crg.pdb")
#prot = plf.Molecule.from_mda(prot)
#prot.n_residues

filename_protein = lines[0]
#filename_protein = "/home/bdebnath_championsoncology_com/MK14/test/2qd9-dock-prep-min.pdb"
print('step1')
u = mda.Universe(filename_protein, guess_bonds=True, vdwradii={"H": 0.98})
print('step2')
elements = mda.topology.guessers.guess_types(u.atoms.names)
print('step3')
u.add_TopologyAttr("elements", elements)
print('step4')
#prot_sel = u.select_atoms("not resname LIG and not resname HOH and not resname DMS")
print('step5')
prot = plf.Molecule.from_mda(u)
print('step6')

path = lines[1]
#path = str("/home/bdebnath_championsoncology_com/MK14/test/ligands-hH.sdf")

df_frame = PandasTools.LoadSDF(path, includeFingerprints=False)
#frame = PandasTools.LoadSDF(my_sdf_file, includeFingerprints=False)
df_ID = df_frame.loc[:,['ID']]
del(df_frame)
#df_ID = df_ID[~df_ID['ID'].duplicated(keep='first')]
print ('Total poses = ' + str(len(df_ID)))

# load ligands
#path = str("/home/bdebnath_championsoncology_com/MK14/test/ligands-hH.sdf")

#path = str('C:/Users/Bikash Debnath/Documents/prolif/MK14/ligands-hH-cleaned.sdf')
#df_ID['docking_score'] = df_ID['REMARK'].str[18:23]
df_ID = df_ID.reset_index(drop=True)

#path = plf.datafiles.datapath / "vina" / "vina_output.sdf"
#lig_suppl = list(plf.mol2_supplier(path))

lig_suppl = list(plf.sdf_supplier(path))
print('step7')
# generate fingerprint
## distance is set at 4 angstrom
class CustomHydrophobic(plf.interactions.Hydrophobic):
    def __init__(self):
        super().__init__(distance=4.0)

#class CustomHydrophobic_4(plf.interactions.Hydrophobic):
#    def __init__(self):
#        super().__init__(distance=4.0)

fp = plf.Fingerprint(["Hydrophobic", "CustomHydrophobic"])
#fp = plf.Fingerprint()
from tqdm.auto import tqdm
progress=True
iterator = tqdm(lig_suppl) if progress else lig_suppl
del(lig_suppl)
ifp = []
for i, lig_mol in enumerate(iterator):
	data = fp.generate(lig_mol, prot, residues=['CYS87.B'], return_atoms=True)  ## change residue number as per your protein
	data["Frame"] = i
	ifp.append(data)

print('step8')
#fp.run_from_iterable(lig_suppl, prot)
print('step9')
fp.ifp = ifp
del(ifp)
df = fp.to_dataframe()
df.columns = df.columns.droplevel(0)
print(df.head())

#path = str('C:/Users/Bikash Debnath/Documents/prolif/MK14/ligands-hH-cleaned.sdf')
#print(df_ID.head())
df_2 = [df_ID, df]
df_3 = pd.concat(df_2, axis=1)
del(df_ID)
del(df_2)
# load the reference

path = lines[2]
#path = str("/home/bdebnath_championsoncology_com/MK14/test/ref.sdf")
#ref = mda.Universe("/home/bdebnath_championsoncology_com/topandbottom20poses/top/ref.sdf")
lig_suppl = list(plf.sdf_supplier(path))
#ref = plf.Molecule.from_mda(ref)
# generate IFP
print('step10')

fp.run_from_iterable(lig_suppl, prot)
print('step11')
del(lig_suppl)
df0 = fp.to_dataframe()
del(fp)
print(df0.head())
df0.rename({0: "ref"}, inplace=True)
# drop the ligand level on both dataframes
df0.columns = df0.columns.droplevel(0)

print(df.head())
# concatenate them
df = (pd.concat([df0, df])
        .fillna(False)
        .sort_index(axis=1, level=0,
                    key=lambda index: [plf.ResidueId.from_string(x) for x in index]))
del(df0)
print(df.head())
#df.to_csv("MK14-docking-prolif_out.csv")

from rdkit import DataStructs

bvs = plf.to_bitvectors(df)
del(df)
#for index,row in df.iterrows():
#   df.loc[index,'d'] = np.random.randint(0, 10)
path = lines[3]
for i, bv in enumerate(bvs[1:]):
    #print(i)
    tc = DataStructs.TanimotoSimilarity(bvs[0], bv)
    df_3.loc[i,'TC'] = tc
    #print(f"{i}: {tc:.3f}")
#df_3.columns = df_3.columns.droplevel(0)
del(bvs)
df_3.to_csv(path)
