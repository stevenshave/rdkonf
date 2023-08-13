#!/usr/bin/env python

"""
Conformer generation using RDKit's ETKDGv3 embedding parameters
"""
# //////////////////////////////
# / THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# / EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# / MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# / IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# / OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# / ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# / OTHER DEALINGS IN THE SOFTWARE.

# / You may not distribute copies of the programs source code in any medium.
# / You may not modify the program or any portion of it, thus forming a work
# / based on the program.
# /-------------------------------------------------------------------------
# / Copyright 2023 Steven Shave (stevenshave@gmail.com)
# / rdkonf_etkdgv3 - v1.00
# //////////////////////////////
from pathlib import Path
import fire
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from typing import Union
import sys
import time
from rdkit.Chem.MolStandardize import rdMolStandardize
import multiprocessing as mp

def shouldWeKeep(keep:list, mol, q:float):
    for k in keep:
        if AllChem.GetBestRMS(mol, mol, refId=k[1], prbId=q) <= 0.35:
            return False
    return True

def _init_worker_confgen(num_conformers_:int, all_smiles_:list[str]):
    global shared_num_conformers, shared_all_smiles
    shared_num_conformers = num_conformers_
    shared_all_smiles = all_smiles_
    
def _work_get_3D_mol(line_number:int, add_rdkit_cansmi_key:bool=True):
    global shared_num_conformers
    line=shared_all_smiles[line_number]
    split_results=line.split()
    smiles=split_results[0]
    if len(split_results)==1:
        id=str(line_number)
    else:
        id=" ".join(split_results[1:])
        
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
        except:
            print(f"Creating mol from {smiles} failed, continuing")
            return (None, None, None, line_number)
    if not mol:
        print("Invalid mol from line:", line)
        return (None, None, None, line_number)
    try:
        mol=rdMolStandardize.Cleanup(mol)
    except:
        print(f"Failed to sanitise {smiles}")
    try:
        if add_rdkit_cansmi_key:
            canonical_smiles = Chem.MolToSmiles(mol)
            mol.SetProp("RDKit_Canonical_SMILES", canonical_smiles)
        mol.SetProp("_Name", id)
        mol_name=id
        
        t0 = time.time()
        if mol:
            e = []
            keep = []
            n = 300
            nrot = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
            if nrot <= 12:
                n = 200
                if nrot <= 7:
                    n = 50
            
            mol=Chem.AddHs(mol)
            Chem.rdmolops.AssignStereochemistryFrom3D(mol)
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
            params = AllChem.ETKDGv3()
            try:
                confIds = AllChem.EmbedMultipleConfs(mol, n, params)
            except:
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                try:
                    confIds = AllChem.EmbedMultipleConfs(mol, n, params)
                except:
                    print(f"Conformer embedding failed for line:", line)
                    return (None, None, None, line_number)

            for confId in confIds:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
                AllChem.UFFOptimizeMolecule(mol, confId=confId)
                e.append(ff.CalcEnergy())
            d = sorted(zip(e, confIds))
            for conf in d:
                if len(keep) < shared_num_conformers:
                    if shouldWeKeep(keep, mol, conf[1]):
                        keep.append(conf)
            return (mol, keep, mol_name, line_number)
    except:
        print("Error making a molecule from:\n" + str(line))
        print(sys.exc_info())
        return (None, None, None, line_number)


def smi_to_sdf(input_file: Path, num_conformers:int, n_jobs:Union[None, int]=None):
    if n_jobs is None:
        n_jobs = mp.cpu_count()
    output_file = str(input_file) + ".sdf"
    t00 = time.time()
    writer = Chem.SDWriter(str(output_file))
    problematic_lines=[]
    with open(str(input_file), mode="rb") as f:
        num_lines_in_file = sum(1 for _ in f)
    with mp.Pool(n_jobs, initializer=_init_worker_confgen, initargs=(num_conformers,[line.strip() for line in open(str(input_file)).readlines() if len(line)>3])) as pool:
            for mol, keep, mol_name, line_number in tqdm(pool.imap_unordered(_work_get_3D_mol, smiles_line_numbers:=range(num_lines_in_file), chunksize=3), total=len(smiles_line_numbers)):
                if mol is None:
                    problematic_lines.append(line_number)
                    continue
                for molout_i, molout in enumerate(keep):
                    mol.SetProp("_Name", f"{mol_name}_{molout_i+1}")
                    writer.write(mol, confId=molout[1])
    writer.close()
    print("\nProcessed", num_lines_in_file, "in", time.time() - t00, "seconds")


if __name__ == "__main__":
    fire.Fire(smi_to_sdf)
    