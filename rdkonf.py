#!/usr/bin/env python

"""
Implementation of high quality RDKit conformer generator as described in:
Ebejer, Jean-Paul, Garrett M. Morris, and Charlotte M. Deane. 'Freely available
conformer generation methods: how good are they?.' Journal of chemical
information and modeling 52.5 (2012): 1146-1158.

"""

import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import time


def should_we_keep_a_mol(keep_list, mol, q):
    for k in keep_list:
        if AllChem.GetBestRMS(mol, mol, refId=k[1], prbId=q) <= 0.35:
            return False
    return True


def smiles_to_3dmol(smiles: str, title: str = ""):

    splitsmiles = smiles.split()
    if len(splitsmiles) > 1:
        smiles = splitsmiles[0]
        title = "".join(splitsmiles[1:])
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        e = []
        keep = []
        n = 300
        mol = Chem.AddHs(mol)
        mol = Chem.RemoveHs(mol)
        nrot = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
        if nrot <= 12:
            n = 200
            if nrot <= 7:
                n = 50
        confIds = AllChem.EmbedMultipleConfs(mol, n)
        for confId in confIds:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
            AllChem.UFFOptimizeMolecule(mol, confId=confId)
            e.append(ff.CalcEnergy())
            d = sorted(zip(e, confIds))
        for conf in d:
            if len(keep) < 1:
                if should_we_keep_a_mol(keep, mol, conf[1]):
                    keep.append(conf)
        if len(keep) == 1:
            mol.SetProp("_Name", title)
            return mol
    return None


def rdkonf_on_sdf_file(input_file: Path):
    output_file = f"{input_file}.sdf"
    numtokeep = 1
    removehydrogens = 1
    t00 = time.time()
    writer = Chem.SDWriter(output_file)
    mol_counter = 0
    smiles_f = open(str(input_file))
    line = smiles_f.readline().strip()
    while line is not None and len(line) > 3:
        try:
            smiles, ids = line.split()
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", ids)
            mol_counter += 1
            t0 = time.time()
            if mol:
                name = mol.GetProp("_Name")
                e = []
                keep = []
                n = 300
                mol = Chem.AddHs(mol)
                if removehydrogens == 1:
                    mol = Chem.RemoveHs(mol)

                nrot = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
                if nrot <= 12:
                    n = 200
                    if nrot <= 7:
                        n = 50
                print(
                    f"Mol {mol_counter}, title = {mol.GetProp('_Name')}, Nrot = {nrot}, generating {numtokeep} low E confs..."
                )
                confIds = AllChem.EmbedMultipleConfs(mol, n)
                for confId in confIds:
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
                    AllChem.UFFOptimizeMolecule(mol, confId=confId)
                    e.append(ff.CalcEnergy())
                    d = sorted(zip(e, confIds))
                for conf in d:
                    if len(keep) < numtokeep:
                        if should_we_keep_a_mol(keep, mol, conf[1]):
                            keep.append(conf)
                counter = 1
                for molout in keep:
                    mol.SetProp("_Name", name)
                    writer.write(mol, confId=molout[1])
                    counter = counter + 1
                print("took", time.time() - t0, "seconds.")
        except:
            print(f"Error making a molecule from:\n {line}")
            print(sys.exc_info())
        line = smiles_f.readline().strip()
    writer.close()
    print(f"\nGenerated {numtokeep} in {time.time() - t00} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_file", help="Smiles file")
    args = parser.parse_args()
    rdkonf_on_sdf_file(Path(args.smiles_file))
