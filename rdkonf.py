#/ Copyright 2019 Steven Shave (stevenshave@gmail.com)
#/ rdkonf - v1.04
#/ Implementation of high quality RDKit conformer generator as described in:
#/ Ebejer, Jean-Paul, Garrett M. Morris, and Charlotte M. Deane.
#/ "Freely available conformer generation methods: how good are they?."
#/ Journal of chemical information and modeling 52.5 (2012): 1146-1158.
#/ Please reference this code/repository : https://github.com/stevenshave/rdkonf
#/ in any publication.
#//////////////////////////////

from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import time
import argparse


def shouldWeKeep(keep, mol, q):
    for k in keep:
        if(AllChem.GetBestRMS(mol, mol, refId=k[1], prbId=q)<=0.35):
            #if(AllChem.GetBestRMS(mol, mol, refConfId=k[1], probeConfId=q)<=0.35):   #  Deprecated - previous argument calls to GetBestRMS
            #print "No! - ", AllChem.GetBestRMS(mol, mol, refConfId=k[1], probeConfId=q)
            return False
    return True

def rdkonf(input_file, output_file, numtokeep, removehydrogens):
    t00=time.time()
    filetype=input_file[len(input_file)-3:]
    suppl=None
    if (filetype=="sdf"):
        print("Reading SDF")
        suppl=Chem.SDMolSupplier(input_file)
    if (filetype=="smi" or filetype=="ism" or filetype=="can" or filetype=="txt"):
        print("Reading SMILES")
        suppl=Chem.SmilesMolSupplier(input_file,titleLine=False)
    if(suppl==None):
        print("Cannot operate on file of type", filetype)
        sys.exit()
    writer=Chem.SDWriter(output_file)
    counter=0
    for it in range(0,len(suppl)):
        counter+=1
        try:
            mol=suppl[it]
            t0=time.time()
            if mol:
                name=str(mol.GetProp("_Name"))
                e=[]
                keep=[]
                n=300
                mol=Chem.AddHs(mol)
                if removehydrogens==1:
                    mol=Chem.RemoveHs(mol)
                
                nrot=Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
                if(nrot<=12):
                    n=200
                    if(nrot<=7):
                        n=50
                print("Mol",counter, "of", len(suppl),"title =",str(mol.GetProp("_Name")), " Nrot = ", nrot, ", generating", numtokeep, "low E confs...")
                confIds=AllChem.EmbedMultipleConfs(mol, n)
                for confId in confIds:
                    ff=AllChem.UFFGetMoleculeForceField(mol, confId=confId)
                    AllChem.UFFOptimizeMolecule(mol, confId=confId)
                    ff.Minimize()
                    e.append(ff.CalcEnergy())
                    d=sorted(zip(e, confIds))
                for conf in d:
                    if(len(keep)<numtokeep):
                        if(shouldWeKeep(keep,mol, conf[1])):
                            keep.append(conf)
                counter=1
                for molout in keep:
                    mol.SetProp("_Name",name+"_"+str(counter))
                    writer.write(mol, confId=molout[1])
                    counter=counter+1
                print("took", time.time()-t0, "seconds.")
        except:
            print("Error reading molecule #"+ str(counter))
            print(sys.exc_info())
    writer.flush()
    writer.close()
    print("\nGenerated", numtokeep, "conformers for each of", len(suppl), "in", time.time()-t00, "seconds")

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="rdkonf conformer generation")
    parser.add_argument('inFile', help="In file (SDF, or SMI)")
    parser.add_argument('outFile', help="Out file (SDF)")
    parser.add_argument('numConformers', help="The number of conformers per molecule to make")
    parser.add_argument('removeHydrogens', help="Should molecules be written using implicit hydrogens (optional, default = 1)", nargs="?", default=1)

    args=parser.parse_args()
    rdkonf(args.inFile, args.outFile, int(args.numConformers), int(args.removeHydrogens))
