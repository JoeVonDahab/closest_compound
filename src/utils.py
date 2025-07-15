import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors

DIM, RADIUS = 2048, 3

def mol_from_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    return mol

def fp_array(mol):
    # Use the modern MorganGenerator instead of deprecated GetMorganFingerprintAsBitVect
    generator = rdMolDescriptors.GetMorganGenerator(radius=RADIUS, fpSize=DIM)
    fp = generator.GetFingerprint(mol)
    arr = np.zeros((DIM,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype("float32")