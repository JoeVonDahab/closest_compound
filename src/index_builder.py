import pandas as pd
import numpy as np
from rdkit import Chem
import faiss, joblib, pathlib
from src.utils import mol_from_smiles, fp_array

def build_ecfp6_index(csv_path: str, output_dir: str = "artifacts"):
    """Build ECFP6 similarity index from CSV file"""
    try:
        print(f"Reading data from: {csv_path}")
        df = pd.read_csv(csv_path)
        
        if 'smiles' not in df.columns:
            raise ValueError("CSV file must contain a 'smiles' column")
        
        smiles = df["smiles"].tolist()
        
        # Validate and filter SMILES
        valid_data = []
        valid_fps = []
        
        print("Processing molecules...")
        for i, smi in enumerate(smiles):
            try:
                mol = mol_from_smiles(smi)
                fp = fp_array(mol)
                valid_data.append(i)
                valid_fps.append(fp)
            except Exception as e:
                print(f"Warning: Skipping invalid SMILES at row {i}: {smi} ({e})")
        
        if not valid_fps:
            raise ValueError("No valid molecules found in the input file")
            
        fps = np.vstack(valid_fps)
        faiss.normalize_L2(fps)

        print(f"Building index for {len(valid_fps)} molecules...")
        index = faiss.IndexFlatIP(fps.shape[1])
        index.add(fps)

        pathlib.Path(output_dir).mkdir(exist_ok=True)
        
        # Save only valid data
        valid_df = df.iloc[valid_data].copy()
        valid_smiles = valid_df["smiles"].tolist()
        
        faiss.write_index(index, f"{output_dir}/ecfp6.index")
        joblib.dump({"smiles": valid_smiles, "meta": valid_df.to_dict("records")}, f"{output_dir}/meta.pkl")
        
        print(f"✓ ECFP6 Index built successfully!")
        print(f"  - {len(valid_fps)} molecules indexed")
        print(f"  - Index saved to: {output_dir}/ecfp6.index")
        print(f"  - Metadata saved to: {output_dir}/meta.pkl")
        
    except Exception as e:
        print(f"Error building ECFP6 index: {e}")
        raise

def build_mol2vec_index(smiles, model_path, output_path):
    from mol2vec.features import mol2alt_sentence, sentences2vec, DfVec
    from gensim.models import word2vec

    model = word2vec.Word2Vec.load(model_path)
    emb = np.vstack([
        DfVec(sentences2vec([mol2alt_sentence(Chem.MolFromSmiles(s), 1)], model)).vec
        for s in smiles
    ]).astype("float32")
    faiss.normalize_L2(emb)

    index = faiss.IndexFlatIP(emb.shape[1])
    index.add(emb)
    faiss.write_index(index, output_path)
    print("✓ Mol2Vec Index built.")
