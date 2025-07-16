import numpy as np
import os
from rdkit import Chem
import faiss, joblib, pathlib
from src.utils import mol_from_smiles, fp_array

def parse_smi_file(smi_path):
    """Parse SMI file and return list of (smiles, name) tuples"""
    data = []
    with open(smi_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            # Split on whitespace - first part is SMILES, rest is compound name
            parts = line.split(None, 1)  # Split on any whitespace, max 1 split
            if len(parts) >= 2:
                smiles, name = parts[0], parts[1]
                data.append({'smiles': smiles, 'drug_name': name})
            elif len(parts) == 1:
                # Only SMILES, no name
                smiles = parts[0]
                data.append({'smiles': smiles, 'drug_name': f'Compound_{line_num}'})
            else:
                print(f"Warning: Skipping empty line {line_num}")
    
    return data

def build_ecfp6_index(smi_path: str, output_dir: str = "artifacts"):
    """Build ECFP6 similarity index from SMI file"""
    try:
        print(f"Reading data from: {smi_path}")
        data = parse_smi_file(smi_path)
        
        if not data:
            raise ValueError("No valid data found in SMI file")
        
        # Validate and filter SMILES
        valid_data = []
        valid_fps = []
        
        print("Processing molecules...")
        for i, entry in enumerate(data):
            try:
                mol = mol_from_smiles(entry['smiles'])
                fp = fp_array(mol)
                valid_data.append(entry)
                valid_fps.append(fp)
            except Exception as e:
                print(f"Warning: Skipping invalid SMILES at line {i+1}: {entry['smiles']} ({e})")
        
        if not valid_fps:
            raise ValueError("No valid molecules found in the input file")
            
        fps = np.vstack(valid_fps)
        faiss.normalize_L2(fps)

        print(f"Building index for {len(valid_fps)} molecules...")
        index = faiss.IndexFlatIP(fps.shape[1])
        index.add(fps)

        pathlib.Path(output_dir).mkdir(exist_ok=True)
        
        faiss.write_index(index, f"{output_dir}/ecfp6.index")
        joblib.dump({"smiles": [d['smiles'] for d in valid_data], "meta": valid_data}, f"{output_dir}/meta.pkl")
        
        print(f"✓ ECFP6 Index built successfully!")
        print(f"  - {len(valid_fps)} molecules indexed")
        print(f"  - Index saved to: {output_dir}/ecfp6.index")
        print(f"  - Metadata saved to: {output_dir}/meta.pkl")
        
    except Exception as e:
        print(f"Error building ECFP6 index: {e}")
        raise

def build_mol2vec_index(smi_path: str, model_path: str, output_dir: str = "artifacts"):
    """Build Mol2Vec similarity index from SMI file and model"""
    try:
        print(f"Building Mol2Vec index from: {smi_path}")
        print(f"Using model: {model_path}")
        
        if not os.path.exists(model_path):
            raise ValueError(f"Mol2Vec model not found: {model_path}")
            
        from mol2vec.features import mol2alt_sentence, sentences2vec, DfVec
        from gensim.models import word2vec

        # Load model
        model = word2vec.Word2Vec.load(model_path)

        data = parse_smi_file(smi_path)
        if not data:
            raise ValueError("No valid data found in SMI file")
        print("Processing molecules with Mol2Vec...")
        valid_embeddings = []
        valid_indices = []
        
        for i, entry in enumerate(data):
            try:
                mol = Chem.MolFromSmiles(entry['smiles'])
                if mol is None:
                    continue
                sentence = mol2alt_sentence(mol, 1)
                embedding = DfVec(sentences2vec([sentence], model)).vec
                
                # Flatten embedding to ensure it's 1D (in case it's (1, 300) instead of (300,))
                embedding = embedding.flatten()
                
                # Debug embedding dimensions
                if len(valid_embeddings) < 5:
                    print(f"  Line {i+1}: {entry['smiles'][:50]}... -> embedding shape: {embedding.shape}")
                
                # Only accept 300-dimensional embeddings
                if embedding.shape == (300,):
                    valid_embeddings.append(embedding)
                    valid_indices.append(i)
                else:
                    print(f"Warning: Skipping SMILES at line {i+1} due to unexpected embedding shape: {embedding.shape}")
            except Exception as e:
                print(f"Warning: Skipping SMILES at line {i+1}: {entry['smiles']} ({e})")
        
        if not valid_embeddings:
            raise ValueError("No valid embeddings generated")
            
        emb = np.vstack(valid_embeddings).astype("float32")
        faiss.normalize_L2(emb)

        print(f"Building Mol2Vec index for {len(valid_embeddings)} molecules...")
        index = faiss.IndexFlatIP(emb.shape[1])
        index.add(emb)
        
        pathlib.Path(output_dir).mkdir(exist_ok=True)
        faiss.write_index(index, f"{output_dir}/mol2vec.index")
        
        print(f"✓ Mol2Vec Index built successfully!")
        print(f"  - {len(valid_embeddings)} molecules indexed")
        print(f"  - Index saved to: {output_dir}/mol2vec.index")
        
    except Exception as e:
        print(f"Error building Mol2Vec index: {e}")
        raise
