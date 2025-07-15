import numpy as np, joblib, faiss, os
from rdkit import Chem
from mol2vec.features import mol2alt_sentence, sentences2vec, DfVec
from gensim.models import word2vec
from src.utils import fp_array

def search(query_smi, top_k, model_path=None):
    try:
        # Validate SMILES
        mol = Chem.MolFromSmiles(query_smi)
        if mol is None:
            print(f"Error: Invalid SMILES string: {query_smi}")
            return
        
        # Check if index files exist
        if not os.path.exists("artifacts/ecfp6.index"):
            print("Error: ECFP6 index not found. Please run build_index.py first.")
            return
        if not os.path.exists("artifacts/meta.pkl"):
            print("Error: Metadata file not found. Please run build_index.py first.")
            return
            
        index = faiss.read_index("artifacts/ecfp6.index")
        data = joblib.load("artifacts/meta.pkl")
        vec = fp_array(Chem.MolFromSmiles(query_smi))
        faiss.normalize_L2(vec.reshape(1, -1))
        D, I = index.search(vec.reshape(1, -1), top_k)

        # Output ECFP6 hits
        results = []
        for rank, (dist, idx) in enumerate(zip(D[0], I[0]), 1):
            hit = data["meta"][idx]
            results.append((rank, hit["drug_name"], hit["smiles"], float(dist)))
        
        if model_path:
            if not os.path.exists(model_path):
                print(f"Error: Mol2Vec model not found at {model_path}")
                return
            if not os.path.exists("artifacts/mol2vec.index"):
                print("Error: Mol2Vec index not found. Please build it first.")
                return
                
            model = word2vec.Word2Vec.load(model_path)
            vec2 = DfVec(sentences2vec([mol2alt_sentence(Chem.MolFromSmiles(query_smi), 1)], model)).vec.astype("float32")
            faiss.normalize_L2(vec2.reshape(1, -1))
            index2 = faiss.read_index("artifacts/mol2vec.index")
            D2, I2 = index2.search(vec2.reshape(1, -1), top_k)

            rank_sum = {}
            for r, idx in enumerate(I[0]): rank_sum[idx] = r
            for r, idx in enumerate(I2[0]): rank_sum[idx] = rank_sum.get(idx, 0) + r

            sorted_hits = sorted(rank_sum.items(), key=lambda x: x[1])[:top_k]
            print("✓ Blended ECFP6 + Mol2Vec:")
            for i, (idx, _) in enumerate(sorted_hits, 1):
                hit = data["meta"][idx]
                print(f"{i:>2}. {hit['drug_name']}  |  {hit['smiles']}")
        else:
            for rank, name, smi, tanimoto in results:
                print(f"{rank:>2}. {name}  |  {smi}  |  Tanimoto: {tanimoto:.3f}")
    
    except Exception as e:
        print(f"Error during search: {e}")
        return
