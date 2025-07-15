import sys
from src.query import search
from rdkit import Chem

def validate_smiles(smiles_str):
    """Validate SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        return mol is not None
    except:
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: nearest.py <query_smiles> [top_k] [mol2vec_model_path]")
        print("Example: python nearest.py 'CCO' 10")
        print("Example with Mol2Vec: python nearest.py 'CCO' 10 /path/to/model.pkl")
        sys.exit(1)
    
    query_smi = sys.argv[1]
    
    # Validate SMILES
    if not validate_smiles(query_smi):
        print(f"Error: Invalid SMILES string: '{query_smi}'")
        print("Please provide a valid SMILES string")
        sys.exit(1)
    
    # Parse top_k parameter
    top_k = 5  # default
    if len(sys.argv) > 2:
        try:
            top_k = int(sys.argv[2])
            if top_k <= 0:
                raise ValueError
        except ValueError:
            print(f"Error: top_k must be a positive integer, got: {sys.argv[2]}")
            sys.exit(1)
    
    # Parse model path parameter
    model_path = None
    if len(sys.argv) > 3:
        model_path = sys.argv[3]
        import os
        if not os.path.exists(model_path):
            print(f"Warning: Mol2Vec model file not found: {model_path}")
            print("Proceeding with ECFP6 only...")
            model_path = None
    
    print(f"Searching for compounds similar to: {query_smi}")
    print(f"Returning top {top_k} results")
    if model_path:
        print(f"Using Mol2Vec model: {model_path}")
    print("-" * 50)
    
    search(query_smi, top_k, model_path)
