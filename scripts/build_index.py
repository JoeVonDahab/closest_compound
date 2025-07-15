from src.index_builder import build_ecfp6_index
import sys, os
import pandas as pd

def parse_smiles_file(file_path):
    """Parse SMILES file (.smi format) and return DataFrame-like structure"""
    data = []
    with open(file_path, 'r') as f:
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
    
    return pd.DataFrame(data)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_path = sys.argv[1]
    else:
        # Check for default files
        if os.path.exists("library/library.smi"):
            input_path = "library/library.smi"
        elif os.path.exists("data/library.csv"):
            input_path = "data/library.csv"
        else:
            print("Error: No input file specified and no default file found.")
            print("Usage: uv run build_index.py <input_file>")
            print("Supported formats: .csv (with 'smiles' column) or .smi (SMILES Name)")
            sys.exit(1)
    
    if not os.path.exists(input_path):
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)
    
    print(f"Building index from: {input_path}")
    
    try:
        if input_path.endswith('.smi'):
            # Convert SMI to temporary CSV for processing
            df = parse_smiles_file(input_path)
            temp_csv = "temp_library.csv"
            df.to_csv(temp_csv, index=False)
            build_ecfp6_index(temp_csv)
            os.remove(temp_csv)  # Clean up
        else:
            build_ecfp6_index(input_path)
    except Exception as e:
        print(f"Error building index: {e}")
        sys.exit(1)
