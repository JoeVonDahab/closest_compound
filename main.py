#!/usr/bin/env python3
"""
Closest Compound - Molecular Similarity Search Tool

This module provides a command-line interface for building molecular similarity
indices and performing similarity searches using ECFP6 fingerprints and 
optionally Mol2Vec embeddings.
"""

import sys
import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Molecular Similarity Search Tool")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Build index command
    build_parser = subparsers.add_parser('build', help='Build similarity index')
    build_parser.add_argument('input_file', help='Input file (.smi or .csv format)')
    build_parser.add_argument('--output-dir', default='artifacts', help='Output directory')
    
    # Search command
    search_parser = subparsers.add_parser('search', help='Search for similar compounds')
    search_parser.add_argument('smiles', help='Query SMILES string')
    search_parser.add_argument('-k', '--top-k', type=int, default=5, help='Number of results to return')
    search_parser.add_argument('--mol2vec-model', help='Path to Mol2Vec model for blended search')
    
    args = parser.parse_args()
    
    if args.command == 'build':
        from src.index_builder import build_ecfp6_index
        from scripts.build_index import parse_smiles_file
        import os
        
        input_path = args.input_file
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
                build_ecfp6_index(temp_csv, args.output_dir)
                os.remove(temp_csv)  # Clean up
            else:
                build_ecfp6_index(input_path, args.output_dir)
        except Exception as e:
            print(f"Error building index: {e}")
            sys.exit(1)
            
    elif args.command == 'search':
        from src.query import search
        search(args.smiles, args.top_k, args.mol2vec_model)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
