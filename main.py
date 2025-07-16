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
    build_parser.add_argument('input_file', help='Input SMI file (SMILES compound_name format)')
    build_parser.add_argument('--output-dir', default='artifacts', help='Output directory')
    build_parser.add_argument('--mol2vec-model', help='Path to Mol2Vec model to also build Mol2Vec index')
    
    # Search command
    search_parser = subparsers.add_parser('search', help='Search for similar compounds')
    search_parser.add_argument('smiles', help='Query SMILES string')
    search_parser.add_argument('-k', '--top-k', type=int, default=5, help='Number of results to return')
    search_parser.add_argument('--mol2vec-model', help='Path to Mol2Vec model for blended search')
    
    args = parser.parse_args()
    
    if args.command == 'build':
        from src.index_builder import build_ecfp6_index, build_mol2vec_index
        import os
        
        input_path = args.input_file
        if not os.path.exists(input_path):
            print(f"Error: Input file not found: {input_path}")
            sys.exit(1)
        
        if not input_path.endswith('.smi'):
            print("Error: Only .smi files are supported")
            print("Expected format: SMILES compound_name (space-separated)")
            sys.exit(1)
        
        print(f"Building index from: {input_path}")
        
        try:
            # Build ECFP6 index
            build_ecfp6_index(input_path, args.output_dir)
            
            # Build Mol2Vec index if model provided
            if args.mol2vec_model:
                if not os.path.exists(args.mol2vec_model):
                    print(f"Warning: Mol2Vec model not found: {args.mol2vec_model}")
                    print("Skipping Mol2Vec index building...")
                else:
                    build_mol2vec_index(input_path, args.mol2vec_model, args.output_dir)
                
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
