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
        from scripts.build_index import main as build_main
        sys.argv = ['build_index.py', args.input_file]
        build_main()
    elif args.command == 'search':
        from scripts.nearest import main as search_main
        sys.argv = ['nearest.py', args.smiles, str(args.top_k)]
        if args.mol2vec_model:
            sys.argv.append(args.mol2vec_model)
        search_main()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
