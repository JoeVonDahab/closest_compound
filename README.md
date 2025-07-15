# Closest Compound

A molecular similarity search tool that uses ECFP6 fingerprints and optionally Mol2Vec embeddings to find chemically similar compounds in a molecular library.

## Features

- **ECFP6 Fingerprints**: Fast similarity search using Extended Connectivity Fingerprints
- **Mol2Vec Integration**: Optional neural embedding-based similarity (when model provided)
- **Blended Search**: Combines both methods for improved results
- **Multiple Input Formats**: Supports both CSV and SMI file formats
- **FAISS Backend**: Efficient similarity search using Facebook's FAISS library

## Installation

```bash
# Clone the repository
git clone <your-repo-url>
cd closest_compound

# Install dependencies using UV
uv sync

# Or install in development mode
uv pip install -e .
```

## Quick Start

### 1. Build the Index

First, build a similarity index from your molecular library:

```bash
# Using the provided library
uv run scripts/build_index.py library/library.smi

# Or using your own CSV file
uv run scripts/build_index.py your_library.csv
```

**Input Formats:**
- **CSV**: Must contain a `smiles` column and preferably a `drug_name` column
- **SMI**: Tab-separated format with SMILES and compound names

### 2. Search for Similar Compounds

```bash
# Basic search for top 5 similar compounds
uv run scripts/nearest.py "CCO"

# Search for top 10 compounds
uv run scripts/nearest.py "CCO" 10

# Blended search using Mol2Vec (if you have a trained model)
uv run scripts/nearest.py "CCO" 10 /path/to/mol2vec_model.pkl
```

### 3. Using the Main Interface

```bash
# Build index
uv run main.py build library/library.smi

# Search
uv run main.py search "CCO" --top-k 10
uv run main.py search "CCO" --top-k 10 --mol2vec-model /path/to/model.pkl
```

## Project Structure

```
closest_compound/
├── main.py                 # Main CLI interface
├── scripts/
│   ├── build_index.py      # Index building script
│   └── nearest.py          # Search script
├── src/
│   ├── index_builder.py    # Index building functions
│   ├── query.py           # Search functions
│   └── utils.py           # Utility functions
├── library/
│   └── library.smi        # Sample compound library
└── artifacts/             # Generated indices (created after build)
    ├── ecfp6.index        # FAISS index file
    └── meta.pkl           # Metadata file
```

## Dependencies

- **RDKit**: Molecular structure handling and fingerprint generation
- **FAISS**: Fast similarity search
- **Mol2Vec**: Optional neural embeddings for molecules
- **NumPy/Pandas**: Data manipulation
- **Gensim**: Word2Vec model handling (for Mol2Vec)

## Output Format

Search results include:
- **Rank**: Position in similarity ranking
- **Compound Name**: Name from the library
- **SMILES**: Chemical structure representation
- **Tanimoto Score**: Similarity score (0-1, higher = more similar)

## Notes

- ECFP6 fingerprints use radius=3 and 2048 bits
- Similarity scores are Tanimoto coefficients
- Index files are stored in the `artifacts/` directory
- The tool validates SMILES strings before processing

## Example Output

```
Searching for compounds similar to: CCO
Returning top 5 results
--------------------------------------------------
 1. Ethanol  |  CCO  |  Tanimoto: 1.000
 2. Methanol  |  CO  |  Tanimoto: 0.400
 3. Propanol  |  CCCO  |  Tanimoto: 0.600
 ...
```