# Closest Compound

A molecular similarity search tool that uses ECFP6 fingerprints and optionally Mol2Vec embeddings to find chemically similar compounds in a molecular library.

## Features

- **ECFP6 Fingerprints + Tanimoto Similarity**: Fast and accurate molecular similarity using Extended Connectivity Fingerprints with Tanimoto coefficient scoring
- **Mol2Vec Neural Embeddings**: Optional molecular embeddings for enhanced similarity search
- **Blended Search**: Combines ECFP6 and Mol2Vec methods for improved results
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

Build a similarity index from your molecular library:

```bash
# Using the provided library
uv run main.py build library/library.smi

# Using your own CSV file
uv run main.py build your_library.csv --output-dir artifacts
```

**Input Formats:**
- **CSV**: Must contain a `smiles` column and preferably a `drug_name` column
- **SMI**: Space-separated format with SMILES and compound names (e.g., "CCO Ethanol")

### 2. Search for Similar Compounds

```bash
# ECFP6 + Tanimoto similarity search (top 5 results)
uv run main.py search "CCO"

# Search for top 10 compounds
uv run main.py search "CCO" --top-k 10

# Blended ECFP6 + Mol2Vec search (if you have a trained model)
uv run main.py search "CCO" --top-k 10 --mol2vec-model /path/to/mol2vec_model.pkl
```

## Project Structure

```
closest_compound/
├── main.py                 # Main CLI interface for build/search
├── src/
│   ├── index_builder.py    # ECFP6 index building functions
│   ├── query.py           # ECFP6 + Tanimoto and Mol2Vec search
│   └── utils.py           # Molecular utility functions
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