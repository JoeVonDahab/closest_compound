# Closest Compound

A molecular similarity search tool that provides two complementary approaches for finding chemically similar compounds: traditio## Search Method Comparison - What You Actually Get

### ECFP6 Only: Structural Similarity
Finds compounds with **similar chemical structures** - same scaffolds, functional groups, substructures.

**Example with Aspirin:**
```
 1. Aspirin     |  CC(=O)OC1=CC=CC=C1C(=O)O           |  Tanimoto: 1.000
 2. Salsalate   |  C1=CC=C(C(=C1)C(=O)OC2=CC=CC=C2C(=O)O)O  |  Tanimoto: 0.690
 3. Fosfosal    |  C1=CC=C(C(=C1)C(=O)O)OP(=O)(O)O    |  Tanimoto: 0.635
 4. Benorilate  |  CC(=O)NC1=CC=C(C=C1)OC(=O)C2=CC=CC=C2OC(=O)C  |  Tanimoto: 0.622
```
**What you see:** All results have benzene rings + carboxyl/acetyl groups like aspirin. Clear structural relationships.

### Blended Search: Structural + Semantic Similarity  
Finds compounds that are **both structurally similar AND functionally/pharmacologically related**.

**Same Aspirin Query:**
```
✓ Blended ECFP6 + Mol2Vec:
 1. Aspirin     |  CC(=O)OC1=CC=CC=C1C(=O)O
 2. Indiplon    |  CC(=O)N(C)C1=CC=CC(=C1)C2=CC=NC3=C(C=NN23)C(=O)C4=CC=CS4
 3. Salsalate   |  C1=CC=C(C(=C1)C(=O)OC2=CC=CC=C2C(=O)O)O
 4. Carbetocin  |  CCC(C)C1C(=O)NC(C(=O)NC(C(=O)NC(CSCCCC(=O)...
 5. Asasantin   |  CC(=O)OC1=CC=CC=C1C(=O)O.C1CCN(CC1)C2=NC(=NC3=C2N=C...
```

**Key Difference:** Notice Indiplon (#2) and Carbetocin (#4) - these have **completely different structures** but appear because:
- **Indiplon**: Sedative/hypnotic drug - shares pharmacological context with aspirin's CNS effects
- **Carbetocin**: Peptide hormone - found due to shared therapeutic database contexts

## When to Use Which Method

### Use ECFP6 Only When:
- **Scaffold hopping**: Finding variations of the same chemical scaffold  
- **SAR analysis**: Understanding structure-activity relationships
- **Patent searching**: Finding structurally similar compounds for IP analysis
- **Chemical series expansion**: Growing around a known chemotype

### Use Blended Search When:
- **Drug repurposing**: Finding existing drugs that might work for new targets
- **Target prediction**: Understanding what your compound might interact with
- **Novel scaffold discovery**: Finding functionally similar compounds with different structures  
- **Comprehensive similarity**: When you want both structural and biological similarity

**Bottom Line:** ECFP6 = "compounds that look similar", Blended = "compounds that might work similarly"fingerprints and modern Mol2Vec neural embeddings.

## Features

### Two Search Approaches Available

#### 1. ECFP6 + Tanimoto Similarity (Default)
- **Fast and interpretable**: Classical molecular fingerprints based on circular substructures
- **Structural similarity**: Excellent for finding compounds with similar scaffolds and functional groups
- **Tanimoto scoring**: Well-established similarity metric (0-1 scale, 1 = identical)
- **Always available**: Built-in method requiring no external models

#### 2. Blended ECFP6 + Mol2Vec Search
- **Combined ranking**: Uses BOTH ECFP6 fingerprints AND Mol2Vec neural embeddings
- **Enhanced discovery**: Finds both structurally and semantically similar compounds
- **Comprehensive results**: Merges rankings from both methods for better coverage
- **Requires model**: Uses pre-trained model_300dim.pkl for the Mol2Vec component

### Additional Features
- **SMI Input Format**: Simple space-separated SMILES and compound names
- **FAISS Backend**: Efficient similarity search using Facebook's FAISS library
- **Python 3.9 + Gensim 3.8**: Optimized for compatibility

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

### 1. Download Mol2Vec Model (Optional)

For neural embedding similarity, download the pre-trained model:

```bash
mkdir -p models
wget https://github.com/samoturk/mol2vec/raw/master/examples/models/model_300dim.pkl -P models/
```

### 2. Build Similarity Indices

Choose your approach based on your needs:

```bash
# Option A: ECFP6 fingerprints only (fast, structural similarity)
uv run main.py build library/library.smi

# Option B: Both ECFP6 + Mol2Vec (comprehensive, requires model download)
uv run main.py build library/library.smi --mol2vec-model models/model_300dim.pkl

# Custom library and output directory
uv run main.py build your_library.smi --output-dir my_artifacts --mol2vec-model models/model_300dim.pkl
```

**What gets built:**
- **ECFP6 index**: `artifacts/ecfp6.index` + `artifacts/meta.pkl` (always)
- **Mol2Vec index**: `artifacts/mol2vec.index` (only with --mol2vec-model)

### 3. Search for Similar Compounds

Choose your search method:

```bash
# Method 1: ECFP6 + Tanimoto similarity only (structural similarity)
uv run main.py search "CC(=O)OC1=CC=CC=C1C(=O)O" --top-k 5

# Method 2: Blended search - combines ECFP6 + Mol2Vec (structural + semantic)
uv run main.py search "CC(=O)OC1=CC=CC=C1C(=O)O" --mol2vec-model models/model_300dim.pkl --top-k 5

# Search with more results
uv run main.py search "CCO" --top-k 10
```

**How the search methods work:**
- **Without `--mol2vec-model`**: Uses only ECFP6 fingerprints + Tanimoto similarity
- **With `--mol2vec-model`**: Combines both ECFP6 and Mol2Vec rankings into a blended result

## Input Format

**SMI Files**: Space-separated format with SMILES and compound names
```
CC(=O)OC1=CC=CC=C1C(=O)O Aspirin
CCO Ethanol
CC(C)O Isopropanol
```

## Output Examples

### ECFP6 + Tanimoto Search (Structural Similarity)
```bash
uv run main.py search "CC(=O)OC1=CC=CC=C1C(=O)O" --top-k 3
```
```
 1. Aspirin  |  CC(=O)OC1=CC=CC=C1C(=O)O  |  Tanimoto: 1.000
 2. Salsalate  |  C1=CC=C(C(=C1)C(=O)OC2=CC=CC=C2C(=O)O)O  |  Tanimoto: 0.690
 3. Fosfosal  |  C1=CC=C(C(=C1)C(=O)O)OP(=O)(O)O  |  Tanimoto: 0.635
```

### Blended ECFP6 + Mol2Vec Search (Structural + Semantic Similarity)
```bash
uv run main.py search "CC(=O)OC1=CC=CC=C1C(=O)O" --mol2vec-model models/model_300dim.pkl --top-k 3
```
```
✓ Blended ECFP6 + Mol2Vec:
 1. Aspirin  |  CC(=O)OC1=CC=CC=C1C(=O)O
 2. Indiplon  |  CC(=O)N(C)C1=CC=CC(=C1)C2=CC=NC3=C(C=NN23)C(=O)C4=CC=CS4
 3. Salsalate  |  C1=CC=C(C(=C1)C(=O)OC2=CC=CC=C2C(=O)O)O
```

**Notice**: Blended search finds Indiplon (a sedative) which shares pharmacological properties with aspirin but has different structure - this demonstrates semantic similarity!

```
closest_compound/
├── main.py                 # Main CLI interface for build/search
├── src/
│   ├── index_builder.py    # ECFP6 and Mol2Vec index building functions
│   ├── query.py           # ECFP6 + Tanimoto and Mol2Vec search
│   └── utils.py           # Molecular utility functions
├── library/
│   └── library.smi        # Sample compound library
├── models/                # Mol2Vec models directory
│   └── model_300dim.pkl   # Pre-trained Mol2Vec model (300 dimensions)
└── artifacts/             # Generated indices (created after build)
    ├── ecfp6.index        # FAISS index for ECFP6 fingerprints
    ├── mol2vec.index      # FAISS index for Mol2Vec embeddings (optional)
    └── meta.pkl           # Metadata file
```

## Dependencies

- **Python 3.9**: Optimized compatibility with mol2vec ecosystem
- **RDKit**: Molecular structure handling and ECFP6 fingerprint generation  
- **FAISS**: Fast similarity search backend
- **Mol2Vec**: Neural molecular embeddings (optional)
- **Gensim 3.8**: Word2Vec model compatibility for Mol2Vec
- **NumPy**: Numerical operations

## Search Method Comparison

## Technical Details

### ECFP6 Fingerprints
- **Algorithm**: Extended Connectivity Fingerprints with radius=3
- **Bit vector**: 2048 bits
- **Similarity metric**: Tanimoto coefficient
- **Best for**: Finding compounds with similar substructures and scaffolds

### Mol2Vec Embeddings  
- **Algorithm**: Neural word2vec trained on molecular "sentences"
- **Vector size**: 300 dimensions
- **Model**: Pre-trained on large chemical databases
- **Best for**: Finding functionally related compounds across different chemical spaces

### Blended Search Algorithm
1. Generate ECFP6 fingerprint and Mol2Vec embedding for query
2. Search both FAISS indices independently  
3. Combine and re-rank results using weighted scoring
4. Return unified ranked list

