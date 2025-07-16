# Closest Compound

A molecular similarity search tool that provides two complementary approaches for finding chemically similar compounds: traditional ECFP6 fingerprints and modern Mol2Vec neural embeddings.

## Installation

```bash
# Clone the repository (includes pre-trained Mol2Vec model)
git clone https://github.com/JoeVonDahab/closest_compound.git
cd closest_compound

# Install dependencies using UV
uv sync
```

## Quick Start

### 1. Build Similarity Indices

```bash
# Option A: ECFP6 fingerprints only (structural similarity)
uv run main.py build library/library.smi

# Option B: Both ECFP6 + Mol2Vec (structural + semantic similarity)
uv run main.py build library/library.smi --mol2vec-model models/model_300dim.pkl
```

### 2. Search for Similar Compounds

```bash
# ECFP6 only - finds structurally similar compounds
uv run main.py search "CC(=O)OC1=CC=CC=C1C(=O)O"

# Blended search - finds both structurally and functionally similar compounds  
uv run main.py search "CC(=O)OC1=CC=CC=C1C(=O)O" --mol2vec-model models/model_300dim.pkl
```
## Difference between two methods: 

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


