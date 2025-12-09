# Comprehensive Lecture Notebooks - Implementation Plan

## ✅ Status: Lecture 1 Complete

**Created:** `notebooks/lecture01_intro_neurogenomics.ipynb`

### Features of Lecture 1:
- ✅ Professional MIT/Harvard-level formatting
- ✅ Real neuroscience data context (Allen Brain Atlas references)
- ✅ Comprehensive explanations of concepts
- ✅ Well-commented, executable code
- ✅ Detailed output interpretations
- ✅ Citations from high-impact journals (Nature, Science)
- ✅ Clinical applications and disease context
- ✅ Homework assignment included
- ✅ Additional resources and references
- ✅ ~400 lines of content with rich markdown

---

## 📚 Complete Lecture Plan (12 Lectures)

### Module 1: Introduction to Single-Cell Neurogenomics

#### **Lecture 1: Introduction to Neurogenomics** ✅ COMPLETE
**File:** `lecture01_intro_neurogenomics.ipynb`
- Genetic basis of brain function
- Brain cell type diversity
- Gene expression patterns
- Disease-associated genes
- Real data exploration

#### **Lecture 2: Introduction to Single-Cell Technology** 🔄 IN PROGRESS
**File:** `lecture02_single_cell_technology.ipynb`
**Content Plan:**
- Principles of single-cell genomics
- 10X Genomics technology
- Droplet microfluidics
- UMIs and barcoding
- Data structure (AnnData)
- **Dataset:** PBMC 3k/10k (10X Genomics, publicly available)
- **Code:** Loading, exploring, visualizing scRNA-seq data
- **Citation:** Zheng et al. (2017) *Nature Comm* 8:14049

#### **Lecture 3: Fundamentals of Python** 🔄 TO CREATE
**File:** `lecture03_python_fundamentals.ipynb`
**Content Plan:**
- Python basics for bioinformatics
- NumPy arrays and operations
- Pandas DataFrames
- Data visualization with matplotlib/seaborn
- Working with biological data
- **Dataset:** Gene expression matrix (simulated/small real dataset)
- **Code:** Practical exercises with increasing complexity
- Functions, loops, data manipulation

#### **Lecture 4: Quantification of Single-Cell RNA-seq Data** 🔄 TO CREATE
**File:** `lecture04_quantification.ipynb**
**Content Plan:**
- FASTQ file structure
- Read alignment concepts
- UMI counting
- Gene-cell matrix generation
- **Dataset:** Small FASTQ example + pre-quantified data
- **Code:** Interpreting CellRanger/kb-python outputs
- **Tools:** Demonstrate with real output files

---

### Module 2: Core Single-Cell RNA-seq Analysis

#### **Lecture 5: Quality Control, Normalization, and Preprocessing** 🔄 TO CREATE
**File:** `lecture05_qc_preprocessing.ipynb`
**Content Plan:**
- QC metrics (UMI counts, genes detected, MT%)
- Cell filtering strategies
- Doublet detection
- Normalization methods (CPM, log-normalization)
- Feature selection (HVGs)
- **Dataset:** Mouse brain cortex (Zeisel et al. 2015, *Science*)
- **Code:** Complete QC pipeline with scanpy
- **Citation:** Zeisel et al. (2015) *Science* 347:1138-1142

#### **Lecture 6: Downstream Analysis and Batch Correction** 🔄 TO CREATE
**File:** `lecture06_analysis_integration.ipynb`
**Content Plan:**
- PCA dimensionality reduction
- UMAP/t-SNE visualization
- Leiden clustering
- Marker gene identification
- Batch effect detection and correction
- **Dataset:** Multi-batch brain data or PBMC multi-sample
- **Code:** Complete analysis workflow
- **Methods:** Harmony, scVI, Seurat integration

---

### Module 3: Advanced Single-Cell Analysis

#### **Lecture 7: Trajectory Inference and Fate Probability** 🔄 TO CREATE
**File:** `lecture07_trajectory_inference.ipynb`
**Content Plan:**
- Pseudotime concepts
- RNA velocity theory
- CellRank for fate mapping
- Lineage tracing
- **Dataset:** Neuronal differentiation (LaManno et al. or pancreas development)
- **Code:** scVelo + CellRank pipeline
- **Citation:** LaManno et al. (2018) *Nature* 560:494-498

#### **Lecture 8: Ligand-Receptor Interactions and Differential Abundance** 🔄 TO CREATE
**File:** `lecture08_cell_communication.ipynb`
**Content Plan:**
- Cell-cell communication theory
- Ligand-receptor databases (CellPhoneDB, LIANA)
- Interaction significance testing
- Network visualization
- Differential abundance methods
- **Dataset:** Brain immune interactions or PBMC
- **Code:** LIANA/CellPhoneDB analysis
- **Citation:** Efremova et al. (2020) *Nat Protoc* 15:1484-1506

#### **Lecture 9: Spatial Analysis and Spatial Structure** 🔄 TO CREATE
**File:** `lecture09_spatial_transcriptomics.ipynb`
**Content Plan:**
- Spatial transcriptomics technologies (Visium, MERFISH, Xenium)
- Spatial data structures
- Spatial clustering
- Niche identification
- Spatially variable genes
- **Dataset:** Mouse brain Visium (Allen Brain Atlas or 10X)
- **Code:** Squidpy analysis pipeline
- **Citation:** Moses & Pachter (2022) *Nat Rev Genet* 23:601-618

---

### Module 4: Deep Learning for Single-Cell Genomics

#### **Lecture 10: Deep Learning Part 1** 🔄 TO CREATE
**File:** `lecture10_deep_learning_part1.ipynb`
**Content Plan:**
- Neural network basics
- Variational autoencoders (VAE) theory
- scVI architecture
- Latent space interpretation
- Feature extraction
- **Dataset:** PBMC 10k or brain dataset
- **Code:** Training scVI model, extracting latent representation
- **Citation:** Lopez et al. (2018) *Nat Methods* 15:1053-1058

#### **Lecture 11: Deep Learning Part 2** 🔄 TO CREATE
**File:** `lecture11_deep_learning_part2.ipynb`
**Content Plan:**
- scANVI for semi-supervised learning
- Transfer learning concepts
- Data integration with deep learning
- Model interpretation
- Comparison with classical methods
- **Dataset:** Multi-batch dataset
- **Code:** scANVI training and prediction
- **Citation:** Gayoso et al. (2021) *Nat Biotechnol* 39:30-34

#### **Lecture 12: Single-Cell Neurogenomics Project** 🔄 TO CREATE
**File:** `lecture12_neurogenomics_project.ipynb`
**Content Plan:**
- Complete workflow from raw data to interpretation
- Brain cell type atlas construction
- Disease gene analysis
- Publication-quality figures
- Biological interpretation
- **Dataset:** Human/mouse brain scRNA-seq (published study)
- **Code:** End-to-end analysis pipeline
- **Focus:** Integration of all learned methods

---

## 🎯 Notebook Structure Template

Each notebook will follow this professional structure:

### 1. **Header Section**
```markdown
# Lecture X: [Title]
**Course:** BRAIN
**Date:** [Date]
**Duration:** 90 minutes
**Learning Objectives:** [List]
```

### 2. **Table of Contents**
- Clickable links to all sections
- Clear organization

### 3. **Introduction** (~500 words)
- Biological/computational background
- Why this topic matters
- Key concepts preview
- Historical context

### 4. **Setup Cell**
```python
# Library imports with version printing
# Configuration settings
# Reproducibility setup (seeds)
```

### 5. **Conceptual Sections** (3-5 major sections)
- Markdown explanation (300-500 words per section)
- Visual diagrams where appropriate
- Key terminology definitions
- Citations to primary literature

### 6. **Code Sections** (Integrated with concepts)
```python
# Detailed comments explaining every step
# Print statements showing intermediate results
# Assertions for data validation
# Error handling where appropriate
```

### 7. **Output Interpretation**
- Markdown cells after each code block
- Explain what the output means biologically
- Common patterns to expect
- What to do if results differ

### 8. **Visualization Sections**
```python
# Professional matplotlib/seaborn plots
# Clear labels, legends, titles
# Color schemes (colorblind-friendly)
# Multiple subplots when appropriate
```

### 9. **Interpretation Blocks**
- Biological meaning of results
- Comparison to literature
- Common pitfalls
- Best practices

### 10. **Disease/Clinical Context**
- How methods apply to disease studies
- Example applications
- Recent publications

### 11. **Summary Section**
- Key takeaways (5-10 bullet points)
- Concepts learned
- Skills acquired
- Preview of next lecture

### 12. **Additional Resources**
- Primary literature (3-5 papers)
- Online tutorials
- Databases
- Software documentation

### 13. **Homework Assignment**
- 3-5 practical tasks
- Grading rubric
- Submission guidelines

---

## 📊 Data Sources (All from High-Impact Publications)

### Mouse Brain Data
1. **Zeisel et al. (2015)** - Mouse cortex/hippocampus
   - *Science* 347:1138-1142
   - 3,005 cells, 19,972 genes
   - Available: https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt

2. **Zeisel et al. (2018)** - Mouse nervous system
   - *Cell* 174:999-1014
   - 560,000 cells
   - Available: http://mousebrain.org/

3. **Yao et al. (2023)** - Allen Mouse Brain Cell Atlas
   - *Nature* 624:317-332
   - 4.3 million cells
   - Available: https://portal.brain-map.org/

### Human Brain Data
4. **Mathys et al. (2019)** - Alzheimer's disease
   - *Nature* 570:332-337
   - 80,660 cells from prefrontal cortex
   - Available: Synapse: syn18485175

5. **Hodge et al. (2019)** - Human cortex
   - *Nature* 573:61-68
   - Multiple cortical areas
   - Available: Allen Brain Map portal

### Development Data
6. **Nowakowski et al. (2017)** - Human cortical development
   - *Science* 358:1318-1323
   - Fetal brain development
   - Available: GEO: GSE99951

### Disease Data
7. **Kamath et al. (2022)** - Parkinson's disease
   - *Nat Neurosci* 25:588-595
   - Human substantia nigra
   - Available: Synapse: syn25921710

### Public Datasets (10X Genomics)
8. **PBMC datasets** - Multiple versions
   - 3k, 10k, 68k cells
   - Freely available
   - Excellent for methods development

---

## 🔧 Technical Requirements

### Software Environment
```bash
# Create conda environment
conda create -n brain-course python=3.9
conda activate brain-course

# Install core packages
pip install scanpy==1.9.3
pip install scvi-tools==1.0.0
pip install squidpy==1.3.0
pip install cellrank==2.0.0
pip install scvelo==0.2.5
pip install liana==0.1.0

# Visualization
pip install matplotlib==3.7.0
pip install seaborn==0.12.2

# Utilities
pip install numpy==1.24.0
pip install pandas==2.0.0
```

### Hardware Requirements
- **RAM:** 16GB minimum (32GB recommended)
- **Storage:** 50GB free space
- **CPU:** Multi-core processor
- **GPU:** Optional (for deep learning lectures)

### Data Storage
```
data/
├── raw/
│   ├── pbmc_3k/
│   ├── mouse_brain/
│   └── human_brain/
├── processed/
│   ├── qc_filtered/
│   └── normalized/
└── results/
    ├── figures/
    └── tables/
```

---

## 📝 Code Quality Standards

### 1. **Documentation**
```python
def process_brain_data(adata, min_genes=200, max_mt=5):
    """
    Perform quality control on brain single-cell data.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations.
    min_genes : int, optional (default: 200)
        Minimum number of genes expressed per cell.
    max_mt : float, optional (default: 5)
        Maximum percentage of mitochondrial gene expression.

    Returns
    -------
    AnnData
        Filtered dataset with QC metrics added to .obs.

    Examples
    --------
    >>> adata = sc.read_h5ad('brain_data.h5ad')
    >>> adata_filtered = process_brain_data(adata, min_genes=300)
    >>> print(f"Kept {adata_filtered.n_obs} cells")

    Notes
    -----
    This function follows best practices from:
    Luecken & Theis (2019). Mol Syst Biol 15:e8746
    """
    # Implementation
    pass
```

### 2. **Comments**
```python
# Calculate mitochondrial gene percentage
# Rationale: High MT% indicates dying/damaged cells
# Threshold: Typically 5-20% depending on tissue
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
```

### 3. **Print Statements**
```python
print("="*70)
print("Quality Control Summary")
print("="*70)
print(f"Cells before filtering: {n_cells_before:,}")
print(f"Cells after filtering:  {n_cells_after:,}")
print(f"Cells removed:          {n_cells_before - n_cells_after:,} ({pct_removed:.1f}%)")
print("="*70)
```

### 4. **Error Handling**
```python
try:
    adata = sc.read_h5ad(file_path)
except FileNotFoundError:
    print(f"Error: Could not find file '{file_path}'")
    print("Please check:")
    print("  1. File path is correct")
    print("  2. File has been downloaded")
    print("  3. You have read permissions")
    sys.exit(1)
```

### 5. **Reproducibility**
```python
# Set random seeds for reproducibility
import numpy as np
import random
np.random.seed(42)
random.seed(42)

# Scanpy settings
sc.settings.verbosity = 2  # Errors and warnings only
sc.settings.n_jobs = -1    # Use all available cores
```

---

## 🎨 Visualization Standards

### Professional Plotting
```python
# Set publication-quality parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.5

# Color schemes (colorblind-friendly)
colors_celltypes = sns.color_palette('Set2', 8)
colors_expression = 'viridis'

# Example multi-panel figure
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Add panel labels
for i, ax in enumerate(axes.flat):
    ax.text(-0.1, 1.1, string.ascii_uppercase[i],
            transform=ax.transAxes, size=14, weight='bold')

# Save high-resolution
plt.savefig('figure.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figure.png', dpi=300, bbox_inches='tight')
```

---

## ✅ Completion Checklist

### For Each Notebook:

- [ ] Professional header with metadata
- [ ] Clear learning objectives
- [ ] Table of contents with links
- [ ] Comprehensive introduction (~500 words)
- [ ] All libraries imported with versions
- [ ] Real neuroscience data loaded
- [ ] Data source cited (paper + accession)
- [ ] Concepts explained before code
- [ ] Code fully commented
- [ ] Outputs interpreted biologically
- [ ] High-quality visualizations
- [ ] Clinical/disease context included
- [ ] Summary with key takeaways
- [ ] Additional resources listed
- [ ] Homework assignment included
- [ ] Notebook runs without errors
- [ ] All plots display correctly
- [ ] Professional formatting throughout

---

## 🚀 Next Steps

### Immediate Actions:
1. ✅ Lecture 1 complete - review and refine if needed
2. 🔄 Create Lecture 2 (single-cell technology)
3. 🔄 Create Lecture 5 (QC/preprocessing) - critical foundation
4. 🔄 Create Lecture 10 (deep learning) - modern methods

### Priority Order for Remaining Lectures:
1. **Lecture 5** (QC) - Most critical, used in all subsequent analyses
2. **Lecture 6** (Clustering) - Core analysis skill
3. **Lecture 2** (Technology) - Foundation understanding
4. **Lecture 7** (Trajectory) - Advanced but increasingly important
5. **Lecture 10-11** (Deep learning) - Cutting-edge methods
6. **Lecture 3** (Python) - Can use existing resources
7. **Lecture 4** (Quantification) - Technical but important
8. **Lecture 8-9** (Communication, Spatial) - Specialized topics
9. **Lecture 12** (Project) - Integrative capstone

### Timeline Estimate:
- Each comprehensive notebook: ~4-6 hours to create
- Testing and refinement: ~2 hours per notebook
- Total for all 12: ~72-96 hours of development

---

## 📚 References for All Lectures

### Key Methodology Papers
1. Luecken & Theis (2019). Current best practices in scRNA-seq analysis. *Mol Syst Biol* 15:e8746
2. Hie et al. (2019). Efficient integration of heterogeneous scRNA-seq datasets. *Nat Biotechnol* 37:685-691
3. Korsunsky et al. (2019). Fast, sensitive and accurate integration. *Nat Methods* 16:1289-1296

### Brain Atlas Papers
4. Yao et al. (2023). High-resolution mouse brain atlas. *Nature* 624:317-332
5. Siletti et al. (2023). Human brain cell atlas. *Science* 382:eadd7046

### Disease Papers
6. Mathys et al. (2019). Alzheimer's disease scRNA-seq. *Nature* 570:332-337
7. Kamath et al. (2022). Parkinson's disease. *Nat Neurosci* 25:588-595

### Technology Papers
8. Zheng et al. (2017). Massively parallel scRNA-seq. *Nat Comm* 8:14049
9. Macosko et al. (2015). Drop-Seq. *Cell* 161:1202-1214

### Methods Papers
10. Bergen et al. (2020). RNA velocity. *Nat Biotechnol* 38:1408-1414
11. Lange et al. (2022). CellRank. *Nat Methods* 19:159-170
12. Palla et al. (2022). Squidpy. *Nat Methods* 19:171-178

---

*Document created: November 2, 2025*
*Status: Lecture 1 complete, 11 lectures planned*
*Next update: Upon completion of next lecture*
