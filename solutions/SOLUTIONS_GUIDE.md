# Solutions Guide for Single-Cell Neurogenomics Assignments

## Overview

This directory contains complete solution notebooks for all 12 course assignments. Solutions include:
- ✓ Fully working code for all tasks
- ✓ Expected outputs and visualizations
- ✓ Detailed explanations and interpretations
- ✓ Answers to reflection questions
- ✓ Best practices and tips

---

## How to Use Solutions

### For Instructors
1. **Reference during grading**: Check expected outputs and approaches
2. **Adapt for demonstrations**: Use code snippets in lectures
3. **Customize for your course**: Modify parameters or datasets
4. **Provide as hints**: Share specific sections with struggling students

### For Students (After Submission)
1. **Compare your approach**: See alternative solutions
2. **Learn best practices**: Study code organization and style
3. **Understand concepts**: Read explanations and interpretations
4. **Prepare for exams**: Review key techniques and methods

### Important Note for Students
⚠️ **Academic Integrity**: Solutions should only be consulted AFTER submitting your own work. Using solutions before completing assignments constitutes plagiarism.

---

## Solution Files

### Complete Solutions Available

| Lecture | Solution File | Key Topics |
|---------|--------------|------------|
| 1 | lecture01_neurogenomics_solution.ipynb | Data structures, visualization, correlations |
| 2 | lecture02_single_cell_technology_solution.ipynb | AnnData, sparsity, QC metrics |
| 3 | lecture03_python_fundamentals_solution.ipynb | Python basics, functions, data structures |
| 4-12 | In development | Core and advanced scRNA-seq analysis |

---

## Solution Approaches by Lecture

### Lecture 1: Introduction to Neurogenomics

**Key Concepts**:
- Creating and manipulating pandas DataFrames
- Gene expression visualization with seaborn
- Calculating correlations and fold-changes
- Biological interpretation of expression patterns

**Common Mistakes**:
- Not using proper indexing for genes/regions
- Forgetting to add labels to plots
- Incorrect fold-change calculations (max/min vs max-min)

**Best Practices**:
```python
# Good: Clear variable names and comments
gene_expression = pd.DataFrame(data, index=genes, columns=regions)

# Good: Informative plot labels
plt.xlabel('Brain Region', fontsize=12)
plt.title('Gene Expression Patterns', fontweight='bold')
```

---

### Lecture 2: Introduction to Single-Cell Technology

**Key Concepts**:
- Understanding AnnData structure (.X, .obs, .var, .uns)
- Working with sparse matrices
- Calculating and interpreting QC metrics
- Comparing single-cell vs bulk approaches

**Common Mistakes**:
- Not converting sparse matrices before calculations
- Forgetting to handle mitochondrial genes separately
- Misinterpreting sparsity (biological vs technical)

**Best Practices**:
```python
# Good: Proper sparse matrix handling
if hasattr(adata.X, 'toarray'):
    X_dense = adata.X.toarray()
else:
    X_dense = adata.X

# Good: Calculating QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
```

---

### Lecture 3: Fundamentals of Python

**Key Concepts**:
- Data types and type conversion
- List comprehensions for efficient operations
- Dictionary operations for biological data
- Writing reusable functions with docstrings

**Common Mistakes**:
- Using mutable default arguments in functions
- Not handling edge cases (empty lists, None values)
- Inefficient loops instead of comprehensions

**Best Practices**:
```python
# Good: List comprehension
normalized = [x / max_val for x in expression]

# Good: Function with docstring and error handling
def calculate_metrics(data):
    """Calculate statistics for gene expression data.

    Parameters:
        data: list of numeric values
    Returns:
        dict with mean, median, min, max
    """
    if not data:
        return None
    # ... rest of function
```

---

### Lecture 4: Quantification of scRNA-seq Data

**Key Concepts**:
- Understanding count matrix structure
- Barcode rank plots and knee detection
- QC metric interpretation
- Sequencing saturation analysis

**Common Mistakes**:
- Not using log scale for barcode rank plots
- Incorrect filtering thresholds
- Missing relationship between metrics

**Best Practices**:
```python
# Good: Barcode rank plot
total_counts = np.array(adata.X.sum(axis=1)).flatten()
plt.plot(range(len(total_counts)),
         np.sort(total_counts)[::-1])
plt.yscale('log')
plt.xlabel('Cells (ranked by total counts)')
```

---

### Lecture 5: QC, Normalization, and Preprocessing

**Key Concepts**:
- Cell and gene filtering strategies
- Library size normalization
- Log transformation rationale
- Highly variable gene selection
- PCA for dimensionality reduction

**Common Mistakes**:
- Filtering too aggressively or too loosely
- Not storing raw counts before normalization
- Using all genes instead of HVGs for PCA

**Best Practices**:
```python
# Good: Store raw counts
adata.raw = adata

# Good: Standard normalization pipeline
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
```

---

### Lecture 6: Downstream Analysis and Batch Correction

**Key Concepts**:
- Neighborhood graph construction
- Leiden clustering algorithm
- Marker gene identification
- Cell type annotation strategies
- Batch effect correction with scVI

**Common Mistakes**:
- Using inappropriate clustering resolution
- Not validating cluster biological meaning
- Over-relying on automatic annotation

**Best Practices**:
```python
# Good: Systematic clustering evaluation
for resolution in [0.4, 0.8, 1.2]:
    sc.tl.leiden(adata, resolution=resolution, key_added=f'leiden_{resolution}')

# Good: Marker gene identification
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20)
```

---

### Lecture 7: Trajectory Inference

**Key Concepts**:
- RNA velocity computation
- Dynamical vs stochastic models
- CellRank for fate mapping
- Pseudotime ordering

**Common Mistakes**:
- Not checking if dataset suitable for trajectory analysis
- Misinterpreting velocity direction
- Over-interpreting rare trajectories

**Best Practices**:
```python
# Good: Velocity workflow
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')
```

---

### Lecture 8: Cell Communication

**Key Concepts**:
- Ligand-receptor interaction databases
- Statistical significance testing
- Interaction specificity vs magnitude
- Differential abundance analysis

**Common Mistakes**:
- Not filtering for significant interactions
- Ignoring biological context
- Wrong statistical test for abundance

**Best Practices**:
```python
# Good: LIANA analysis
li.mt.rank_aggregate.by_sample(
    adata,
    groupby='cell_type',
    min_cells=10,
    use_raw=False
)
```

---

### Lecture 9: Spatial Analysis

**Key Concepts**:
- Spatial coordinates and neighborhoods
- Spatially variable genes (SVGs)
- Neighborhood enrichment
- Spatial domain identification

**Common Mistakes**:
- Not accounting for tissue structure
- Confusing spatial and expression similarity
- Ignoring batch effects in spatial data

**Best Practices**:
```python
# Good: Spatial analysis workflow
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata, mode='moran')
sq.gr.nhood_enrichment(adata, cluster_key='clusters')
```

---

### Lecture 10-11: Deep Learning

**Key Concepts**:
- VAE architecture (scVI)
- Semi-supervised learning (scANVI)
- Latent space interpretation
- Integration quality assessment

**Common Mistakes**:
- Insufficient training epochs
- Not checking model convergence
- Over-fitting to batch structure

**Best Practices**:
```python
# Good: scVI training
scvi.model.SCVI.setup_anndata(adata, layer="counts")
model = scvi.model.SCVI(adata, n_latent=30)
model.train(max_epochs=400, early_stopping=True)

# Check convergence
model.history["elbo_train"].plot()
```

---

### Lecture 12: Neurogenomics Project

**Key Concepts**:
- Complete analysis workflow
- Brain cell type identification
- Biological interpretation
- Reporting findings

**Common Mistakes**:
- Skipping quality control steps
- Not validating cell type assignments
- Missing biological context

**Best Practices**:
- Follow systematic workflow from Lectures 5-6
- Use known brain markers for validation
- Interpret results in biological context
- Document all analysis decisions

---

## Grading Rubrics

### Code Quality (50%)
- **Excellent (45-50)**: Clean, well-commented, follows best practices
- **Good (40-44)**: Works correctly, minor style issues
- **Adequate (35-39)**: Works but inefficient or poorly organized
- **Poor (<35)**: Errors, incomplete, or incorrect results

### Completeness (30%)
- **Excellent (27-30)**: All tasks complete with detailed outputs
- **Good (24-26)**: All tasks complete, outputs adequate
- **Adequate (21-23)**: Most tasks complete
- **Poor (<21)**: Multiple tasks incomplete

### Interpretation (20%)
- **Excellent (18-20)**: Deep biological insight, correct interpretation
- **Good (16-17)**: Correct interpretation, adequate explanation
- **Adequate (14-15)**: Basic interpretation present
- **Poor (<14)**: Missing or incorrect interpretation

---

## Common Technical Issues

### Issue: Package Import Errors
```python
# Solution: Check installation
import sys
!{sys.executable} -m pip install scanpy --upgrade
```

### Issue: Memory Errors
```python
# Solution: Subsample data
adata = adata[np.random.choice(adata.n_obs, 5000, replace=False)]
```

### Issue: Plotting Problems
```python
# Solution: Set matplotlib parameters
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 100
```

---

## Additional Resources for Solutions

### Key Papers
1. **Luecken & Theis (2019)**: Current best practices in scRNA-seq analysis
2. **Gayoso et al. (2022)**: A Python library for probabilistic analysis of scRNA-seq data
3. **Bergen et al. (2020)**: Generalizing RNA velocity to transient cell states

### Advanced Topics Not Covered
- Multi-modal integration (CITE-seq, ATAC-seq)
- Spatial deconvolution methods
- Perturbation response modeling
- Large-scale atlas integration

---

## Solution Updates

Solutions are regularly updated to:
- Incorporate latest package versions
- Add alternative approaches
- Improve explanations
- Fix any errors

Check the repository for the latest versions.

---

**For instructors**: Contact [email] for editable solution files
**For students**: Use solutions responsibly to enhance learning

Last updated: 2025-11-02
