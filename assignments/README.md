# Single-Cell Neurogenomics Course Assignments

## Overview

This directory contains 12 Jupyter notebook assignments for the Single-Cell Neurogenomics course. Each assignment is designed to be completed in approximately 60 minutes and covers specific learning objectives aligned with the course curriculum.

## Course Structure

### Module 1: Introduction (Lectures 1-4)
- **Lecture 1**: Introduction to Neurogenomics
- **Lecture 2**: Introduction to Single-Cell Technology
- **Lecture 3**: Fundamentals of Python
- **Lecture 4**: Quantification of Single-Cell RNA-seq Data

### Module 2: Core scRNA-seq Analysis (Lectures 5-6)
- **Lecture 5**: Quality Control, Normalization, and Preprocessing
- **Lecture 6**: Downstream Analysis and Batch Correction

### Module 3: Advanced Analysis (Lectures 7-9)
- **Lecture 7**: Trajectory Inference and Fate Probability
- **Lecture 8**: Ligand-Receptor Interactions and Differential Abundance
- **Lecture 9**: Spatial Analysis and Spatial Structure

### Module 4: Deep Learning and Projects (Lectures 10-12)
- **Lecture 10**: Deep Learning in Single-Cell Genomics - Part 1
- **Lecture 11**: Deep Learning in Single-Cell Genomics - Part 2
- **Lecture 12**: Single-Cell Neurogenomics Project

## Assignment Details

Each assignment includes:
- **Learning Objectives**: Clear goals for what students will learn
- **Instructions**: Step-by-step guidance for each task
- **Tasks**: 3-5 tasks with point values totaling 100-110 points
- **Hints**: Technical tips to help students complete tasks
- **Expected Outputs**: Description of what results should look like
- **Reflection Questions**: (Bonus) Conceptual questions for deeper understanding

## Required Software

### Python Environment
```bash
# Create conda environment
conda create -n scverse python=3.9
conda activate scverse

# Install core packages
pip install scanpy scvi-tools scvelo cellrank liana squidpy
pip install numpy pandas matplotlib seaborn jupyter
```

### Package Versions
- Python: 3.9+
- scanpy: 1.9+
- anndata: 0.9+
- scvi-tools: 1.0+
- scvelo: 0.3+
- cellrank: 2.0+
- squidpy: 1.3+
- liana: 0.1+

## Datasets Used

### Built-in Datasets (from scanpy/scvi/squidpy)
- PBMC 3k (10X Genomics)
- PBMC 10k (10X Genomics)
- Pancreas development dataset
- Mouse brain Visium (spatial)
- Allen Brain Atlas data

### Downloading Additional Data
```python
import scanpy as sc

# PBMC datasets
adata = sc.datasets.pbmc3k()
adata = sc.datasets.pbmc68k_reduced()

# Other datasets
adata = sc.datasets.paul15()  # Mouse blood development
```

## Assignment Workflow

### For Students
1. **Download** the assignment notebook
2. **Read** learning objectives and task instructions
3. **Complete** each task by filling in code cells
4. **Run** all cells to generate outputs
5. **Save** notebook with outputs included
6. **Answer** reflection questions (bonus points)
7. **Submit** completed notebook

### Grading Criteria
Each task is graded on:
- **Correctness** (50%): Code runs without errors and produces correct results
- **Completeness** (30%): All required components are included
- **Interpretation** (20%): Demonstrates understanding through visualizations and explanations

## Tips for Success

### General Best Practices
1. **Read instructions carefully** before starting each task
2. **Use hints** provided in each task
3. **Consult documentation** for scanpy, scvi-tools, etc.
4. **Visualize results** to check if they make biological sense
5. **Comment your code** to show understanding
6. **Ask questions** if stuck

### Common Issues and Solutions

#### Issue: ModuleNotFoundError
```bash
# Solution: Install missing package
pip install <package_name>
```

#### Issue: Memory errors with large datasets
```python
# Solution: Use subsampling
adata = adata[np.random.choice(adata.n_obs, 5000, replace=False), :]
```

#### Issue: Plotting doesn't work
```python
# Solution: Ensure matplotlib backend is set
import matplotlib
matplotlib.use('Agg')  # For non-interactive
# Or
%matplotlib inline  # For Jupyter notebooks
```

## Learning Resources

### Official Documentation
- **Scanpy**: https://scanpy.readthedocs.io/
- **scvi-tools**: https://docs.scvi-tools.org/
- **CellRank**: https://cellrank.readthedocs.io/
- **Squidpy**: https://squidpy.readthedocs.io/
- **scVelo**: https://scvelo.readthedocs.io/

### Tutorials
- **Scanpy Tutorial**: https://scanpy-tutorials.readthedocs.io/
- **scvi-tools Tutorial**: https://docs.scvi-tools.org/en/stable/tutorials/
- **Best Practices**: Luecken & Theis (2019) - Current best practices in scRNA-seq

### Community Support
- **Scanpy Discourse**: https://discourse.scverse.org/
- **scvi-tools Discussions**: https://github.com/scverse/scvi-tools/discussions
- **Biostars**: https://www.biostars.org/

## Assignment Schedule

| Lecture | Date | Assignment | Topics |
|---------|------|------------|--------|
| 1 | Dec 5, 2025 | lecture01_neurogenomics.ipynb | Gene expression, brain regions |
| 2 | Dec 6, 2025 | lecture02_single_cell_technology.ipynb | AnnData, sparsity, PBMC data |
| 3 | Dec 12, 2025 | lecture03_python_fundamentals.ipynb | Python basics, functions |
| 4 | Dec 13, 2025 | lecture04_quantification.ipynb | Count matrices, QC metrics |
| 5 | Dec 19, 2025 | lecture05_qc_preprocessing.ipynb | Filtering, normalization, PCA |
| 6 | Dec 20, 2025 | lecture06_clustering_integration.ipynb | Clustering, markers, batch correction |
| 7 | Dec 26, 2025 | lecture07_trajectory_inference.ipynb | RNA velocity, CellRank |
| 8 | Dec 27, 2025 | lecture08_cell_communication.ipynb | Ligand-receptor, LIANA |
| 9 | Jan 2, 2026 | lecture09_spatial_analysis.ipynb | Spatial transcriptomics, Squidpy |
| 10 | Jan 3, 2026 | lecture10_deep_learning_part1.ipynb | scVI, VAE, latent space |
| 11 | Jan 9, 2026 | lecture11_deep_learning_part2.ipynb | scANVI, integration |
| 12 | Jan 10, 2026 | lecture12_neurogenomics_project.ipynb | Brain cell types, final project |

## Submission Guidelines

### File Naming
```
<studentID>_<lectureNumber>_<lastName>.ipynb

Example: 12345_lecture05_smith.ipynb
```

### What to Submit
1. Completed Jupyter notebook (.ipynb)
2. All code cells executed with outputs visible
3. All plots and visualizations generated
4. Reflection questions answered (if applicable)

### Submission Checklist
- [ ] All tasks completed
- [ ] All code cells run without errors
- [ ] All plots properly labeled and visible
- [ ] Reflection questions answered
- [ ] Notebook saved with outputs
- [ ] File named correctly

## Getting Help

### During Class
- Ask instructor or TA during lecture
- Use office hours for detailed questions

### Online Resources
- Course discussion forum
- Slack/Discord channel (if available)
- Email instructor: [contact info]

### Debugging Tips
1. Read error messages carefully
2. Check variable types and shapes
3. Print intermediate results
4. Use Google/Stack Overflow for common errors
5. Check package documentation

## Academic Integrity

- You may discuss concepts with classmates
- You may consult documentation and online resources
- You **must** write your own code
- You **must** generate your own outputs
- Copying code from others is prohibited
- Cite any external resources used

## Additional Notes

### Computational Requirements
- **RAM**: 8GB minimum, 16GB recommended
- **Storage**: 10GB free space
- **Processor**: Modern multi-core CPU
- **GPU**: Optional, helpful for deep learning lectures

### Time Management
- Start assignments early
- Don't wait until the last minute
- Budget 60-90 minutes per assignment
- Seek help if stuck for >30 minutes

### Extensions and Late Policy
- Refer to course syllabus
- Contact instructor for extensions
- Late penalties may apply

---

**Good luck with your assignments!**

For questions about assignments, contact: [instructor email]

Last updated: 2025-11-02
