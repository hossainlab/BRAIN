# Bioinformatics Research with AI in Neurosciences (BRAIN)

**Bioinformatics Research with AI in Neurosciences (BRAIN)** is a comprehensive research-oriented training program focusing on the intersection of Neurogenomics, Bioinformatics, and Artificial Intelligence. Participants will explore data-driven neuroscience—covering genomics, transcriptomics, and machine learning methods applied to brain diseases and neurological research.

## Course Website

This repository contains the source code for the BRAIN course website built with [Quarto](https://quarto.org/).

### Features

- **Professional Design**: Neuroscience-inspired color scheme (purples, pinks, neural network aesthetic)
- **16 Comprehensive Modules**: From neurogenomics basics to advanced deep learning
- **Interactive Navigation**: Top navbar with dropdown menu for all modules
- **Responsive Layout**: Works on desktop, tablet, and mobile devices
- **Search Functionality**: Built-in search across all content

## Course Modules

1. **Introduction to Neurogenomics** - Understanding the genetic basis of brain function
2. **Introduction to Single-Cell Technology** - Exploring cellular heterogeneity
3. **Linux Command Line Fundamentals** - Essential Unix skills for bioinformatics
4. **CellRanger Pipeline** - Processing 10x Genomics single-cell data
5. **Fundamentals of Python: Part 1** - Introduction to Python programming
6. **Fundamentals of Python: Part 2** - Advanced Python and scientific computing
7. **Quality Control** - Ensuring data quality in single-cell analysis
8. **Dimensionality Reduction** - Visualizing high-dimensional data
9. **Data Integration** - Combining multiple single-cell datasets
10. **Clustering** - Identifying cell populations
11. **Differential Expression** - Finding genes that define cell types
12. **Cell Type Prediction** - Automated cell type annotation
13. **Trajectory Inference** - Analyzing cellular differentiation
14. **Spatial Transcriptomics** - Mapping gene expression in tissue context
15. **Deep Learning: Part 1** - Neural networks and autoencoders
16. **Deep Learning: Part 2** - Graph neural networks and transformers


## Setup and Installation

### Prerequisites

1. **Install Quarto**
   - Download from: https://quarto.org/docs/get-started/
   - Or use package managers:
     ```bash
     # macOS
     brew install quarto

     # Linux (Ubuntu/Debian)
     sudo apt install quarto-cli
     ```

2. **Clone the Repository**
   ```bash
   git clone https://github.com/yourusername/BRAIN.git
   cd BRAIN
   ```

### Building the Website

#### Preview Locally

```bash
quarto preview
```

This will start a local server (usually at http://localhost:4200) with live reload.

#### Render the Website

```bash
quarto render
```

This generates the static website in the `docs/` directory.

## Customization

### Update Configuration

Edit `_quarto.yml` to customize:
- Website title and description
- Repository URL
- Navigation structure
- Color theme

### Add Your Logo

1. Place logo files in `assets/` directory:
   - `logo.png` - Main logo (recommended: 200x200px)
   - `favicon.png` - Browser icon (32x32px)

2. Update paths in `_quarto.yml`:
   ```yaml
   website:
     logo: assets/logo.png
     favicon: assets/favicon.png
   ```

### Customize Colors

Edit `styles.css` to change the color scheme. Current theme uses:
- Primary purple: `#8B5CF6`
- Accent pink: `#EC4899`
- Neural blue: `#60A5FA`

### Update Content

- **Instructor Information**: Edit `about.qmd`
- **Schedule/Dates**: Edit `schedule.qmd`
- **Module Content**: Edit individual files in `modules/`

## Deployment

### GitHub Pages

1. **Configure GitHub Pages**
   - Go to repository Settings → Pages
   - Source: Deploy from a branch
   - Branch: `main` → `/docs` folder
   - Save

2. **Update `_quarto.yml`**
   ```yaml
   website:
     site-url: https://yourusername.github.io/BRAIN
     repo-url: https://github.com/yourusername/BRAIN
   ```

3. **Render and Push**
   ```bash
   quarto render
   git add .
   git commit -m "Update website"
   git push
   ```

4. **Access Your Site**
   - URL: `https://yourusername.github.io/BRAIN`

### Netlify

1. **Connect Repository**
   - Go to https://netlify.com
   - Import from Git

2. **Build Settings**
   - Build command: `quarto render`
   - Publish directory: `docs`

3. **Deploy**
   - Netlify will automatically build and deploy

### Custom Domain

Add a `CNAME` file in the `docs/` directory:
```
brain-course.yourdomain.com
```

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## Contact

For questions or feedback, please contact:
- Email: [your-email@example.com]
- GitHub Issues: [Repository Issues](https://github.com/yourusername/BRAIN/issues)

## Acknowledgments

- Built with [Quarto](https://quarto.org/)
- Inspired by modern single-cell genomics courses
- Neuroscience-themed design

---

**Note**: Replace placeholder URLs and contact information with your actual details before deployment.