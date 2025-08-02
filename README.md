# 🧬 Peptide Tag Generator

Convert any text into a unique peptide sequence and visualize its 2D molecular structure!

## ✨ Features

- **Text to Peptide Conversion** - Transform names, words, or any text into amino acid sequences
- **2D Molecular Visualization** - Interactive molecular structures using RDKit
- **Protein Properties** - Calculate molecular weight, isoelectric point, and more
- **Function Analysis** - AI-powered peptide function prediction and database searches
- **Downloadable SVGs** - High-quality vector graphics for presentations and publications
- **Scientific APIs** - Integration with PubChem, UniProt, PubMed, and specialized peptide databases

## 🚀 Live Demo

**[Deploy your own instance](https://railway.app) or visit the live demo!**

## 🛠️ Local Development

### Prerequisites
- Python 3.11+
- pip

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/your-username/peptide-tag.git
   cd peptide-tag
   ```

2. **Create virtual environment**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Run the application**
   ```bash
   python app.py
   ```

5. **Open your browser**
   Visit `http://localhost:5000`

## 🧪 How It Works

### Text Conversion
- Input text is converted to uppercase
- Each letter is mapped to amino acids using the standard genetic code
- Ambiguous amino acids (B, Z, X, J) are resolved to common alternatives

### Molecular Structure
- **RDKit** generates 2D molecular structures from SMILES notation
- **PubChem API** provides accurate molecular data
- **Interactive visualization** with copyable SMILES strings

### Function Analysis
- **AI-powered analysis** of peptide composition and patterns
- **Database searches** across PubChem, UniProt, PubMed
- **Specialized peptide databases** (APD3, CAMP)
- **Known peptide matching** for functional peptides

## 📊 Example Output

Input: `"ROOSA"` → Output: `"ROOSA"` → Peptide: `Arg-Gln-Gln-Ser-Ala`

**Properties:**
- Molecular Weight: 573.6 Da
- Isoelectric Point: 10.2
- Predicted Function: Cell-penetrating peptide

## 🏗️ Architecture

- **Flask** - Web framework
- **Biopython** - Protein analysis and properties
- **RDKit** - Molecular cheminformatics and 2D structures
- **PubChem API** - Chemical database integration
- **UniProt API** - Protein database searches
- **SVG Generation** - Vector graphics for downloads

## 🚀 Deployment

This app is ready for deployment on:
- **Railway** (recommended) - Free tier, 5-minute setup
- **Render** - Free tier, good performance
- **Heroku** - $7/month, established platform
- **PythonAnywhere** - Free tier, Python-focused

See [DEPLOYMENT.md](DEPLOYMENT.md) for detailed instructions.

## 🔬 Scientific Features

### Molecular Analysis
- Complete 2D molecular structure visualization
- Amino acid composition analysis
- Secondary structure prediction
- Biological pathway identification

### Database Integration
- **PubChem** - Chemical structures and properties
- **UniProt** - Protein function and organism data
- **PubMed** - Scientific literature search
- **Specialized databases** - Antimicrobial peptides, therapeutic peptides

### AI-Powered Analysis
- Sequence pattern recognition
- Function prediction based on composition
- Therapeutic potential assessment
- Structural feature identification

## 📁 Project Structure

```
peptide-tag/
├── app.py                 # Main Flask application
├── requirements.txt       # Python dependencies
├── Procfile              # Deployment configuration
├── runtime.txt           # Python version specification
├── templates/            # HTML templates
│   ├── index.html       # Landing page
│   └── result.html      # Results page
├── static/              # Static files
│   └── style.css        # Styling
└── DEPLOYMENT.md        # Deployment guide
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is open source and available under the [MIT License](LICENSE).

## 🙏 Acknowledgments

- **RDKit** - Molecular cheminformatics toolkit
- **Biopython** - Biological computation library
- **PubChem** - Chemical database
- **UniProt** - Protein database
- **Flask** - Web framework

## 📞 Contact

Made by: [@https://robin-gustafsson.com/](https://robin-gustafsson.com/)

---

**Turn any text into a protein! 🧬✨** 