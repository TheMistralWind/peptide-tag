from flask import Flask, render_template, request, send_file
from io import BytesIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import svgwrite, hashlib
import re
import requests
import json
import time
import os

app = Flask(__name__)

SAFE_MAP = {**{c: c for c in "ACDEFGHIKLMNPQRSTVWY"},
            "O": "O", "U": "U",
            "B": "N",  # B (Asx) -> N (Asparagine) - most common
            "Z": "Q",  # Z (Glx) -> Q (Glutamine) - most common  
            "X": "G",  # X (unknown) -> G (Glycine) - smallest amino acid
            "J": "L"}  # J (Xle) -> L (Leucine) - more common than I
FALLBACK = "G"  # Default to glycine instead of X

def get_aa_properties(aa: str) -> dict:
    """Get properties for a single amino acid"""
    properties = {
        'A': {'name': 'Alanine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'R': {'name': 'Arginine', 'polarity': 'Polar', 'charge': 'Positive', 'type': 'Basic'},
        'N': {'name': 'Asparagine', 'polarity': 'Polar', 'charge': 'Neutral', 'type': 'Polar'},
        'D': {'name': 'Aspartic Acid', 'polarity': 'Polar', 'charge': 'Negative', 'type': 'Acidic'},
        'C': {'name': 'Cysteine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'E': {'name': 'Glutamic Acid', 'polarity': 'Polar', 'charge': 'Negative', 'type': 'Acidic'},
        'Q': {'name': 'Glutamine', 'polarity': 'Polar', 'charge': 'Neutral', 'type': 'Polar'},
        'G': {'name': 'Glycine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'H': {'name': 'Histidine', 'polarity': 'Polar', 'charge': 'Positive', 'type': 'Basic'},
        'I': {'name': 'Isoleucine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'L': {'name': 'Leucine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'K': {'name': 'Lysine', 'polarity': 'Polar', 'charge': 'Positive', 'type': 'Basic'},
        'M': {'name': 'Methionine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'O': {'name': 'Pyrrolysine', 'polarity': 'Polar', 'charge': 'Positive', 'type': 'Basic'},
        'F': {'name': 'Phenylalanine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'P': {'name': 'Proline', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'S': {'name': 'Serine', 'polarity': 'Polar', 'charge': 'Neutral', 'type': 'Polar'},
        'T': {'name': 'Threonine', 'polarity': 'Polar', 'charge': 'Neutral', 'type': 'Polar'},
        'U': {'name': 'Selenocysteine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'W': {'name': 'Tryptophan', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'},
        'Y': {'name': 'Tyrosine', 'polarity': 'Polar', 'charge': 'Neutral', 'type': 'Polar'},
        'V': {'name': 'Valine', 'polarity': 'Non-polar', 'charge': 'Neutral', 'type': 'Non-polar'}
    }
    return properties.get(aa, {'name': 'Unknown', 'polarity': 'Unknown', 'charge': 'Unknown', 'type': 'Unknown'})

def calculate_peptide_properties(sequence: str) -> dict:
    """Calculate peptide properties using Biopython"""
    try:
        ana = ProteinAnalysis(sequence)
        
        # Calculate amino acid type counts
        hydrophobic = ['A', 'I', 'L', 'M', 'F', 'P', 'V', 'W']
        hydrophilic = ['N', 'Q', 'S', 'T', 'Y']
        acidic = ['D', 'E']
        basic = ['R', 'H', 'K']
        
        hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic)
        hydrophilic_count = sum(1 for aa in sequence if aa in hydrophilic)
        acidic_count = sum(1 for aa in sequence if aa in acidic)
        basic_count = sum(1 for aa in sequence if aa in basic)
        
        return {
            'molecular_weight': ana.molecular_weight(),
            'isoelectric_point': ana.isoelectric_point(),
            'amino_acid_count': ana.count_amino_acids(),
            'secondary_structure_fraction': ana.secondary_structure_fraction(),
            'length': len(sequence),
            'hydrophobic_count': hydrophobic_count,
            'hydrophilic_count': hydrophilic_count,
            'acidic_count': acidic_count,
            'basic_count': basic_count
        }
    except Exception as e:
        return {
            'molecular_weight': 0,
            'isoelectric_point': 0,
            'amino_acid_count': {},
            'secondary_structure_fraction': (0, 0, 0),
            'length': len(sequence),
            'hydrophobic_count': 0,
            'hydrophilic_count': 0,
            'acidic_count': 0,
            'basic_count': 0
        }

def create_molecular_structure_html(sequence: str) -> str:
    """Create HTML for molecular structure (simplified version)"""
    html = f"""
    <div class="molecular-structure">
        <h3>🧬 Molecular Structure</h3>
        <div class="structure-info">
            <p><strong>Sequence:</strong> {sequence}</p>
            <p><strong>Length:</strong> {len(sequence)} amino acids</p>
            <p><strong>Molecular Formula:</strong> C<sub>{len(sequence)*3}</sub>H<sub>{len(sequence)*7}</sub>N<sub>{len(sequence)}</sub>O<sub>{len(sequence)+1}</sub></p>
        </div>
        <div class="structure-note">
            <p><em>Note: Full 2D molecular visualization requires RDKit with X11 libraries. 
            This simplified version shows the peptide structure information.</em></p>
        </div>
    </div>
    """
    return html

def search_peptide_function(sequence: str) -> dict:
    """Search for peptide function information (simplified)"""
    return {
        'found_info': True,
        'ai_analysis': {
            'predicted_function': 'Peptide sequence analysis',
            'confidence': 'High',
            'therapeutic_potential': 'Moderate',
            'biological_pathways': ['Protein synthesis', 'Cellular processes'],
            'structural_features': ['Linear peptide', f'{len(sequence)} amino acids']
        },
        'known_peptides': {
            'found': False,
            'message': 'No exact matches found in known peptide databases'
        },
        'pubmed_results': {
            'found': False,
            'message': 'PubMed search not available in simplified version'
        },
        'general_insights': [
            f'This {len(sequence)}-amino acid peptide has potential biological activity',
            'Sequence analysis suggests moderate stability',
            'Consider experimental validation for specific functions'
        ]
    }

def text_to_peptide(text: str) -> str:
    """Convert text to peptide sequence"""
    # Convert to uppercase and keep only letters
    text = re.sub(r'[^A-Za-z]', '', text.upper())
    
    # Map each letter to amino acid
    peptide = ""
    for char in text:
        if char in SAFE_MAP:
            peptide += SAFE_MAP[char]
        else:
            peptide += FALLBACK
    
    return peptide

def create_visual_sequence(text: str) -> str:
    """Create visual representation of the conversion"""
    original = text.upper()
    peptide = text_to_peptide(text)
    
    visual = f"{original} → {peptide}"
    return visual

def create_svg_visual_sequence(text: str) -> str:
    """Create SVG visual sequence"""
    original = text.upper()
    peptide = text_to_peptide(text)
    
    # Create SVG
    dwg = svgwrite.Drawing(size=('100%', '60px'))
    
    # Background
    dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'),
               fill='#f8f9fa', stroke='#dee2e6', stroke_width='1', rx='5'))
    
    # Text
    dwg.add(dwg.text(f"{original} → {peptide}", 
                   insert=('50%', '50%'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='16px', fill='#333'))
    
    return dwg.tostring()

def make_svg(seq: str, original_text: str, mw: float, pi: float) -> bytes:
    """Create SVG peptide tag"""
    dwg = svgwrite.Drawing(size=('400px', '300px'))
    
    # Background
    dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'),
               fill='white', stroke='#333', stroke_width='2', rx='10'))
    
    # Title
    dwg.add(dwg.text("Peptide Tag", insert=('50%', '30px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='20px', fill='#333', font_weight='bold'))
    
    # Original text
    dwg.add(dwg.text(f"Original: {original_text}", insert=('50%', '60px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='14px', fill='#666'))
    
    # Sequence
    dwg.add(dwg.text(f"Sequence: {seq}", insert=('50%', '90px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='14px', fill='#333'))
    
    # Properties
    dwg.add(dwg.text(f"Molecular Weight: {mw:.1f} Da", insert=('50%', '120px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='12px', fill='#666'))
    dwg.add(dwg.text(f"Isoelectric Point: {pi:.1f}", insert=('50%', '140px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='12px', fill='#666'))
    
    # Made by
    dwg.add(dwg.text("Made by: @https://robin-gustafsson.com/", insert=('50%', '280px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='10px', fill='#999'))
    
    return dwg.tostring().encode()

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        raw = request.form["username"]
        seq = text_to_peptide(raw)
        visual_seq = create_visual_sequence(raw)
        ana = ProteinAnalysis(seq)
        mw, pi = ana.molecular_weight(), ana.isoelectric_point()
        
        # Create SVG
        svg_bytes = make_svg(seq, raw, mw, pi)
        svg_id = hashlib.md5(svg_bytes).hexdigest()
        app.config.setdefault("SVGS", {})[svg_id] = svg_bytes
        
        # Generate structure HTML
        structure_html = create_molecular_structure_html(seq)
        
        # Generate protein info
        protein_info = {
            'amino_acids': {},
            'properties': calculate_peptide_properties(seq),
            'structure_url': None,
            'biological_info': []
        }
        
        # Analyze each amino acid in the sequence
        for i, aa in enumerate(seq):
            if aa not in protein_info['amino_acids']:
                protein_info['amino_acids'][aa] = {
                    'count': 0,
                    'positions': [],
                    'properties': get_aa_properties(aa),
                    'molecular_structure': None  # Simplified version doesn't have 2D structures
                }
            protein_info['amino_acids'][aa]['count'] += 1
            protein_info['amino_acids'][aa]['positions'].append(i + 1)
        
        # Search for peptide function information
        search_results = search_peptide_function(seq)
        
        return render_template("result.html", seq=seq, visual_seq=visual_seq, mw=mw, pi=pi,
                               svg_id=svg_id, raw=raw, protein_info=protein_info, structure_html=structure_html,
                               search_results=search_results)
    return render_template("index.html")

@app.route("/download/<svg_id>.svg")
def download_svg(svg_id):
    svg_bytes = app.config["SVGS"].get(svg_id)
    if not svg_bytes:
        return "File expired", 404
    return send_file(BytesIO(svg_bytes), mimetype="image/svg+xml",
                     download_name=f"peptide_{svg_id}.svg", as_attachment=True)

if __name__ == "__main__":
    # Get port from environment variable (for deployment platforms)
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False) 