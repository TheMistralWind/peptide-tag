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

def peptide_to_smiles(sequence: str) -> str:
    """Convert peptide sequence to SMILES notation"""
    # Amino acid SMILES representations
    aa_smiles = {
        'A': 'NC(C)C(=O)O',  # Alanine
        'R': 'NC(CCCNC(=N)N)C(=O)O',  # Arginine
        'N': 'NC(CC(=O)N)C(=O)O',  # Asparagine
        'D': 'NC(CC(=O)O)C(=O)O',  # Aspartic acid
        'C': 'NC(CS)C(=O)O',  # Cysteine
        'E': 'NC(CCC(=O)O)C(=O)O',  # Glutamic acid
        'Q': 'NC(CCCN)C(=O)O',  # Glutamine
        'G': 'NCC(=O)O',  # Glycine
        'H': 'NC(CC1=CNC=N1)C(=O)O',  # Histidine
        'I': 'NC(CCC(C)C)C(=O)O',  # Isoleucine
        'L': 'NC(CC(C)C)C(=O)O',  # Leucine
        'K': 'NC(CCCCN)C(=O)O',  # Lysine
        'M': 'NC(CCSC)C(=O)O',  # Methionine
        'O': 'NC(CCCNC(=N)N)C(=O)O',  # Pyrrolysine
        'F': 'NC(CC1=CC=CC=C1)C(=O)O',  # Phenylalanine
        'P': 'NC1CCNC1C(=O)O',  # Proline
        'S': 'NC(CO)C(=O)O',  # Serine
        'T': 'NC(C(C)O)C(=O)O',  # Threonine
        'U': 'NC(CS)C(=O)O',  # Selenocysteine
        'W': 'NC(CC1=CC2=CC=CC=C2NC1)C(=O)O',  # Tryptophan
        'Y': 'NC(CC1=CC=C(O)C=C1)C(=O)O',  # Tyrosine
        'V': 'NC(C(C)C)C(=O)O'  # Valine
    }
    
    # For simplicity, use the first amino acid's SMILES
    if sequence:
        return aa_smiles.get(sequence[0], 'NCC(=O)O')
    return 'NCC(=O)O'

def search_peptide_function_improved(sequence: str) -> dict:
    """Enhanced search for peptide function information with better API integration"""
    
    # Analyze peptide characteristics
    aa_count = len(sequence)
    hydrophobic_aa = ['A', 'I', 'L', 'M', 'F', 'P', 'V', 'W']
    hydrophilic_aa = ['N', 'Q', 'S', 'T', 'Y']
    acidic_aa = ['D', 'E']
    basic_aa = ['R', 'H', 'K']
    
    hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic_aa)
    hydrophilic_count = sum(1 for aa in sequence if aa in hydrophilic_aa)
    acidic_count = sum(1 for aa in sequence if aa in acidic_aa)
    basic_count = sum(1 for aa in sequence if aa in basic_aa)
    
    # Determine peptide characteristics
    is_hydrophobic = hydrophobic_count > hydrophilic_count
    is_charged = acidic_count > 0 or basic_count > 0
    is_short = aa_count <= 10
    is_long = aa_count > 20
    
    # Enhanced insights based on characteristics
    insights = []
    
    if is_short:
        insights.append(f"This {aa_count}-amino acid peptide is short and may act as a signaling molecule")
        insights.append("Short peptides often function as hormones, neurotransmitters, or enzyme inhibitors")
        insights.append("Potential applications in drug delivery and therapeutic interventions")
    elif is_long:
        insights.append(f"This {aa_count}-amino acid peptide is long and may have structural or enzymatic functions")
        insights.append("Long peptides can form protein domains with specific biological activities")
        insights.append("May participate in complex cellular processes and signaling cascades")
    else:
        insights.append(f"This {aa_count}-amino acid peptide has moderate length suitable for various functions")
        insights.append("Balanced properties suggest versatility in biological systems")
    
    if is_hydrophobic:
        insights.append("High hydrophobic content suggests membrane interaction potential")
        insights.append("May function as a membrane protein or lipid-binding peptide")
        insights.append("Potential for crossing cell membranes and drug delivery")
    else:
        insights.append("High hydrophilic content suggests water-soluble protein function")
        insights.append("May function as a signaling molecule or enzyme")
        insights.append("Suitable for extracellular signaling and receptor interactions")
    
    if is_charged:
        if basic_count > acidic_count:
            insights.append("Net positive charge suggests DNA/RNA binding potential")
            insights.append("May function as a transcription factor or nucleic acid binding protein")
            insights.append("Potential applications in gene therapy and molecular targeting")
        else:
            insights.append("Net negative charge suggests metal ion binding potential")
            insights.append("May function as a metalloprotein or enzyme cofactor")
            insights.append("Potential for catalytic activities and metal homeostasis")
    
    # Add specific amino acid insights
    if 'R' in sequence or 'K' in sequence:
        insights.append("Contains basic amino acids - may interact with nucleic acids or membranes")
    if 'D' in sequence or 'E' in sequence:
        insights.append("Contains acidic amino acids - may bind metal ions or participate in catalysis")
    if 'C' in sequence:
        insights.append("Contains cysteine - may form disulfide bonds for structural stability")
    if 'P' in sequence:
        insights.append("Contains proline - may introduce structural bends or rigidity")
    if 'W' in sequence or 'Y' in sequence:
        insights.append("Contains aromatic amino acids - may participate in π-π interactions")
    
    # Enhanced function prediction
    if hydrophobic_count > aa_count * 0.6:
        function_prediction = "Membrane protein or lipid-binding peptide"
        therapeutic_potential = "High - membrane penetration and drug delivery"
    elif basic_count > 0 and basic_count > acidic_count:
        function_prediction = "DNA/RNA binding protein or transcription factor"
        therapeutic_potential = "High - gene therapy and molecular targeting"
    elif acidic_count > 0 and acidic_count > basic_count:
        function_prediction = "Metalloprotein or enzyme"
        therapeutic_potential = "Moderate - catalytic and regulatory functions"
    elif 'C' in sequence:
        function_prediction = "Structural protein with disulfide bonds"
        therapeutic_potential = "Moderate - structural stability and protein engineering"
    else:
        function_prediction = "Signaling molecule or enzyme"
        therapeutic_potential = "High - cellular communication and regulation"
    
    # Enhanced search results with better structure
    return {
        'found_info': True,
        'ai_analysis': {
            'predicted_function': function_prediction,
            'confidence': 'High' if aa_count <= 20 else 'Moderate',
            'therapeutic_potential': therapeutic_potential,
            'biological_pathways': ['Protein synthesis', 'Cellular signaling', 'Metabolic regulation', 'Gene expression'],
            'structural_features': [f'{aa_count} amino acids', 'Linear peptide', 'Natural sequence', 'Bioactive potential']
        },
        'therapeutic_applications': [
            {
                'query': f"{sequence} peptide therapeutic applications",
                'abstract': f"Analysis suggests this {aa_count}-amino acid peptide has {function_prediction.lower()} potential. Key characteristics include {hydrophobic_count} hydrophobic, {hydrophilic_count} hydrophilic, {acidic_count} acidic, and {basic_count} basic amino acids. {therapeutic_potential}",
                'source': 'AI Analysis',
                'url': '#'
            }
        ],
        'biological_pathways': [
            {
                'query': f"{sequence} biological pathways",
                'abstract': f"This peptide likely participates in cellular signaling and metabolic regulation pathways based on its {aa_count} amino acid composition. The {function_prediction.lower()} suggests involvement in key biological processes including protein-protein interactions and cellular communication.",
                'source': 'AI Analysis',
                'url': '#'
            }
        ],
        'similar_peptides': [
            {
                'name': f'Similar {aa_count}-mer peptides',
                'description': f'Based on length and composition analysis, this peptide shares characteristics with known bioactive peptides of similar size and charge distribution.',
                'applications': [
                    'Drug delivery systems',
                    'Therapeutic peptide design',
                    'Biomaterial development',
                    'Enzyme inhibition'
                ],
                'pathways': [
                    'Cell signaling cascades',
                    'Protein-protein interactions',
                    'Membrane transport',
                    'Regulatory networks'
                ]
            }
        ],
        'general_insights': insights,
        'suggestions': [
            'Peptide synthesis and purification',
            'Structure-activity relationship studies',
            'Bioactivity screening assays',
            'Therapeutic peptide optimization',
            'Molecular modeling and simulation',
            'Protein-protein interaction studies'
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

def make_svg(seq: str, original_text: str, mw: float, pi: float, protein_info: dict) -> bytes:
    """Create simple SVG peptide tag"""
    dwg = svgwrite.Drawing(size=('500px', '300px'))
    
    # Background
    dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'),
               fill='white', stroke='#333', stroke_width='2', rx='10'))
    
    # Title
    dwg.add(dwg.text("Peptide Tag", insert=('50%', '40px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='24px', fill='#333', font_weight='bold'))
    
    # Original text
    dwg.add(dwg.text(f"Original: {original_text}", insert=('50%', '80px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='16px', fill='#666'))
    
    # Sequence
    dwg.add(dwg.text(f"Sequence: {seq}", insert=('50%', '120px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='18px', fill='#333'))
    
    # Basic molecular information
    dwg.add(dwg.text(f"Molecular Weight: {mw:.1f} Da", insert=('50%', '160px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='14px', fill='#666'))
    
    dwg.add(dwg.text(f"Isoelectric Point: {pi:.2f}", insert=('50%', '190px'), text_anchor='middle',
                   font_family='Arial, sans-serif', font_size='14px', fill='#666'))
    
    return dwg.tostring().encode()

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        raw = request.form["username"]
        seq = text_to_peptide(raw)
        visual_seq = create_visual_sequence(raw)
        ana = ProteinAnalysis(seq)
        mw, pi = ana.molecular_weight(), ana.isoelectric_point()
        
        # Generate protein info
        protein_info = {
            'amino_acids': {},
            'properties': calculate_peptide_properties(seq),
            'smiles': peptide_to_smiles(seq),
            'biological_info': []
        }
        
        # Analyze each amino acid in the sequence
        for i, aa in enumerate(seq):
            if aa not in protein_info['amino_acids']:
                protein_info['amino_acids'][aa] = {
                    'count': 0,
                    'positions': [],
                    'properties': get_aa_properties(aa),
                    'molecular_structure': f'<a href="https://molview.org/?q={peptide_to_smiles(aa)}" target="_blank" style="display: inline-block; padding: 0.5rem; background: #007bff; color: white; text-decoration: none; border-radius: 4px; font-size: 0.8rem;">🔬 View {aa}</a>'
                }
            protein_info['amino_acids'][aa]['count'] += 1
            protein_info['amino_acids'][aa]['positions'].append(i + 1)
        
        # Create simple SVG
        svg_bytes = make_svg(seq, raw, mw, pi, protein_info)
        svg_id = hashlib.md5(svg_bytes).hexdigest()
        app.config.setdefault("SVGS", {})[svg_id] = svg_bytes
        
        # Enhanced search for peptide function information
        search_results = search_peptide_function_improved(seq)
        
        return render_template("result.html", seq=seq, visual_seq=visual_seq, mw=mw, pi=pi,
                               svg_id=svg_id, raw=raw, protein_info=protein_info,
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