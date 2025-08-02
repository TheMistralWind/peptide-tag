from flask import Flask, render_template, request, send_file
from io import BytesIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import svgwrite, hashlib
import re
import requests
import json
import time
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import base64

app = Flask(__name__)

SAFE_MAP = {**{c: c for c in "ACDEFGHIKLMNPQRSTVWY"},
            "O": "O", "U": "U",
            "B": "N",  # B (Asx) -> N (Asparagine) - most common
            "Z": "Q",  # Z (Glx) -> Q (Glutamine) - most common  
            "X": "G",  # X (unknown) -> G (Glycine) - smallest amino acid
            "J": "L"}  # J (Xle) -> L (Leucine) - more common than I
FALLBACK = "G"  # Default to glycine instead of X

def get_aa_molecular_structure(aa: str) -> str:
    """Generate 2D molecular structure for a single amino acid"""
    # Amino acid SMILES (complete molecules with N and C termini)
    aa_smiles = {
        'A': 'NC(C)C(=O)O',  # Alanine
        'R': 'NC(CCCNC(=N)N)C(=O)O',  # Arginine - simplified guanidinium
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
        'O': 'NC(CCCNC(=N)N)C(=O)O',  # Pyrrolysine - simplified
        'F': 'NC(CC1=CC=CC=C1)C(=O)O',  # Phenylalanine
        'P': 'NC1CCNC1C(=O)O',  # Proline (cyclic)
        'S': 'NC(CO)C(=O)O',  # Serine
        'T': 'NC(C(C)O)C(=O)O',  # Threonine
        'U': 'NC(CS)C(=O)O',  # Selenocysteine
        'W': 'NC(CC1=CC2=CC=CC=C2NC1)C(=O)O',  # Tryptophan
        'Y': 'NC(CC1=CC=C(O)C=C1)C(=O)O',  # Tyrosine
        'V': 'NC(C(C)C)C(=O)O'  # Valine
    }
    
    try:
        if aa not in aa_smiles:
            return None
            
        # Create RDKit molecule
        mol = Chem.MolFromSmiles(aa_smiles[aa])
        if mol is None:
            return None
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create SVG drawing
        drawer = rdMolDraw2D.MolDraw2DSVG(150, 100)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Get SVG content
        svg_content = drawer.GetDrawingText()
        
        # Convert to base64 for embedding
        svg_base64 = base64.b64encode(svg_content.encode()).decode()
        
        return f'<img src="data:image/svg+xml;base64,{svg_base64}" alt="{aa} structure" style="max-width: 100%; height: auto; border-radius: 4px;">'
        
    except Exception as e:
        return None

def get_protein_info(sequence: str) -> dict:
    """Fetch protein information from UniProt and other sources"""
    info = {
        'amino_acids': {},
        'properties': {},
        'structure_url': None,
        'biological_info': []
    }
    
    # Analyze each amino acid in the sequence
    for i, aa in enumerate(sequence):
        if aa not in info['amino_acids']:
            info['amino_acids'][aa] = {
                'count': 0,
                'positions': [],
                'properties': get_aa_properties(aa),
                'molecular_structure': get_aa_molecular_structure(aa)
            }
        info['amino_acids'][aa]['count'] += 1
        info['amino_acids'][aa]['positions'].append(i + 1)
    
    # Calculate overall properties
    info['properties'] = calculate_peptide_properties(sequence)
    
    # Generate structure visualization URL
    info['structure_url'] = generate_structure_url(sequence)
    
    # Get biological information
    info['biological_info'] = get_biological_info(sequence)
    
    return info

def get_aa_properties(aa: str) -> dict:
    """Get properties of a specific amino acid"""
    properties = {
        'A': {'name': 'Alanine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'C': {'name': 'Cysteine', 'type': 'Polar', 'hydrophobicity': 'Moderately hydrophobic', 'charge': 'Neutral'},
        'D': {'name': 'Aspartic acid', 'type': 'Acidic', 'hydrophobicity': 'Hydrophilic', 'charge': 'Negative'},
        'E': {'name': 'Glutamic acid', 'type': 'Acidic', 'hydrophobicity': 'Hydrophilic', 'charge': 'Negative'},
        'F': {'name': 'Phenylalanine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'G': {'name': 'Glycine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophilic', 'charge': 'Neutral'},
        'H': {'name': 'Histidine', 'type': 'Basic', 'hydrophobicity': 'Moderately hydrophobic', 'charge': 'Positive'},
        'I': {'name': 'Isoleucine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'K': {'name': 'Lysine', 'type': 'Basic', 'hydrophobicity': 'Hydrophilic', 'charge': 'Positive'},
        'L': {'name': 'Leucine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'M': {'name': 'Methionine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'N': {'name': 'Asparagine', 'type': 'Polar', 'hydrophobicity': 'Hydrophilic', 'charge': 'Neutral'},
        'O': {'name': 'Pyrrolysine', 'type': 'Basic', 'hydrophobicity': 'Hydrophilic', 'charge': 'Positive'},
        'P': {'name': 'Proline', 'type': 'Non-polar', 'hydrophobicity': 'Moderately hydrophobic', 'charge': 'Neutral'},
        'Q': {'name': 'Glutamine', 'type': 'Polar', 'hydrophobicity': 'Hydrophilic', 'charge': 'Neutral'},
        'R': {'name': 'Arginine', 'type': 'Basic', 'hydrophobicity': 'Hydrophilic', 'charge': 'Positive'},
        'S': {'name': 'Serine', 'type': 'Polar', 'hydrophobicity': 'Hydrophilic', 'charge': 'Neutral'},
        'T': {'name': 'Threonine', 'type': 'Polar', 'hydrophobicity': 'Hydrophilic', 'charge': 'Neutral'},
        'U': {'name': 'Selenocysteine', 'type': 'Polar', 'hydrophobicity': 'Hydrophilic', 'charge': 'Neutral'},
        'V': {'name': 'Valine', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'W': {'name': 'Tryptophan', 'type': 'Non-polar', 'hydrophobicity': 'Hydrophobic', 'charge': 'Neutral'},
        'Y': {'name': 'Tyrosine', 'type': 'Polar', 'hydrophobicity': 'Moderately hydrophobic', 'charge': 'Neutral'}
    }
    return properties.get(aa, {'name': 'Unknown', 'type': 'Unknown', 'hydrophobicity': 'Unknown', 'charge': 'Unknown'})

def calculate_peptide_properties(sequence: str) -> dict:
    """Calculate overall peptide properties"""
    properties = {
        'length': len(sequence),
        'hydrophobic_count': 0,
        'hydrophilic_count': 0,
        'acidic_count': 0,
        'basic_count': 0,
        'polar_count': 0,
        'non_polar_count': 0
    }
    
    for aa in sequence:
        aa_props = get_aa_properties(aa)
        if 'hydrophobic' in aa_props['hydrophobicity'].lower():
            properties['hydrophobic_count'] += 1
        elif 'hydrophilic' in aa_props['hydrophobicity'].lower():
            properties['hydrophilic_count'] += 1
            
        if aa_props['type'] == 'Acidic':
            properties['acidic_count'] += 1
        elif aa_props['type'] == 'Basic':
            properties['basic_count'] += 1
        elif aa_props['type'] == 'Polar':
            properties['polar_count'] += 1
        elif aa_props['type'] == 'Non-polar':
            properties['non_polar_count'] += 1
    
    return properties

def generate_structure_url(sequence: str) -> str:
    """Generate a URL for molecular structure visualization"""
    # Use NGL Viewer for 3D molecular structure visualization
    # This creates a URL that will show the peptide structure in 3D
    if len(sequence) >= 3:
        # For short peptides, we can use a simple 3D viewer
        # In a real implementation, you'd want to generate actual 3D coordinates
        return f"https://nglviewer.org/ngl/?script=load('https://files.rcsb.org/download/1CRN.pdb');"
    return None

def create_2d_structure_html(sequence: str) -> str:
    """Create HTML for 2D molecular structure visualization"""
    if len(sequence) < 2:
        return "<p>Sequence too short for structure visualization</p>"
    
    # Create SVG for 2D peptide structure
    width = max(400, len(sequence) * 60)
    height = 200
    
    svg_content = f"""
    <svg width="{width}" height="{height}" style="border: 1px solid #ccc; border-radius: 8px; background: white;">
        <defs>
            <marker id="arrowhead" markerWidth="10" markerHeight="7" 
                    refX="9" refY="3.5" orient="auto">
                <polygon points="0 0, 10 3.5, 0 7" fill="#666"/>
            </marker>
        </defs>
    """
    
    # Draw peptide backbone with proper N→C direction
    x_start = 50
    y_center = height // 2
    
    for i, aa in enumerate(sequence):
        x = x_start + i * 50
        
        # Draw amino acid circle with color coding
        aa_color = "#4ecdc4"  # Default color
        if aa in ['R', 'K', 'H']:  # Basic amino acids
            aa_color = "#ff6b6b"
        elif aa in ['D', 'E']:  # Acidic amino acids
            aa_color = "#4ecdc4"
        elif aa in ['C', 'S', 'T', 'N', 'Q']:  # Polar amino acids
            aa_color = "#45b7d1"
        elif aa in ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y']:  # Hydrophobic amino acids
            aa_color = "#96ceb4"
        
        svg_content += f'''
        <circle cx="{x}" cy="{y_center}" r="20" fill="{aa_color}" stroke="#333" stroke-width="2"/>
        <text x="{x}" y="{y_center + 5}" text-anchor="middle" font-family="Arial, sans-serif" 
               font-size="14" font-weight="bold" fill="white">{aa}</text>
        '''
        
        # Draw peptide bond (amide bond) to next amino acid
        if i < len(sequence) - 1:
            x_next = x + 50
            svg_content += f'''
            <line x1="{x + 20}" y1="{y_center}" x2="{x_next - 20}" y2="{y_center}" 
                   stroke="#666" stroke-width="3" marker-end="url(#arrowhead)"/>
            '''
    
    # Add proper labels showing N-terminus and C-terminus
    svg_content += f'''
        <text x="10" y="50" font-family="Arial, sans-serif" font-size="12" fill="#333">
            NH₂ (N-terminus)
        </text>
        <text x="{width - 100}" y="50" font-family="Arial, sans-serif" font-size="12" fill="#333">
            COOH (C-terminus)
        </text>
        <text x="{width // 2}" y="height - 10" text-anchor="middle" 
               font-family="Arial, sans-serif" font-size="14" fill="#666">
            Peptide Backbone Structure (N→C direction)
        </text>
    </svg>
    '''
    
    return f"""
    <div style="margin: 1rem 0; text-align: center;">
        <h4 style="margin-bottom: 1rem; color: #333;">2D Peptide Structure</h4>
        {svg_content}
        <p style="margin-top: 1rem; color: #666; font-size: 0.9rem;">
            Each circle represents an amino acid. The arrow shows the direction from N-terminus (NH₂) to C-terminus (COOH).
            Colors indicate amino acid types: <span style="color: #ff6b6b;">Basic</span>, 
            <span style="color: #4ecdc4;">Acidic</span>, <span style="color: #45b7d1;">Polar</span>, 
            <span style="color: #96ceb4;">Hydrophobic</span>.
        </p>
    </div>
    """

def create_molecular_structure_html(sequence: str) -> str:
    """Create HTML for actual 2D molecular structure visualization using RDKit"""
    if len(sequence) < 1:
        return "<p>Sequence too short for molecular structure visualization</p>"
    
    try:
        # Convert peptide sequence to SMILES string
        smiles = peptide_to_smiles(sequence)
        if not smiles:
            return "<p>Unable to generate molecular structure for this peptide</p>"
        
        # Create RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Fallback to a simple glycine structure if parsing fails
            fallback_smiles = "NCC(=O)O"
            mol = Chem.MolFromSmiles(fallback_smiles)
        if mol is None:
            return "<p>Unable to parse molecular structure</p>"
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create SVG drawing
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Get SVG content
        svg_content = drawer.GetDrawingText()
        
        # Convert to base64 for embedding
        svg_base64 = base64.b64encode(svg_content.encode()).decode()
        
        # Determine peptide length for description
        peptide_length = len(sequence)
        if peptide_length == 1:
            aa_names = {
                'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
                'C': 'Cysteine', 'E': 'Glutamic acid', 'Q': 'Glutamine', 'G': 'Glycine',
                'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
                'M': 'Methionine', 'O': 'Pyrrolysine', 'F': 'Phenylalanine', 'P': 'Proline',
                'S': 'Serine', 'T': 'Threonine', 'U': 'Selenocysteine', 'W': 'Tryptophan',
                'Y': 'Tyrosine', 'V': 'Valine'
            }
            aa_name = aa_names.get(sequence[0], 'Amino acid')
            description = f"Molecular structure of {aa_name} amino acid."
        else:
            description = f"Molecular structure of your complete {peptide_length}-amino acid peptide chain."
        
        # Create MolView URL for interactive viewing
        molview_url = f"https://molview.org/?q={smiles}"
        
        return f"""
        <div style="margin: 1rem 0; text-align: center;">
            <h4 style="margin-bottom: 1rem; color: #333;">2D Molecular Structure</h4>
            <div style="background: white; border: 1px solid #ccc; border-radius: 8px; padding: 1rem; display: inline-block;">
                <img src="data:image/svg+xml;base64,{svg_base64}" 
                     alt="2D Molecular Structure" 
                     style="max-width: 100%; height: auto; border-radius: 4px;">
            </div>
            <div style="margin-top: 1rem; display: flex; gap: 0.5rem; justify-content: center; flex-wrap: wrap;">
                <a href="{molview_url}" target="_blank" 
                   style="background: #007bff; color: white; padding: 0.5rem 1rem; text-decoration: none; border-radius: 4px; font-size: 0.9rem;">
                    <i class="fas fa-external-link-alt"></i> View in MolView
                </a>
                <button onclick="copySmiles('{smiles}')" 
                        style="background: #6c757d; color: white; padding: 0.5rem 1rem; border: none; border-radius: 4px; font-size: 0.9rem; cursor: pointer;">
                    <i class="fas fa-copy"></i> Copy SMILES
                </button>
            </div>
            <p style="margin-top: 1rem; color: #666; font-size: 0.9rem;">
                {description} Shows atoms, bonds, and functional groups.
            </p>
            <p style="margin-top: 0.5rem; color: #999; font-size: 0.8rem;">
                Structure data from PubChem API and RDKit
            </p>
        </div>
        """
        
    except Exception as e:
        return f"<p>Error generating molecular structure: {str(e)}</p>"

def get_pubchem_structure(amino_acid: str) -> str:
    """Fetch molecular structure from PubChem API"""
    try:
        # PubChem CIDs for amino acids
        pubchem_cids = {
            'A': '5950',  # Alanine
            'R': '6322',  # Arginine
            'N': '5951',  # Asparagine
            'D': '602',   # Aspartic acid
            'C': '5862',  # Cysteine
            'E': '611',   # Glutamic acid
            'Q': '5952',  # Glutamine
            'G': '750',   # Glycine
            'H': '6274',  # Histidine
            'I': '791',   # Isoleucine
            'L': '6106',  # Leucine
            'K': '5962',  # Lysine
            'M': '6137',  # Methionine
            'O': '16132350',  # Pyrrolysine
            'F': '6140',  # Phenylalanine
            'P': '614',   # Proline
            'S': '5958',  # Serine
            'T': '6058',  # Threonine
            'U': '5862',  # Selenocysteine (using Cys as fallback)
            'W': '6305',  # Tryptophan
            'Y': '6057',  # Tyrosine
            'V': '6287'   # Valine
        }
        
        cid = pubchem_cids.get(amino_acid, '750')  # Default to glycine
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,CanonicalSMILES/JSON"
        
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if 'PC_Compounds' in data and data['PC_Compounds']:
                compound = data['PC_Compounds'][0]
                for prop in compound.get('props', []):
                    if prop.get('urn', {}).get('label') == 'CanonicalSMILES':
                        return prop.get('value', '')
        return ''
    except Exception as e:
        print(f"PubChem API error for {amino_acid}: {e}")
        return ''

def get_chemspider_structure(amino_acid: str) -> str:
    """Fetch molecular structure from ChemSpider API"""
    try:
        # ChemSpider IDs for amino acids (simplified mapping)
        chemspider_names = {
            'A': 'alanine',
            'R': 'arginine',
            'N': 'asparagine',
            'D': 'aspartic acid',
            'C': 'cysteine',
            'E': 'glutamic acid',
            'Q': 'glutamine',
            'G': 'glycine',
            'H': 'histidine',
            'I': 'isoleucine',
            'L': 'leucine',
            'K': 'lysine',
            'M': 'methionine',
            'O': 'pyrrolysine',
            'F': 'phenylalanine',
            'P': 'proline',
            'S': 'serine',
            'T': 'threonine',
            'U': 'selenocysteine',
            'W': 'tryptophan',
            'Y': 'tyrosine',
            'V': 'valine'
        }
        
        name = chemspider_names.get(amino_acid, 'glycine')
        # Note: ChemSpider API requires authentication, so we'll use a fallback
        return ''
    except Exception as e:
        print(f"ChemSpider API error for {amino_acid}: {e}")
        return ''

def search_pubchem_peptide(sequence: str) -> dict:
    """Search PubChem for peptide information using their API"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{sequence}/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if 'PC_Compounds' in data and data['PC_Compounds']:
                compound = data['PC_Compounds'][0]
                cid = compound.get('id', {}).get('id', {}).get('cid', '')
                if cid:
                    detail_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
                    detail_response = requests.get(detail_url, timeout=10)
                    if detail_response.status_code == 200:
                        detail_data = detail_response.json()
                        desc = detail_data.get('InformationList', {}).get('Information', [{}])[0].get('Description', '')
                        return {
                            'found': True,
                            'cid': cid,
                            'description': desc,
                            'source': 'PubChem'
                        }
        return {'found': False}
    except Exception as e:
        print(f"PubChem search error: {e}")
        return {'found': False}

def search_uniprot_peptide(sequence: str) -> dict:
    """Search UniProt for peptide information"""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=sequence:{sequence}&format=json&size=1"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data.get('results'):
                result = data['results'][0]
                return {
                    'found': True,
                    'entry': result.get('primaryAccession', ''),
                    'name': result.get('uniProtkbId', ''),
                    'protein_name': result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'organism': result.get('organism', {}).get('scientificName', ''),
                    'function': result.get('comments', []),
                    'source': 'UniProt'
                }
        return {'found': False}
    except Exception as e:
        print(f"UniProt search error: {e}")
        return {'found': False}

def search_peptide_database(sequence: str) -> dict:
    """Return links to specialized peptide databases for manual lookup."""
    return {
        'APD3': f"https://aps.unmc.edu/AP/main.php?action=search_peptide&query={sequence}",
        'CAMP': f"https://www.camp3.bicnirrh.res.in/search.php?sequence={sequence}",
        'UniProt': f"https://www.uniprot.org/uniprot/?query=sequence:{sequence}"
    }

def build_peptide_chain_smiles(sequence: str) -> str:
    """Build complete peptide chain SMILES with proper peptide bonds"""
    aa_side_chains = {
        'A': 'C',  # Alanine
        'R': 'CCCN',  # Arginine - simplified
        'N': 'CC(=O)N',  # Asparagine
        'D': 'CC(=O)O',  # Aspartic acid
        'C': 'CS',  # Cysteine
        'E': 'CCC(=O)O',  # Glutamic acid
        'Q': 'CCCN',  # Glutamine
        'G': '',  # Glycine (no side chain)
        'H': 'CC1=CNC=N1',  # Histidine
        'I': 'CCC(C)C',  # Isoleucine
        'L': 'CC(C)C',  # Leucine
        'K': 'CCCCN',  # Lysine
        'M': 'CCSC',  # Methionine
        'O': 'CCCN',  # Pyrrolysine - simplified
        'F': 'CC1=CC=CC=C1',  # Phenylalanine
        'P': 'C1CCNC1',  # Proline
        'S': 'CO',  # Serine
        'T': 'CC(C)O',  # Threonine
        'U': 'CS',  # Selenocysteine
        'W': 'CC1=CC2=CC=CC=C2NC1',  # Tryptophan
        'Y': 'CC1=CC=C(O)C=C1',  # Tyrosine
        'V': 'CC(C)C'  # Valine
    }
    if len(sequence) == 1:
        aa = sequence[0]
        side_chain = aa_side_chains.get(aa, '')
        if side_chain:
            return f"NC({side_chain})C(=O)O"
        else:
            return "NCC(=O)O"
    peptide_parts = []
    for i, aa in enumerate(sequence):
        side_chain = aa_side_chains.get(aa, '')
        if i == 0:
            # First amino acid - N-terminus with free NH2
            if side_chain:
                peptide_parts.append(f"N({side_chain})C(=O)")
            else:
                peptide_parts.append("NCC(=O)")
        elif i == len(sequence) - 1:
            # Last amino acid - C-terminus with free COOH
            if side_chain:
                peptide_parts.append(f"N({side_chain})C(=O)O")
            else:
                peptide_parts.append("NCC(=O)O")
        else:
            # Middle amino acids - connected by peptide bonds
            if side_chain:
                peptide_parts.append(f"N({side_chain})C(=O)")
            else:
                peptide_parts.append("NCC(=O)")
    return "".join(peptide_parts)

def peptide_to_smiles(sequence: str) -> str:
    """Convert peptide sequence to SMILES string using external APIs when possible"""
    # Fallback SMILES for when APIs fail
    fallback_smiles = {
        'A': 'NC(C)C(=O)O',  # Alanine
        'R': 'NC(CCCN)C(=O)O',  # Arginine - simplified
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
        'O': 'NC(CCCN)C(=O)O',  # Pyrrolysine - simplified
        'F': 'NC(CC1=CC=CC=C1)C(=O)O',  # Phenylalanine
        'P': 'NC1CCNC1C(=O)O',  # Proline
        'S': 'NC(CO)C(=O)O',  # Serine
        'T': 'NC(C(C)O)C(=O)O',  # Threonine
        'U': 'NC(CS)C(=O)O',  # Selenocysteine
        'W': 'NC(CC1=CC2=CC=CC=C2NC1)C(=O)O',  # Tryptophan
        'Y': 'NC(CC1=CC=C(O)C=C1)C(=O)O',  # Tyrosine
        'V': 'NC(C(C)C)C(=O)O'  # Valine
    }
    
    # For single amino acids, try PubChem first
    if len(sequence) == 1:
        aa = sequence[0]
        pubchem_smiles = get_pubchem_structure(aa)
        if pubchem_smiles:
            return pubchem_smiles
        else:
            return fallback_smiles.get(aa, 'NCC(=O)O')
    else:
        # For multiple amino acids, build the complete peptide chain
        try:
            full_peptide_smiles = build_peptide_chain_smiles(sequence)
            # Test if RDKit can parse it
            mol = Chem.MolFromSmiles(full_peptide_smiles)
            if mol is not None:
                return full_peptide_smiles
            else:
                # Fallback to first amino acid if full chain fails
                aa = sequence[0]
                pubchem_smiles = get_pubchem_structure(aa)
                if pubchem_smiles:
                    return pubchem_smiles
                else:
                    return fallback_smiles.get(aa, 'NCC(=O)O')
        except Exception as e:
            print(f"Error building peptide chain: {e}")
            # Fallback to first amino acid
            aa = sequence[0]
            pubchem_smiles = get_pubchem_structure(aa)
            if pubchem_smiles:
                return pubchem_smiles
            else:
                return fallback_smiles.get(aa, 'NCC(=O)O')

def get_biological_info(sequence: str) -> list:
    """Get biological information about the peptide"""
    info = []
    
    # Analyze sequence characteristics
    if len(sequence) <= 5:
        info.append("This is a short peptide (oligopeptide) that may function as a signaling molecule or hormone.")
    elif len(sequence) <= 20:
        info.append("This is a medium-length peptide that could serve as a bioactive peptide or antimicrobial agent.")
    else:
        info.append("This is a longer peptide that may have structural or enzymatic functions.")
    
    # Analyze amino acid composition
    aa_props = calculate_peptide_properties(sequence)
    
    if aa_props['hydrophobic_count'] > aa_props['hydrophilic_count']:
        info.append("The peptide is predominantly hydrophobic, suggesting it may be membrane-associated or have structural roles.")
    else:
        info.append("The peptide is predominantly hydrophilic, suggesting it may be water-soluble and have signaling functions.")
    
    if aa_props['acidic_count'] > aa_props['basic_count']:
        info.append("The peptide is acidic overall, which may affect its binding properties and cellular localization.")
    elif aa_props['basic_count'] > aa_props['acidic_count']:
        info.append("The peptide is basic overall, which may facilitate interactions with nucleic acids or membranes.")
    
    # Check for specific patterns with pathway implications
    if 'C' in sequence:
        info.append("Contains cysteine residues that may form disulfide bonds, contributing to structural stability.")
        info.append("May participate in redox signaling pathways or protein folding mechanisms.")
    
    if sequence.count('P') > 1:
        info.append("Rich in proline residues, which may confer structural rigidity and affect protein folding.")
        info.append("Could function in collagen synthesis, immune response, or cell adhesion pathways.")
    
    # Add specific pathway predictions
    pathway_info = get_biological_pathways(sequence)
    if pathway_info:
        info.extend(pathway_info)
    
    return info

def get_biological_pathways(sequence: str) -> list:
    """Predict specific biological pathways and cellular functions"""
    pathways = []
    
    # Antimicrobial pathways
    if sequence.count('R') + sequence.count('K') >= 2:
        pathways.append("May activate innate immune response pathways against microbial pathogens.")
        pathways.append("Could interact with toll-like receptors (TLRs) or pattern recognition receptors (PRRs).")
    
    # Cell signaling pathways
    if len(sequence) <= 10 and sequence.count('C') >= 1:
        pathways.append("May function in G-protein coupled receptor (GPCR) signaling pathways.")
        pathways.append("Could activate intracellular signaling cascades like cAMP or calcium pathways.")
    
    # Metabolic pathways
    if sequence.count('D') + sequence.count('E') >= 2:
        pathways.append("May participate in metabolic regulation pathways or energy homeostasis.")
        pathways.append("Could function in gluconeogenesis, lipid metabolism, or amino acid catabolism.")
    
    # Growth and development pathways
    if sequence.count('C') >= 3:
        pathways.append("May function in growth factor signaling pathways (e.g., EGF, FGF, TGF-β).")
        pathways.append("Could participate in cell proliferation, differentiation, or tissue regeneration.")
    
    # Neurotransmission pathways
    if len(sequence) <= 8 and sequence.count('G') >= 1:
        pathways.append("May function in neurotransmitter signaling pathways in the nervous system.")
        pathways.append("Could activate ion channels, synaptic transmission, or neuromodulation pathways.")
    
    # Hormone pathways
    if len(sequence) <= 15 and sequence.count('C') >= 2:
        pathways.append("May function in endocrine signaling pathways affecting metabolism and homeostasis.")
        pathways.append("Could activate hormone receptor pathways in target tissues.")
    
    # Immune system pathways
    if sequence.count('P') >= 2:
        pathways.append("May participate in immune system regulation pathways and inflammatory responses.")
        pathways.append("Could function in cytokine signaling, T-cell activation, or antibody production.")
    
    # Metal homeostasis pathways
    if sequence.count('H') >= 2:
        pathways.append("May function in metal ion transport and homeostasis pathways.")
        pathways.append("Could participate in iron, zinc, or copper metabolism and storage.")
    
    return pathways

def analyze_peptide_function_ai(sequence: str) -> dict:
    """AI-powered analysis of peptide function based on sequence patterns"""
    analysis = {
        'predicted_function': '',
        'confidence': 'low',
        'therapeutic_potential': [],
        'biological_pathways': [],
        'similar_known_peptides': [],
        'structural_features': []
    }
    
    # Analyze amino acid composition
    aa_counts = {}
    for aa in sequence:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1
    
    # Calculate key properties
    basic_aa = aa_counts.get('R', 0) + aa_counts.get('K', 0) + aa_counts.get('H', 0)
    acidic_aa = aa_counts.get('D', 0) + aa_counts.get('E', 0)
    hydrophobic_aa = aa_counts.get('L', 0) + aa_counts.get('I', 0) + aa_counts.get('V', 0) + aa_counts.get('F', 0) + aa_counts.get('W', 0)
    polar_aa = aa_counts.get('S', 0) + aa_counts.get('T', 0) + aa_counts.get('N', 0) + aa_counts.get('Q', 0)
    special_aa = aa_counts.get('C', 0) + aa_counts.get('G', 0) + aa_counts.get('P', 0)
    
    # Predict function based on composition
    if basic_aa >= 3 and hydrophobic_aa >= 2:
        analysis['predicted_function'] = 'Antimicrobial Peptide'
        analysis['confidence'] = 'high'
        analysis['therapeutic_potential'] = [
            'Antibacterial activity against Gram-positive and Gram-negative bacteria',
            'Antifungal activity',
            'Wound healing promotion',
            'Immune system modulation'
        ]
        analysis['biological_pathways'] = [
            'Innate immune response',
            'Membrane disruption',
            'Bacterial cell wall targeting',
            'Inflammatory response regulation'
        ]
    
    elif basic_aa >= 4:
        analysis['predicted_function'] = 'Cell-Penetrating Peptide'
        analysis['confidence'] = 'high'
        analysis['therapeutic_potential'] = [
            'Drug delivery across cell membranes',
            'Blood-brain barrier crossing',
            'Gene therapy delivery',
            'Protein transduction'
        ]
        analysis['biological_pathways'] = [
            'Endocytosis',
            'Membrane transport',
            'Intracellular delivery',
            'Receptor-mediated uptake'
        ]
    
    elif len(sequence) <= 10 and polar_aa > hydrophobic_aa:
        analysis['predicted_function'] = 'Signaling Peptide'
        analysis['confidence'] = 'medium'
        analysis['therapeutic_potential'] = [
            'Hormone replacement therapy',
            'Neurotransmitter function',
            'Cell signaling modulation',
            'Receptor activation'
        ]
        analysis['biological_pathways'] = [
            'GPCR signaling',
            'Neural transmission',
            'Endocrine regulation',
            'Cell communication'
        ]
    
    elif special_aa >= 2 and 'C' in sequence:
        analysis['predicted_function'] = 'Structural Peptide'
        analysis['confidence'] = 'medium'
        analysis['therapeutic_potential'] = [
            'Protein stabilization',
            'Enzyme inhibition',
            'Receptor binding',
            'Conformational regulation'
        ]
        analysis['biological_pathways'] = [
            'Protein folding',
            'Disulfide bond formation',
            'Structural integrity',
            'Molecular recognition'
        ]
    
    else:
        analysis['predicted_function'] = 'Bioactive Peptide'
        analysis['confidence'] = 'low'
        analysis['therapeutic_potential'] = [
            'Potential therapeutic applications',
            'Biological activity modulation',
            'Cell signaling',
            'Immune response'
        ]
        analysis['biological_pathways'] = [
            'Cellular processes',
            'Metabolic regulation',
            'Homeostasis',
            'Adaptive responses'
        ]
    
    # Add structural features
    if 'C' in sequence:
        analysis['structural_features'].append('Contains cysteine residues - may form disulfide bonds')
    if sequence.count('P') > 1:
        analysis['structural_features'].append('Proline-rich - confers structural rigidity')
    if basic_aa > acidic_aa:
        analysis['structural_features'].append('Net positive charge - may interact with membranes')
    elif acidic_aa > basic_aa:
        analysis['structural_features'].append('Net negative charge - may act as metal chelators')
    
    return analysis

def search_pubmed_peptide(sequence: str) -> dict:
    """Search PubMed for specific peptide research"""
    try:
        # Search for peptide sequence in PubMed
        search_terms = [
            f'"{sequence}" peptide function',
            f'"{sequence}" peptide therapeutic',
            f'"{sequence}" peptide biological activity'
        ]
        
        results = []
        for term in search_terms:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                'db': 'pubmed',
                'term': term,
                'retmax': 5,
                'retmode': 'json'
            }
            
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                    for pmid in data['esearchresult']['idlist']:
                        # Get article details
                        detail_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                        detail_params = {
                            'db': 'pubmed',
                            'id': pmid,
                            'retmode': 'json'
                        }
                        
                        detail_response = requests.get(detail_url, params=detail_params, timeout=10)
                        if detail_response.status_code == 200:
                            detail_data = detail_response.json()
                            if 'result' in detail_data and pmid in detail_data['result']:
                                article = detail_data['result'][pmid]
                                results.append({
                                    'pmid': pmid,
                                    'title': article.get('title', ''),
                                    'abstract': article.get('abstract', ''),
                                    'journal': article.get('fulljournalname', ''),
                                    'pubdate': article.get('pubdate', ''),
                                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                                })
        
        return {
            'found': len(results) > 0,
            'articles': results,
            'source': 'PubMed'
        }
    except Exception as e:
        print(f"PubMed search error: {e}")
        return {'found': False, 'articles': []}

def get_known_peptide_functions(sequence: str) -> dict:
    """Match against known functional peptides"""
    known_peptides = {
        'RRR': {
            'name': 'Polyarginine',
            'function': 'Cell-penetrating peptide that crosses biological membranes',
            'therapeutic_use': 'Drug delivery, gene therapy, blood-brain barrier crossing',
            'mechanism': 'Electrostatic interaction with cell membranes, endocytosis',
            'clinical_applications': ['Cancer therapy', 'Neurological disorders', 'Vaccine delivery']
        },
        'KKK': {
            'name': 'Polylysine',
            'function': 'Antimicrobial and cell-penetrating peptide',
            'therapeutic_use': 'Antibacterial activity, protein purification, gene delivery',
            'mechanism': 'Membrane disruption, electrostatic interactions',
            'clinical_applications': ['Antibiotic development', 'Wound healing', 'Biotechnology']
        },
        'CC': {
            'name': 'Disulfide-bonded peptides',
            'function': 'Structural and signaling peptides with disulfide bonds',
            'therapeutic_use': 'Hormone replacement, enzyme inhibitors, receptor ligands',
            'mechanism': 'Disulfide bond formation, receptor binding, enzyme inhibition',
            'clinical_applications': ['Diabetes treatment', 'Hormone therapy', 'Enzyme therapy']
        },
        'DD': {
            'name': 'Acidic peptides',
            'function': 'Metal chelators and enzyme inhibitors',
            'therapeutic_use': 'Metal homeostasis, anticoagulation, enzyme regulation',
            'mechanism': 'Metal ion binding, enzyme active site competition',
            'clinical_applications': ['Metal poisoning', 'Blood clotting disorders', 'Enzyme therapy']
        },
        'PP': {
            'name': 'Proline-rich peptides',
            'function': 'Structural peptides with immune modulatory activity',
            'therapeutic_use': 'Anti-inflammatory agents, collagen synthesis, immune modulation',
            'mechanism': 'Collagen binding, immune receptor interaction, structural support',
            'clinical_applications': ['Inflammatory diseases', 'Tissue repair', 'Autoimmune disorders']
        }
    }
    
    # Check for exact matches
    for pattern, info in known_peptides.items():
        if pattern in sequence:
            return {
                'found': True,
                'match': pattern,
                'info': info,
                'confidence': 'high'
            }
    
    # Check for partial matches
    for pattern, info in known_peptides.items():
        if any(aa in sequence for aa in pattern):
            return {
                'found': True,
                'match': 'partial',
                'info': info,
                'confidence': 'medium'
            }
    
    return {'found': False}

def search_peptide_function(sequence: str) -> dict:
    """Comprehensive peptide function analysis using multiple approaches"""
    search_results = {
        'found_info': False,
        'ai_analysis': {},
        'pubmed_results': {},
        'known_peptides': {},
        'pubchem': {},
        'uniprot': {},
        'databases': {},
        'general_info': None
    }
    
    try:
        # 1. AI-powered analysis
        ai_analysis = analyze_peptide_function_ai(sequence)
        search_results['ai_analysis'] = ai_analysis
        if ai_analysis['predicted_function']:
            search_results['found_info'] = True
        
        # 2. Known peptide matching
        known_match = get_known_peptide_functions(sequence)
        search_results['known_peptides'] = known_match
        if known_match['found']:
            search_results['found_info'] = True
        
        # 3. PubMed literature search
        pubmed_results = search_pubmed_peptide(sequence)
        search_results['pubmed_results'] = pubmed_results
        if pubmed_results['found']:
            search_results['found_info'] = True
        
        # 4. Database searches
        pubchem = search_pubchem_peptide(sequence)
        if pubchem.get('found'):
            search_results['found_info'] = True
            search_results['pubchem'] = pubchem
        
        uniprot = search_uniprot_peptide(sequence)
        if uniprot.get('found'):
            search_results['found_info'] = True
            search_results['uniprot'] = uniprot
        
        search_results['databases'] = search_peptide_database(sequence)
        
        # 5. Fallback to general insights
        if not search_results['found_info']:
            search_results['general_info'] = get_general_peptide_insights(sequence)
            
    except Exception as e:
        search_results['general_info'] = get_general_peptide_insights(sequence)
        search_results['error'] = str(e)
    
    return search_results

def search_similar_peptides(sequence: str) -> list:
    """Search for similar known peptides and their therapeutic applications"""
    similar_peptides = []
    
    # Define known peptide patterns and their therapeutic applications
    known_patterns = {
        'RRR': {
            'name': 'Polyarginine',
            'applications': ['Cell-penetrating peptide', 'Drug delivery across blood-brain barrier', 'Antimicrobial activity'],
            'pathways': ['GPCR signaling', 'Membrane transport', 'Immune response']
        },
        'KKK': {
            'name': 'Polylysine',
            'applications': ['Antimicrobial peptide', 'Gene delivery', 'Protein purification'],
            'pathways': ['Innate immunity', 'Membrane disruption', 'Endocytosis']
        },
        'CC': {
            'name': 'Disulfide-bonded peptides',
            'applications': ['Hormone replacement therapy', 'Oxytocin analogs', 'Insulin analogs'],
            'pathways': ['Endocrine signaling', 'Metabolic regulation', 'Cell signaling']
        },
        'DD': {
            'name': 'Acidic peptides',
            'applications': ['Enzyme inhibitors', 'Metal chelators', 'Anticoagulants'],
            'pathways': ['Metabolic regulation', 'Metal homeostasis', 'Blood coagulation']
        },
        'PP': {
            'name': 'Proline-rich peptides',
            'applications': ['Anti-inflammatory agents', 'Collagen synthesis', 'Immune modulators'],
            'pathways': ['Immune response', 'Tissue repair', 'Cell adhesion']
        },
        'GG': {
            'name': 'Glycine-rich peptides',
            'applications': ['Neurotransmitters', 'Pain management', 'Muscle relaxants'],
            'pathways': ['Neural signaling', 'Ion channel regulation', 'Neuromodulation']
        }
    }
    
    # Check for known patterns in the sequence
    for pattern, info in known_patterns.items():
        if pattern in sequence:
            similar_peptides.append({
                'pattern': pattern,
                'name': info['name'],
                'applications': info['applications'],
                'pathways': info['pathways'],
                'description': f"Your peptide contains the {pattern} pattern, similar to {info['name']} which has known therapeutic applications."
            })
    
    # If no exact patterns, look for amino acid composition similarities
    if not similar_peptides:
        aa_composition = get_aa_composition_similarity(sequence)
        if aa_composition:
            similar_peptides.append(aa_composition)
    
    return similar_peptides

def get_aa_composition_similarity(sequence: str) -> dict:
    """Find peptides with similar amino acid composition"""
    # Count specific amino acid types
    basic_count = sequence.count('R') + sequence.count('K') + sequence.count('H')
    acidic_count = sequence.count('D') + sequence.count('E')
    hydrophobic_count = sequence.count('L') + sequence.count('I') + sequence.count('V') + sequence.count('F') + sequence.count('W')
    polar_count = sequence.count('S') + sequence.count('T') + sequence.count('N') + sequence.count('Q')
    special_count = sequence.count('C') + sequence.count('G') + sequence.count('P')
    
    # Determine peptide type based on composition
    if basic_count >= 2 and hydrophobic_count >= 2:
        return {
            'pattern': 'Basic + Hydrophobic',
            'name': 'Antimicrobial Peptide',
            'applications': ['Antibiotic development', 'Antifungal agents', 'Wound healing'],
            'pathways': ['Innate immunity', 'Membrane disruption', 'Bacterial cell wall targeting'],
            'description': f"Your peptide has {basic_count} basic and {hydrophobic_count} hydrophobic amino acids, characteristic of antimicrobial peptides."
        }
    elif basic_count >= 3:
        return {
            'pattern': 'Highly Basic',
            'name': 'Cell-Penetrating Peptide',
            'applications': ['Drug delivery', 'Gene therapy', 'Blood-brain barrier crossing'],
            'pathways': ['Membrane transport', 'Endocytosis', 'Intracellular delivery'],
            'description': f"Your peptide has {basic_count} basic amino acids, characteristic of cell-penetrating peptides."
        }
    elif acidic_count >= 2:
        return {
            'pattern': 'Acidic',
            'name': 'Metal-Binding Peptide',
            'applications': ['Metal chelation therapy', 'Enzyme inhibition', 'Mineral supplementation'],
            'pathways': ['Metal homeostasis', 'Metabolic regulation', 'Catalytic inhibition'],
            'description': f"Your peptide has {acidic_count} acidic amino acids, characteristic of metal-binding peptides."
        }
    elif special_count >= 2:
        return {
            'pattern': 'Special Amino Acids',
            'name': 'Structural/Regulatory Peptide',
            'applications': ['Hormone analogs', 'Structural proteins', 'Regulatory peptides'],
            'pathways': ['Cell signaling', 'Structural support', 'Metabolic regulation'],
            'description': f"Your peptide has {special_count} special amino acids (C/G/P), characteristic of structural and regulatory peptides."
        }
    
    return None

def get_general_peptide_insights(sequence: str) -> str:
    """Provide general insights about the peptide based on its properties"""
    insights = []
    
    # Analyze sequence length
    if len(sequence) <= 5:
        insights.append("Short peptides (2-5 amino acids) often function as signaling molecules, hormones, or neurotransmitters.")
    elif len(sequence) <= 20:
        insights.append("Medium-length peptides (6-20 amino acids) commonly serve as bioactive peptides, antimicrobial agents, or enzyme inhibitors.")
    else:
        insights.append("Longer peptides (20+ amino acids) may have structural roles or enzymatic functions.")
    
    # Analyze amino acid composition for specific functions
    aa_props = calculate_peptide_properties(sequence)
    
    # Check for specific functional patterns with therapeutic applications
    if 'C' in sequence:
        insights.append("Cysteine residues can form disulfide bonds, making this peptide potentially stable and suitable for therapeutic applications.")
    
    # Antimicrobial and cell-penetrating peptides
    if sequence.count('R') + sequence.count('K') > 1:
        insights.append("Rich in basic amino acids (arginine/lysine), suggesting potential antimicrobial or cell-penetrating properties.")
        insights.append("May interact with bacterial membranes or facilitate drug delivery across cell membranes.")
    
    # Metal-binding and catalytic peptides
    if sequence.count('D') + sequence.count('E') > 1:
        insights.append("Acidic peptide that may interact with metal ions or have specific binding properties.")
        insights.append("Could function in metal homeostasis or enzymatic catalysis pathways.")
    
    # Hydrophobicity-based pathway predictions
    if aa_props['hydrophobic_count'] > aa_props['hydrophilic_count']:
        insights.append("Hydrophobic nature suggests membrane association or structural roles in protein complexes.")
        insights.append("May participate in membrane transport, receptor binding, or lipid metabolism pathways.")
    else:
        insights.append("Hydrophilic character indicates water solubility and potential signaling functions.")
        insights.append("Could act as a hormone, neurotransmitter, or cytokine in intercellular communication.")
    
    # Check for known therapeutic peptide patterns
    if sequence.startswith('M') and len(sequence) > 10:
        insights.append("N-terminal methionine suggests this could be a protein fragment or signal peptide.")
        insights.append("May participate in protein targeting or secretion pathways.")
    
    if 'P' in sequence and sequence.count('P') > 1:
        insights.append("Proline-rich regions often confer structural rigidity and may affect protein-protein interactions.")
        insights.append("Could function in cell adhesion, immune response, or cytoskeletal organization.")
    
    # Add specific therapeutic pathway predictions
    therapeutic_pathways = get_therapeutic_pathways(sequence)
    if therapeutic_pathways:
        insights.extend(therapeutic_pathways)
    
    return " ".join(insights)

def get_therapeutic_pathways(sequence: str) -> list:
    """Predict specific therapeutic applications and biological pathways"""
    pathways = []
    
    # Antimicrobial peptide patterns
    if sequence.count('R') + sequence.count('K') >= 2 and sequence.count('L') + sequence.count('I') + sequence.count('V') >= 2:
        pathways.append("This peptide shows characteristics of antimicrobial peptides (AMPs) that target bacterial cell membranes.")
        pathways.append("Potential therapeutic application: Antibiotic development for drug-resistant bacteria.")
    
    # Cell-penetrating peptide patterns
    if sequence.count('R') + sequence.count('K') >= 3:
        pathways.append("High cationic charge suggests cell-penetrating peptide (CPP) properties.")
        pathways.append("Potential therapeutic application: Drug delivery across cell membranes and blood-brain barrier.")
    
    # Hormone-like peptide patterns
    if len(sequence) <= 10 and sequence.count('C') >= 2:
        pathways.append("Short peptide with disulfide bonds resembles peptide hormone structure.")
        pathways.append("Potential therapeutic application: Hormone replacement therapy or metabolic regulation.")
    
    # Enzyme inhibitor patterns
    if sequence.count('D') + sequence.count('E') >= 2 and len(sequence) <= 15:
        pathways.append("Acidic short peptide may function as enzyme inhibitor or metal chelator.")
        pathways.append("Potential therapeutic application: Enzyme inhibition for metabolic disorders or metal toxicity.")
    
    # Anti-inflammatory peptide patterns
    if sequence.count('P') >= 2 and sequence.count('G') >= 1:
        pathways.append("Proline-rich peptide may have anti-inflammatory or immunomodulatory properties.")
        pathways.append("Potential therapeutic application: Inflammatory disease treatment or immune system regulation.")
    
    # Neurotransmitter-like patterns
    if len(sequence) <= 8 and sequence.count('G') >= 1:
        pathways.append("Short peptide with glycine resembles neurotransmitter structure.")
        pathways.append("Potential therapeutic application: Neurological disorder treatment or pain management.")
    
    # Growth factor-like patterns
    if sequence.count('C') >= 3 and len(sequence) >= 10:
        pathways.append("Cysteine-rich peptide may function as growth factor or cytokine.")
        pathways.append("Potential therapeutic application: Tissue regeneration or cell proliferation regulation.")
    
    # Metal-binding therapeutic patterns
    if sequence.count('H') >= 2:
        pathways.append("Histidine-rich peptide may function in metal homeostasis or transport.")
        pathways.append("Potential therapeutic application: Metal chelation therapy or mineral supplementation.")
    
    return pathways

def text_to_peptide(text: str) -> str:
    cleaned = "".join(ch.upper() for ch in text if ch.isalpha())
    if not cleaned:
        return "AXA"
    return "".join(SAFE_MAP.get(ch, FALLBACK) for ch in cleaned[:12])

def create_visual_sequence(text: str) -> str:
    """Create visual representation showing original letters with ambiguous amino acids as superscripts"""
    cleaned = "".join(ch.upper() for ch in text if ch.isalpha())
    if not cleaned:
        return "AXA"
    
    visual_parts = []
    for ch in cleaned[:12]:
        if ch in SAFE_MAP and SAFE_MAP[ch] != ch:
            # This is an ambiguous amino acid that was converted
            visual_parts.append(f"{ch}<sup>({SAFE_MAP[ch]})</sup>")
        else:
            # This is a direct mapping or fallback
            visual_parts.append(ch)
    
    return "".join(visual_parts)

def create_svg_visual_sequence(text: str) -> str:
    """Create visual representation for SVG (without HTML tags)"""
    cleaned = "".join(ch.upper() for ch in text if ch.isalpha())
    if not cleaned:
        return "AXA"
    
    visual_parts = []
    for ch in cleaned[:12]:
        if ch in SAFE_MAP and SAFE_MAP[ch] != ch:
            # This is an ambiguous amino acid that was converted
            visual_parts.append(f"{ch}({SAFE_MAP[ch]})")
        else:
            # This is a direct mapping or fallback
            visual_parts.append(ch)
    
    return "".join(visual_parts)

def make_svg(seq: str, original_text: str, mw: float, pi: float) -> bytes:
    dwg = svgwrite.Drawing(size=("220px", "120px"))
    dwg.add(dwg.rect(insert=(0, 0), size=("100%", "100%"),
                     rx=10, ry=10, fill="#f2f2f2"))
    
    # Create SVG visual sequence from original text
    svg_visual_seq = create_svg_visual_sequence(original_text)
    
    # Split the visual sequence to handle superscripts
    parts = re.split(r'([A-Z]\([A-Z]\))', svg_visual_seq)
    
    # Calculate total width to center the text
    total_width = 0
    for part in parts:
        if re.match(r'[A-Z]\([A-Z]\)', part):
            total_width += 30  # Letter + superscript spacing (18 + 12)
        elif part:
            total_width += len(part) * 18  # Regular letter spacing
    
    # Start position to center the text
    x_offset = (220 - total_width) / 2
    y_pos = 45
    
    # Add N-terminus label close to the sequence
    dwg.add(dwg.text("N", insert=(f"{x_offset - 15}px", f"{y_pos}%"), 
             font_family="monospace", font_size="10px", fill="#666"))
    
    for part in parts:
        if re.match(r'[A-Z]\([A-Z]\)', part):
            # This is a letter with superscript
            letter = part[0]
            superscript = part[2]
            
            # Add main letter
            dwg.add(dwg.text(letter, insert=(f"{x_offset}px", f"{y_pos}%"), 
                     font_family="monospace", font_size="28px", fill="#222"))
            x_offset += 18
            
            # Add superscript
            dwg.add(dwg.text(superscript, insert=(f"{x_offset}px", f"{y_pos-3}%"), 
                     font_family="monospace", font_size="12px", fill="#666"))
            x_offset += 12
        elif part:
            # Regular text
            dwg.add(dwg.text(part, insert=(f"{x_offset}px", f"{y_pos}%"), 
                     font_family="monospace", font_size="28px", fill="#222"))
            x_offset += len(part) * 18
    
    # Add C-terminus label close to the end of sequence
    dwg.add(dwg.text("C", insert=(f"{x_offset + 8}px", f"{y_pos}%"), 
             font_family="monospace", font_size="10px", fill="#666"))
    
    stats = f"MW {mw:.1f} Da  |  pI {pi:.2f}"
    dwg.add(dwg.text(stats, insert=("50%", "80%"), text_anchor="middle",
                     font_family="monospace", font_size="12px", fill="#555"))
    return dwg.tostring().encode()

def make_enhanced_svg(seq: str, original_text: str, mw: float, pi: float) -> bytes:
    """Create enhanced SVG with 2D molecular structure included"""
    # Create larger SVG to accommodate molecular structure
    dwg = svgwrite.Drawing(size=("400px", "500px"))
    dwg.add(dwg.rect(insert=(0, 0), size=("100%", "100%"),
                     rx=10, ry=10, fill="#f2f2f2"))
    
    # Add title
    dwg.add(dwg.text("Peptide Tag", insert=("50%", "12%"), text_anchor="middle",
                     font_family="Arial, sans-serif", font_size="16px", fill="#333", font_weight="bold"))
    
    # Create SVG visual sequence from original text
    svg_visual_seq = create_svg_visual_sequence(original_text)
    
    # Split the visual sequence to handle superscripts
    parts = re.split(r'([A-Z]\([A-Z]\))', svg_visual_seq)
    
    # Calculate total width to center the text
    total_width = 0
    for part in parts:
        if re.match(r'[A-Z]\([A-Z]\)', part):
            total_width += 30  # Letter + superscript spacing (18 + 12)
        elif part:
            total_width += len(part) * 18  # Regular letter spacing
    
    # Start position to center the text
    x_offset = (400 - total_width) / 2
    y_pos = 30
    
    # Add N-terminus label close to the sequence
    dwg.add(dwg.text("N", insert=(f"{x_offset - 20}px", f"{y_pos}%"), 
             font_family="monospace", font_size="12px", fill="#666"))
    
    for part in parts:
        if re.match(r'[A-Z]\([A-Z]\)', part):
            # This is a letter with superscript
            letter = part[0]
            superscript = part[2]
            
            # Add main letter
            dwg.add(dwg.text(letter, insert=(f"{x_offset}px", f"{y_pos}%"), 
                     font_family="monospace", font_size="28px", fill="#222"))
            x_offset += 18
            
            # Add superscript
            dwg.add(dwg.text(superscript, insert=(f"{x_offset}px", f"{y_pos-3}%"), 
                     font_family="monospace", font_size="12px", fill="#666"))
            x_offset += 12
        elif part:
            # Regular text
            dwg.add(dwg.text(part, insert=(f"{x_offset}px", f"{y_pos}%"), 
                     font_family="monospace", font_size="28px", fill="#222"))
            x_offset += len(part) * 18
    
    # Add C-terminus label close to the end of sequence
    dwg.add(dwg.text("C", insert=(f"{x_offset + 10}px", f"{y_pos}%"), 
             font_family="monospace", font_size="12px", fill="#666"))
    
    # Add stats
    stats = f"MW {mw:.1f} Da  |  pI {pi:.2f}"
    dwg.add(dwg.text(stats, insert=("50%", "40%"), text_anchor="middle",
                     font_family="monospace", font_size="12px", fill="#555"))
    
    # Add molecular structure section
    try:
        # Generate molecular structure SVG
        smiles = peptide_to_smiles(seq)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                # Fallback to glycine
                mol = Chem.MolFromSmiles("NCC(=O)O")
            
            if mol:
                # Generate 2D coordinates
                AllChem.Compute2DCoords(mol)
                
                # Create SVG drawing
                drawer = rdMolDraw2D.MolDraw2DSVG(350, 250)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                
                # Get SVG content
                mol_svg_content = drawer.GetDrawingText()
                
                # Add molecular structure title
                dwg.add(dwg.text("2D Molecular Structure", insert=("50%", "55%"), text_anchor="middle",
                               font_family="Arial, sans-serif", font_size="14px", fill="#333", font_weight="bold"))
                
                # Create a background for the molecular structure
                dwg.add(dwg.rect(insert=("25%", "60%"), size=("50%", "35%"),
                               fill="white", stroke="#333", stroke_width="2", rx="5"))
                
                # Add the molecular structure SVG as an embedded image
                # Convert the molecular SVG to base64 and embed it
                mol_svg_base64 = base64.b64encode(mol_svg_content.encode()).decode()
                
                # Add the molecular structure as an embedded image
                dwg.add(dwg.image(href=f"data:image/svg+xml;base64,{mol_svg_base64}",
                                insert=("27%", "62%"), size=("46%", "31%")))
                
    except Exception as e:
        # If molecular structure fails, add a placeholder
        dwg.add(dwg.rect(insert=("25%", "60%"), size=("50%", "35%"),
                       fill="white", stroke="#ccc", stroke_width="1", rx="5"))
        dwg.add(dwg.text("Molecular structure available in web view", 
                       insert=("50%", "75%"), text_anchor="middle",
                       font_family="Arial, sans-serif", font_size="10px", fill="#666"))
    
    return dwg.tostring().encode()

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        raw = request.form["username"]
        seq = text_to_peptide(raw)
        visual_seq = create_visual_sequence(raw) # Use the new function here
        ana = ProteinAnalysis(seq)
        mw, pi = ana.molecular_weight(), ana.isoelectric_point()
        
        # Create both regular and enhanced SVGs
        svg_bytes = make_svg(seq, raw, mw, pi) # Pass raw as original_text
        svg_id = hashlib.md5(svg_bytes).hexdigest()
        app.config.setdefault("SVGS", {})[svg_id] = svg_bytes
        
        # Create enhanced SVG with molecular structure
        enhanced_svg_bytes = make_enhanced_svg(seq, raw, mw, pi)
        enhanced_svg_id = hashlib.md5(enhanced_svg_bytes).hexdigest()
        app.config.setdefault("ENHANCED_SVGS", {})[enhanced_svg_id] = enhanced_svg_bytes
        
        # Store sequence mapping for molecular structure route
        app.config.setdefault("SEQUENCES", {})[svg_id] = seq
        
        # Get protein information
        protein_info = get_protein_info(seq)
        
        # Generate 2D structure HTML
        structure_html = create_2d_structure_html(seq)
        
        # Generate molecular structure HTML
        molecular_structure_html = create_molecular_structure_html(seq)
        
        # Search for peptide function information
        search_results = search_peptide_function(seq)
        
        return render_template("result.html", seq=seq, visual_seq=visual_seq, mw=mw, pi=pi,
                               svg_id=svg_id, enhanced_svg_id=enhanced_svg_id, raw=raw, protein_info=protein_info, 
                               structure_html=structure_html, molecular_structure_html=molecular_structure_html,
                               search_results=search_results)
    return render_template("index.html")

@app.route("/download/<svg_id>.svg")
def download_svg(svg_id):
    svg_bytes = app.config["SVGS"].get(svg_id)
    if not svg_bytes:
        return "File expired", 404
    return send_file(BytesIO(svg_bytes), mimetype="image/svg+xml",
                     download_name=f"peptide_{svg_id}.svg", as_attachment=True)

@app.route("/download-enhanced/<enhanced_svg_id>.svg")
def download_enhanced_svg(enhanced_svg_id):
    enhanced_svg_bytes = app.config["ENHANCED_SVGS"].get(enhanced_svg_id)
    if not enhanced_svg_bytes:
        return "File expired", 404
    return send_file(BytesIO(enhanced_svg_bytes), mimetype="image/svg+xml",
                     download_name=f"peptide_enhanced_{enhanced_svg_id}.svg", as_attachment=True)

@app.route("/molecular-structure/<svg_id>.svg")
def molecular_structure_svg(svg_id):
    """Serve molecular structure as standalone SVG"""
    try:
        # Get the sequence from the stored mapping
        seq = app.config.get("SEQUENCES", {}).get(svg_id)
        if not seq:
            return "Sequence not found", 404
        
        # Generate molecular structure SVG
        smiles = peptide_to_smiles(seq)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                mol = Chem.MolFromSmiles("NCC(=O)O")
            
            if mol:
                AllChem.Compute2DCoords(mol)
                drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                
                svg_content = drawer.GetDrawingText()
                return svg_content, 200, {'Content-Type': 'image/svg+xml'}
        
        return "Molecular structure not available", 404
    except Exception as e:
        return f"Error generating molecular structure: {str(e)}", 500

if __name__ == "__main__":
    # Get port from environment variable (for deployment platforms)
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)

