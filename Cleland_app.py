import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import requests
from itertools import combinations
from scipy.special import comb
from io import BytesIO

def fetch_protein_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    try:
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)
        response.raise_for_status()
        fasta = response.text
        sequence = ''.join(fasta.split('\n')[1:]).replace(' ', '')
        return sequence
    except requests.RequestException as e:
        st.error(f"Failed to fetch data for UniProt ID {uniprot_id}: {e}")
        return None

def calculate_molecular_mass(sequence):
    """Calculate molecular mass of the protein based on the sequence."""
    # Approximate average mass of an amino acid residue in kDa
    avg_residue_mass = 0.110  # kDa (110 Da)
    molecular_mass = len(sequence) * avg_residue_mass
    return molecular_mass

def generate_proteoforms(num_cysteines):
    """Generate all unique proteoforms and group them by oxidation state."""
    proteoforms = []
    grouped_proteoforms = [[] for _ in range(num_cysteines + 1)]
    
    for i in range(num_cysteines + 1):
        for combi in combinations(range(num_cysteines), i):
            proteoform = np.zeros(num_cysteines, dtype=int)
            proteoform[list(combi)] = 1
            proteoforms.append(proteoform)
            grouped_proteoforms[i].append(proteoform)
    
    return proteoforms, grouped_proteoforms

def plot_immunoblot(molecular_mass, grouped_proteoforms, num_cysteines):
    """Plot the positions of the redox proteoforms on a simulated immunoblot."""
    redox_grades = [(100 * (num_cysteines - k) / num_cysteines) for k in range(num_cysteines + 1)]
    
    # Calculate band positions on the blot
    band_positions = [molecular_mass + (5 * i) for i in range(len(redox_grades))]
    
    # Simulate the appearance on a 4-15% gradient gel
    fig, ax = plt.subplots(figsize=(5, 8))
    
    for i, pos in enumerate(band_positions):
        band_intensity = (num_cysteines - i + 1) / (num_cysteines + 1)  # Intensity decreases with oxidation
        ax.plot([0.3, 0.7], [pos, pos], linewidth=10 * band_intensity, color='black')
        ax.text(0.75, pos, f'{redox_grades[i]:.1f}%', verticalalignment='center', fontsize=12)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, band_positions[-1] + 10)
    ax.set_yticks(np.arange(0, band_positions[-1] + 10, 5))
    ax.set_yticklabels([f'{tick:.0f} kDa' for tick in ax.get_yticks()])
    ax.set_xticks([])
    ax.set_xlabel('Protein Redox States', fontsize=15)
    ax.set_ylabel('Molecular Mass (kDa)', fontsize=15)
    ax.set_title('Simulated Immunoblot', fontsize=20)
    
    # Invert the y-axis to match the appearance of a real blot
    ax.invert_yaxis()
    
    plt.tight_layout()
    
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    return buf

# Streamlit app
st.title('Cysteine Redox Proteoforms Immunoblot Simulation')

uniprot_id = st.text_input("Enter UniProt Accession Number:", "P04406")

if uniprot_id:
    sequence = fetch_protein_sequence(uniprot_id)
    if sequence:
        num_cysteines = sequence.count('C')  # Count the number of cysteines
        molecular_mass = calculate_molecular_mass(sequence)

        st.write(f"Protein Sequence Length: {len(sequence)} amino acids")
        st.write(f"Molecular Mass: {molecular_mass:.2f} kDa")
        st.write(f"Number of Cysteines: {num_cysteines}")

        if num_cysteines > 0:
            _, grouped_proteoforms = generate_proteoforms(num_cysteines)
            buf = plot_immunoblot(molecular_mass, grouped_proteoforms, num_cysteines)

            if buf:
                # Display the plot
                st.image(buf, use_column_width=True, caption='Simulated Immunoblot')

                # Download button
                buf.seek(0)
                st.download_button(
                    label="Download Immunoblot",
                    data=buf,
                    file_name="Simulated_Immunoblot.png",
                    mime="image/png"
                )
