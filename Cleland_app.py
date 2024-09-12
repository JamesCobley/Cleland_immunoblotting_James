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

def plot_immunoblot(molecular_mass, grouped_proteoforms, num_cysteines, coefficients):
    """Plot the positions of the redox proteoforms on a scale-invariant simulated immunoblot."""
    
    # Fixed molecular weight markers from 10 kDa to 250 kDa
    marker_positions = np.array([10, 25, 37, 50, 75, 100, 150, 250])
    
    # Calculate band positions using molecular weights
    band_positions = [molecular_mass + (5 * i) for i in range(len(grouped_proteoforms))]
    
    fig, ax = plt.subplots(figsize=(5, 8))

    # Plot the molecular weight markers (10 to 250 kDa)
    for marker in marker_positions:
        y_pos = predict_band_position(marker, coefficients)  # Use the standard curve to place marker
        ax.plot([0.2, 0.8], [y_pos, y_pos], 'r--', lw=2)  # Plot red dashed line for marker
        ax.text(0.85, y_pos, f'{marker} kDa', verticalalignment='center', fontsize=12, color='red')

    # Plot the redox proteoforms at their corresponding molecular weights
    for i, pos in enumerate(band_positions):
        y_pos = predict_band_position(pos, coefficients)  # Position the band using the scaling
        band_intensity = (num_cysteines - i + 1) / (num_cysteines + 1)  # Intensity decreases with oxidation
        ax.plot([0.3, 0.7], [y_pos, y_pos], linewidth=10 * band_intensity, color='black')
        ax.text(0.75, y_pos, f'{(100 * (num_cysteines - i) / num_cysteines):.1f}%', verticalalignment='center', fontsize=12)

    # Set fixed y-axis limits (scale invariant)
    ax.set_ylim(predict_band_position(250, coefficients), predict_band_position(10, coefficients))
    
    ax.set_xlim(0, 1)
    ax.set_yticks([predict_band_position(mw, coefficients) for mw in marker_positions])
    ax.set_yticklabels([f'{mw:.0f} kDa' for mw in marker_positions])
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
        
        # Calculate the molecular mass of the 100%-oxidised form
        oxidised_mass = molecular_mass + (num_cysteines * 5)

        st.write(f"Protein Sequence Length: {len(sequence)} amino acids")
        st.write(f"Molecular Mass (Reduced Form): {molecular_mass:.2f} kDa")
        st.write(f"Molecular Mass (100%-Oxidised Form): {oxidised_mass:.2f} kDa")
        st.write(f"Number of Cysteines: {num_cysteines}")

        # Cleland Immunoblot suitability
        if oxidised_mass < 152:
            st.success("Yes, this protein is a good candidate for Cleland immunoblotting.")
        else:
            st.warning("No, this protein is not a good candidate for Cleland immunoblotting.")

        if num_cysteines > 0:
            # Coefficients from the standard curve
            coefficients = [-0.00501309, 2.38407094]  # Replace with real coefficients
            _, grouped_proteoforms = generate_proteoforms(num_cysteines)
            buf = plot_immunoblot(molecular_mass, grouped_proteoforms, num_cysteines, coefficients)

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
