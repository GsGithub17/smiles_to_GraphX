import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.colors as mcolors
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd

# Function to plot 2D structures
def plot_2D_structures(df):
    mols = [Chem.MolFromSmiles(smiles) for smiles in df['SMILES'][:10]]
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200), legends=df['ID'][:10].tolist(), returnPNG=False)
    img.save('first_10_compounds.png')
    #img.show()

# Function to compute fingerprints from SMILES
def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(2048)  # Return a zero vector if the molecule can't be parsed
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

# Function to compute and plot Tanimoto similarity matrix
def compute_and_plot_similarity(df):
    fingerprints = df['SMILES'].apply(smiles_to_fp)
    num_molecules = len(fingerprints)
    similarity_matrix = np.zeros((num_molecules, num_molecules))
    for i in range(num_molecules):
        for j in range(i+1, num_molecules):
            similarity = DataStructs.FingerprintSimilarity(fingerprints[i], fingerprints[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    mask = np.triu(np.ones_like(similarity_matrix, dtype=bool))
    sns.heatmap(similarity_matrix, mask=mask, cmap='viridis', square=True)
    plt.title('Tanimoto Similarity Half Triangle')
    plt.xlabel('Molecule Index')
    plt.ylabel('Molecule Index')
    plt.savefig('similarity_matrix.png')
    plt.show()

def create_and_plot_graph(df):
    # Create a dictionary mapping index to fingerprint
    index_to_fp = {i: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(row['SMILES']), 2, nBits=2048) 
                   for i, row in df.iterrows()}
    
    # Create a dictionary mapping index to pIC50 score
    index_to_score = df['pIC50'].to_dict()
    
    # Create the graph
    similarity_graph = nx.Graph()
    for idx1, fp1 in index_to_fp.items():
        for idx2, fp2 in index_to_fp.items():
            if idx1 != idx2:
                similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
                if similarity > 0.95:
                    similarity_graph.add_edge(idx1, idx2, weight=similarity)

    # Set up colors for nodes based on pIC50 score
    plt.figure(figsize=(15, 15))
    min_score = min(index_to_score.values())
    max_score = max(index_to_score.values())
    norm = mcolors.Normalize(vmin=min_score, vmax=max_score, clip=True)
    colors = [plt.cm.viridis(norm(index_to_score[node])) for node in similarity_graph.nodes()]

    # Draw the graph using kamada_kawai_layout
    fig, ax = plt.subplots()
    pos = nx.kamada_kawai_layout(similarity_graph)
    nx.draw(similarity_graph, pos, with_labels=True, node_color=colors, ax=ax)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('pIC50')
    plt.savefig('similarity_graph.png')


    
# Main function to analyze compounds
def analyze_compounds(df):
    plot_2D_structures(df)
    compute_and_plot_similarity(df)
    create_and_plot_graph(df)
    # Call your graph plotting function here (Part 3)

# Example usage
df = pd.read_csv('path_to_your_input_file.csv')
analyze_compounds(df)
