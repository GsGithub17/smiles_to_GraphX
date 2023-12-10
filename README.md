# Molecular SMILES to Graph

## Overview
This project involves analyzing molecular compounds using their SMILES (Simplified Molecular Input Line Entry System) representations. The analysis consists of three major steps:
1. Plotting the 2D structures of the compounds with IDs.
   ![plot](first_10_compounds.png?raw=true)
2. Computing and plotting a Tanimoto similarity matrix.
   ![plot](similarity_matrix.png?raw=true)
4. Creating and plotting a graph based on the similarity matrix with nodes color-coded by pIC50 values.
   ![plot](similarity_graph.png?raw=true)

## Requirements
- Python
- RDKit
- Pandas
- NumPy
- Seaborn
- Matplotlib
- NetworkX

Install the necessary libraries using the following command:

```bash
pip install pandas numpy seaborn matplotlib networkx rdkit-pypi
```
## Usage

## Step 1: Plotting 2D Structures
This step involves generating and saving images of the first 10 compounds from the input dataset. Each image includes the compound's 2D structure along with its corresponding ID.

## Step 2: Computing and Plotting Similarity Matrix
The Tanimoto similarity matrix is computed using the ECFP (Extended-Connectivity Fingerprints) fingerprints derived from the SMILES strings. A heatmap of this similarity matrix is then plotted and saved.

## Step 3: Creating and Plotting the Graph
A graph is created based on the similarity matrix, where each node represents a compound, and edges represent similarity between compounds. Nodes are color-coded according to their pIC50 values, and the graph is saved as an image.

## Example

To run the analysis, ensure your data file is in CSV format with columns 'ID', 'SMILES', and 'pIC50'. Replace 'your_data.csv' with the path to your data file.

```python
df = pd.read_csv('your_data.csv')
analyze_compounds(df)
```

## Output

The output includes:

- An image of the first 10 compounds' 2D structures.
- A heatmap of the Tanimoto similarity matrix.
- A graph image with nodes color-coded by pIC50 values.

## Tips
To Change the layout ot the Graph, refere to the following site: https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.kamada_kawai_layout.html

## License

This project is licensed under the MIT License - see the LICENSE file for details.
