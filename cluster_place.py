from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import umap.umap_ as umap

def arguments():
    parser = argparse.ArgumentParser(description="Finds Pantoea plasmids based on homology")
    parser.add_argument('-i', '--infile', required=True, type=str, help="Fasta with contigs to be classified fasta")
    args = parser.parse_args()
    return args

def main():
    args=arguments()
    with open(args.infile, 'r') as file:
        lines = file.readlines()
    lines=lines[1:]
    dimension_sq=len(lines[1:])
    arr = np.zeros((dimension_sq, dimension_sq), dtype=float)

    line_num=0

    for line in lines:
        line_no_enter=line.replace('\n','')
        line_arr=line.split("\t")[1:]
        entry_num=0

        for entry in line_arr:
            try:
                arr[line_num, entry_num]=entry
                try:
                    arr[entry_num, line_num]=entry
                except:
                    print(line_num, entry_num)
                entry_num+=1
            except:
                print(line_num, entry_num)
        line_num+=1
    print(arr)

    condensed_distance = squareform(arr)

    # Calculate the linkage matrix using single linkage
    linkage_matrix = linkage(condensed_distance, method='single')

    # Replace 'threshold_value' with your desired threshold for clustering
    threshold_value = 0.01

    # Get the clustered elements using fcluster
    clusters = fcluster(linkage_matrix, threshold_value, criterion='distance')

    print("Clusters:", clusters)

    reducer = umap.UMAP(metric='precomputed', n_neighbors=5, min_dist=0.0, random_state=42)
    embedding = reducer.fit_transform(arr)

    # Assuming `clusters` is an array of cluster labels
    plt.scatter(embedding[:, 0], embedding[:, 1], c=clusters, cmap='Spectral', s=5)
    plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP projection of the DNA dataset', fontsize=24)
    plt.show()
main()

