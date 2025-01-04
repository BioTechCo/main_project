import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import squareform


def check_distance_matrix(distance_matrix):
    is_square = distance_matrix.shape[0] == distance_matrix.shape[1]
    print("Matrix is square:", is_square)

    diagonal_elements = np.diag(distance_matrix)
    are_diagonals_zero = np.all(diagonal_elements == 0)
    print("All diagonal elements are zero:", are_diagonals_zero)

    nan_count = distance_matrix.isna().sum().sum()
    print("Number of NaN values in the matrix:", nan_count)


def hierarchical_clustering(
    distance_matrix, range_min=2, range_max=31, cluster_number=None, out_path=None
):
    condensed_matrix = squareform(distance_matrix, force="tovector", checks=False)
    Z = linkage(condensed_matrix, method="ward")
    range_n_clusters = range(range_min, range_max + 1)
    silhouette_avg = []

    for n_clusters in range_n_clusters:
        labels = fcluster(Z, n_clusters, criterion="maxclust")
        silhouette_avg.append(
            silhouette_score(distance_matrix, labels, metric="precomputed")
        )
    if cluster_number:
        best_n_clusters = cluster_number
        print("chosen number of clusters:", best_n_clusters)
    else:
        best_n_clusters = range_n_clusters[np.argmax(silhouette_avg)]
        print("Best number of clusters:", best_n_clusters)

    # Create a figure with two subplots: one for the dendrogram and one for the silhouette scores
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Plot dendrogram on the left
    dendrogram(Z, ax=ax1)
    ax1.set_title("Dendrogram")
    ax1.set_xlabel("Sample index")
    ax1.set_ylabel("Distance")

    # Plot silhouette scores on the right
    ax2.plot(range_n_clusters, silhouette_avg, marker="o")
    ax2.set_title("Silhouette Score for Various Number of Clusters")
    ax2.set_xlabel("Number of Clusters")
    ax2.set_ylabel("Average Silhouette Score")
    ax2.set_xticks(range(2, range_max))

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(fig)

    labels = fcluster(Z, best_n_clusters, criterion="maxclust")
    dis_col = distance_matrix
    gene_names = dis_col.columns
    gene_clusters = pd.DataFrame({"gene": gene_names, "cluster": labels})
    return gene_clusters


def hierarchical_clustering_compare(
    distance_matrices, names, range_min=2, range_max=5, out_path=None
):
    # Check if exactly three distance matrices and three names are provided
    if len(distance_matrices) != 3 or len(names) != 3:
        raise ValueError(
            "Please provide exactly three distance matrices and three corresponding names."
        )

    # Prepare the plot
    plt.figure(figsize=(8, 6))

    # Loop through each distance matrix and calculate silhouette scores
    for idx, (distance_matrix, name) in enumerate(zip(distance_matrices, names)):
        condensed_matrix = squareform(distance_matrix, force="tovector", checks=False)
        Z = linkage(condensed_matrix, method="ward")
        range_n_clusters = range(range_min, range_max + 1)
        silhouette_avg = []

        for n_clusters in range_n_clusters:
            labels = fcluster(Z, n_clusters, criterion="maxclust")
            silhouette_avg.append(
                silhouette_score(distance_matrix, labels, metric="precomputed")
            )

        best_n_clusters = range_n_clusters[np.argmax(silhouette_avg)]
        print(f"Best number of clusters for {name}:", best_n_clusters)

        # Plot silhouette scores for the current distance matrix
        plt.plot(range_n_clusters, silhouette_avg, marker="o", label=name)

    # Finalize the plot
    plt.title("Silhouette Score for Various Number of Clusters")
    plt.xlabel("Number of Clusters")
    plt.ylabel("Average Silhouette Score")
    plt.xticks(range(range_min, range_max + 1))
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
