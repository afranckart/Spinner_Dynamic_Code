import matplotlib.pyplot as plt
import numpy as np

def plot_from_txt(path):
    # Charger la matrice à partir du fichier texte avec NumPy
    data_matrix = np.loadtxt(path, delimiter='\t')

    # Tracer l'array plot (heatmap)
    plt.imshow(data_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='dist')
    plt.xlabel('states')
    plt.ylabel('states')
    plt.title('dist EG')

    # Afficher le graphique
    plt.show()

dirc = "C:/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/"

plot_from_txt(dirc + 'ppm_10x10_T00.780000_L0.025000_distH.txt')
