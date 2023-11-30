import matplotlib.pyplot as plt
import numpy as np

def plot_from_txt_and_save_png(path):
    # Charger la matrice à partir du fichier texte avec NumPy
    data_matrix = np.loadtxt(path + '.txt', delimiter='\t')

    # Tracer l'array plot (heatmap)
    plt.imshow(data_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='dist')
    plt.xlabel('states')
    plt.ylabel('states')
    plt.title('dist EG')

    # Sauvegarder le graphique en tant qu'image PNG
    plt.savefig(path+ '.png', format='png')

    # Fermer la figure pour libérer la mémoire
    plt.close()

# Chemin vers le dossier contenant les fichiers texte
dirc = "C:/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/"

for i in range(67,68,1):
    file_path = dirc + 'ppm_10x10_T00.' + str(i) + '0000_L0.025000_distEG'
    plot_from_txt_and_save_png(file_path)

    
 
            
