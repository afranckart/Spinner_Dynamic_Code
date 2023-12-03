import matplotlib.pyplot as plt
import numpy as np

def plot_from_txt_and_save_png(path):
    # Charger la matrice à partir du fichier texte avec NumPy
    data_matrix = np.loadtxt(path + '.txt', delimiter='\t')
    
    if  data_matrix.ndim == 0:
        return

    # Tracer l'array plot (heatmap)
    plt.imshow(data_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='dist')
    plt.xlabel('states')
    plt.ylabel('states')
    plt.title('dist HI')

    # Sauvegarder le graphique en tant qu'image PNG
    plt.savefig(path + '.png', format='png', dpi=1000)

    # Fermer la figure pour libérer la mémoire
    plt.close()

# Chemin vers le dossier contenant les fichiers texte
dirc = "C:/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/allmeta/"

file_path = dirc + 'ppm_3x3_L0.025000_distHI_allmeta'
plot_from_txt_and_save_png(file_path)
    
 
            
