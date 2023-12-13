import matplotlib.pyplot as plt
import numpy as np

def plot_from_txt_and_save_png(path, name):
    # Charger la matrice à partir du fichier texte avec NumPy
    data_matrix = np.loadtxt(path + '.txt', delimiter='\t')
    
    if  data_matrix.ndim == 0:
        return

    # Tracer l'array plot (heatmap)
    plt.imshow(data_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='dist')
    plt.xlabel('states')
    plt.ylabel('states')
    plt.title('dist' + name)

    # Sauvegarder le graphique en tant qu'image PNG
    plt.savefig(path + '.png', format='png', dpi=1000)

    # Fermer la figure pour libérer la mémoire
    plt.close()

# Chemin vers le dossier contenant les fichiers texte
dirc = "C:/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/recuit_of_T0/run3/distEL/"

a = 'EL'
    
for i in range(10,100):
    file_path = dirc + 'ppm_10x10_T00.' + str(i) + '0000_L0.025000_dist'+ a +'_allmeta'
    plot_from_txt_and_save_png(file_path, a)
    file_path = dirc + 'ppm_10x10_T00.' + str(i) + '0000_L0.025000_dist'+ a +'_allmeta_ultra'
    plot_from_txt_and_save_png(file_path, a + ' ultra')
 
for i in range(1,10):
    file_path = dirc + 'ppm_10x10_T00.0' + str(i) + '0000_L0.025000_dist'+ a +'_allmeta'
    plot_from_txt_and_save_png(file_path, a)
    file_path = dirc + 'ppm_10x10_T00.0' + str(i) + '0000_L0.025000_dist'+ a +'_allmeta_ultra'
    plot_from_txt_and_save_png(file_path, a + ' ultra')

file_path = dirc + 'ppm_10x10_T01.000000_L0.025000_dist'+ a +'_allmeta'
plot_from_txt_and_save_png(file_path, a)
file_path = dirc + 'ppm_10x10_T01.000000_L0.025000_dist'+ a +'_allmeta_ultra'
plot_from_txt_and_save_png(file_path, a + ' ultra')
            
