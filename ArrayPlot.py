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
    
dirc = "C:/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/UM_demo/UM_demo_run3/UM_demo_run3_plot/"

a = 'EG'

file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist'+ a +'_allmeta'
plot_from_txt_and_save_png(file_path, a)
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist' + a +'_allmeta_ultra'
plot_from_txt_and_save_png(file_path, a + ' ultra')
a = 'EL'
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist' + a +'_allmeta'
plot_from_txt_and_save_png(file_path, a)
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist' + a +'_allmeta_ultra'
plot_from_txt_and_save_png(file_path, a + ' ultra')
a = 'H'
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist' + a +'_allmeta'
plot_from_txt_and_save_png(file_path, a)
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist' + a +'_allmeta_ultra'
plot_from_txt_and_save_png(file_path, a + ' ultra')
a = 'HI'
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) +'_T00.010000_TF0.0000001_L0.130000_dist'+ a +'_allmeta'
plot_from_txt_and_save_png(file_path, a)
file_path = dirc + 'ppm_' + str(20) +'x'+ str(20) + '_T00.010000_TF0.0000001_L0.130000_dist'+ a +'_allmeta_ultra'
plot_from_txt_and_save_png(file_path, a + ' ultra')
    
for i in range(20,21):
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta'
    plot_from_txt_and_save_png(file_path, a)
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta_ultra'
    plot_from_txt_and_save_png(file_path, a + ' ultra')

a = 'EL'
    
for i in range(2,30):
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta'
    plot_from_txt_and_save_png(file_path, a)
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta_ultra'
    plot_from_txt_and_save_png(file_path, a + ' ultra')

a = 'H'
    
for i in range(2,30):
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta'
    plot_from_txt_and_save_png(file_path, a)
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta_ultra'
    plot_from_txt_and_save_png(file_path, a + ' ultra')

a = 'EH'
    
for i in range(2,30):
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta'
    plot_from_txt_and_save_png(file_path, a)
    file_path = dirc + 'ppm_' + str(i) +'x'+ str(i) + '_T00.900000_L0.025000_dist'+ a +'_allmeta_ultra'
    plot_from_txt_and_save_png(file_path, a + ' ultra')
 
            
