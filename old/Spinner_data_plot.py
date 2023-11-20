# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm



def plot_data_from_files(path, pathout, name):
    # Charger les données
    probE = np.loadtxt(path + "_PDE.txt", delimiter='\t')
    prob = np.loadtxt(path + "_PD.txt", delimiter='\t')
    En = np.loadtxt(path + "_E.txt", delimiter='\t')
    matrice = np.loadtxt(path + "_MD.txt", delimiter='\t')
    matriceE = np.loadtxt(path + "_MDE.txt", delimiter='\t')
    
    # Tracer les données
    plt.figure()
    plt.scatter(En[:, 0],En[:, 1], marker='o',s=0.5)
    plt.xlabel("#E")
    plt.ylabel("E")
    plt.title(name)
    plt.savefig(pathout +'_E.png')  # Exporter le graphique En

    plt.figure()
    plt.scatter(probE[:, 0],probE[:, 1], marker='o',s=0.01)
    plt.xlabel("DE_alpha_beta")
    plt.ylabel("P(DE_alpha_beta)")
    plt.title(name)
    plt.savefig(pathout + '_PDE.png')  # Exporter le graphique probE

    plt.figure()
    plt.imshow(matriceE, aspect='auto', interpolation='none',cmap='gray')
    plt.colorbar()
    plt.title(name)
    plt.savefig(pathout +'_MDE.png')  # Exporter le graphique matriceE

    plt.figure()
    plt.scatter(prob[:, 0],prob[:, 1], marker='o')
    plt.xlabel("D_alpha_beta")
    plt.ylabel("P(D_alpha_beta)")
    plt.title(name)
    plt.savefig(pathout +'_PD.png')  # Exporter le graphique prob

    plt.figure()
    plt.imshow(matrice, aspect='auto', interpolation='none')
    plt.colorbar()
    plt.title(name)
    plt.savefig(pathout + '_MD.png')  # Exporter le graphique matrice


path = 'C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\++-2\\'
pathout = 'C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\++-2\\plot\\'
name = "5x5"
file = 'G500_'+name+'_ppm_meta'

plot_data_from_files(path + file, pathout + file, name)
