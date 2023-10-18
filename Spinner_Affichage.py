# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def print_matrice (n, direc, add):
    angles = np.loadtxt(direc + add + '.txt')
    matrice = np.loadtxt( direc + 'G200_10x10' + '.txt')

    a = 1  # paramètre de maille
    L = 1 * (n + 1)  # taille du graph
    largeur = 0.4
    hauteur = largeur / 5
    
    for s in range(angles.shape[0]):
        # Création de la figure
        fig, ax = plt.subplots()
        
        # on trace les Spinners
        y = a / 2

        for k in range(n):
            x = a / 2
            if (k % 2 == 1):
                x += a / 2

            for j in range(n):
                # Création du Spinner
                sites = matrice[k * n + j][1], matrice[k * n + j][2], matrice[k * n + j][3]
                angle = angles[s][k * n + j]

                for i in range(3):
                    if (sites[i] == 1):
                        rectangle = plt.Rectangle((x, y), largeur, hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Red', facecolor='Red')
                        ax.add_patch(rectangle)
                        rectangle = plt.Rectangle((x, y), largeur, -hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Red', facecolor='Red')
                        ax.add_patch(rectangle)
                    else:
                        rectangle = plt.Rectangle((x, y), largeur, hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Blue', facecolor='Blue')
                        ax.add_patch(rectangle)
                        rectangle = plt.Rectangle((x, y), largeur, -hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Blue', facecolor='Blue')
                        ax.add_patch(rectangle)
                x += a
            y += a * np.cos(30 * 3.14159 / 180)

        # Réglages des axes pour que le rectangle soit correctement visible
        ax.set_xlim(0, L)
        ax.set_ylim(L * np.cos(30 * 3.14259 / 180), 0)

        # Supprimer les marges
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

        # Afficher la figure
        plt.gca().set_aspect('equal', adjustable='box')
        plt.gca().set_axis_off()

        # Sauvegarder la figure au format PDF
        plt.savefig(direc + add + '_' + str(s) + '.pdf', format='pdf', bbox_inches='tight', dpi=300)

        # plt.show()
        plt.close()


def print_spinner( n, direc, add):

    matrice = np.loadtxt(direc + add + '.txt')
    print(matrice)
    a = 1  # paramètre de maille
    L = 1 * (n + 1)  # taille du graph
    largeur = 0.4
    hauteur = largeur / 5

    # Création de la figure
    fig, ax = plt.subplots()
    # on trace les Spinners
    y = a / 2

    for k in range(n):
        x = a / 2
        if (k % 2 == 1):
            x += a / 2

        for j in range(n):
            # Création du Spinner
            sites = matrice[k * n + j][1], matrice[k * n + j][2], matrice[k * n + j][3]
            angle = matrice[k * n + j][4]

            for i in range(3):
                if (sites[i] == 1):
                    rectangle = plt.Rectangle((x, y), largeur, hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Red', facecolor='Red')
                    ax.add_patch(rectangle)
                    rectangle = plt.Rectangle((x, y), largeur, -hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Red', facecolor='Red')
                    ax.add_patch(rectangle)
                else:
                    rectangle = plt.Rectangle((x, y), largeur, hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Blue', facecolor='Blue')
                    ax.add_patch(rectangle)
                    rectangle = plt.Rectangle((x, y), largeur, -hauteur,
                                                  angle=360 - (angle / 2 + i) * 120, edgecolor='Blue', facecolor='Blue')
                    ax.add_patch(rectangle)
            x += a
        y += a * np.cos(30 * 3.14159 / 180)

    # Réglages des axes pour que le rectangle soit correctement visible
    ax.set_xlim(0, L)
    ax.set_ylim(L * np.cos(30 * 3.14259 / 180), 0)

    # Supprimer les marges
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Afficher la figure
    plt.gca().set_aspect('equal', adjustable='box')
    plt.gca().set_axis_off()

    # Sauvegarder la figure au format PDF
    plt.savefig(direc + add  + '.pdf', format='pdf', bbox_inches='tight', dpi=300)

    # plt.show()
    plt.close()


direc = 'C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\from_meta_1\\'
add = 'G200_10x10_meta_1_T010_meta_meta_RC1_T01.000000'

n = 10

print_matrice(10, direc, add)