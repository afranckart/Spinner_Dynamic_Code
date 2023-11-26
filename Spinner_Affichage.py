# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import cairo
from matplotlib.patches import FancyArrow
from PIL import Image

def print_matrice(nx, ny, direc, add):
    angles = np.loadtxt(direc + add + '.txt')
    ##matrice = np.loadtxt(direc + addspin + '.txt')
    
    a = 1  # paramètre de maille
    Lx = 1 * (nx + 1)  # taille du graph
    Ly = 1 * (ny + 1)
    largeur = 0.4
    hauteur = largeur / 2.5

    resolution = 300  # DPI (points par pouce)

    for s in range(angles.shape[0]):
        # Création de la surface
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(Lx*resolution), int(Ly*resolution))
        #surface = cairo.PDF2urface(direc + add + '_' + str(s) + '.pdf', Lx*resolution, Ly*resolution)
        ctx = cairo.Context(surface)
        ctx.scale(resolution, resolution)

        # on trace les Spinners
        y = a / 2
        
        for k in range(ny):
            x = a / 2
            if k % 2 == 1:
                x += a / 2

            for j in range(nx):
                
                # Création du Spinner
                #sites = matrice[k * n + j][1], matrice[k * n + j][2], matrice[k * n + j][3]
                angle = angles[s][k * ny+ j]

                for i in range(3):
                    if i != 2:
                        ctx.set_source_rgb(1, 0, 0)  # Rouge
                    else:
                        ctx.set_source_rgb(0, 0, 1)  # Bleu

                    ctx.save()
                    ctx.translate(x, y)
                    ctx.rotate((360 - (angle / 2 + i) * 120) * 3.14259 / 180)

                    ctx.rectangle(0, -hauteur / 2, largeur, hauteur)
                    ctx.fill()
                    ctx.restore()

                x += a
            y += a * np.cos(30 * 3.14159 / 180)

        surface.write_to_png(direc + add + '_' + str(s) + '.png')

        # Conversion en PNG si nécessaire
        img = Image.open(direc + add + '_' + str(s) + '.png')
        img.save(direc + add + '_' + str(s) + '.png', 'PNG')

        surface.finish()


def print_spinner( nx, ny, direc, add):

    angles = np.loadtxt(direc + add + '.txt')
    ##matrice = np.loadtxt(direc + addspin + '.txt')
    
    a = 1  # paramètre de maille
    Lx = 1 * (nx + 1)  # taille du graph
    Ly = 1 * (ny + 1)
    largeur = 0.4
    hauteur = largeur / 2.5

    resolution = 300  # DPI (points par pouce)

    # Création de la surface
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(Lx*resolution), int(Ly*resolution))
    #surface = cairo.PDF2urface(direc + add + '_' + str(s) + '.pdf', Lx*resolution, Ly*resolution)
    ctx = cairo.Context(surface)
    ctx.scale(resolution, resolution)

    # on trace les Spinners
    y = a / 2
        
    for k in range(ny):
        x = a / 2
        if k % 2 == 1:
            x += a / 2

        for j in range(nx):
                
            # Création du Spinner
            #sites = matrice[k * n + j][1], matrice[k * n + j][2], matrice[k * n + j][3]
            angle = angles[k * ny+ j]

            for i in range(3):
                if i != 2:
                    ctx.set_source_rgb(1, 0, 0)  # Rouge
                else:
                    ctx.set_source_rgb(0, 0, 1)  # Bleu

                ctx.save()
                ctx.translate(x, y)
                ctx.rotate((360 - (angle / 2 + i) * 120) * 3.14259 / 180)

                ctx.rectangle(0, -hauteur / 2, largeur, hauteur)
                ctx.fill()
                ctx.restore()

            x += a
        y += a * np.cos(30 * 3.14159 / 180)

    surface.write_to_png(direc + add +  '.png')

    # Conversion en PNG si nécessaire
    img = Image.open(direc + add + '.png')
    img.save(direc + add + '.png', 'PNG')

    surface.finish()
        

def print_spinner_arrow(n, direc, add): #pour un ppm
    matrice = np.loadtxt(direc + add + '.txt')
  
    a = 1  # paramètre de maille
    L = 1 * (n + 1)  # taille du graph
    largeur = 0.32
    hauteur = largeur / 5

    # Création de la figure
    fig, ax = plt.subplots()
    # on trace les Spinners
    y = a / 2

    for k in range(n):
        x = a / 2
        if k % 2 == 1:
            x += a / 2

        for j in range(n):
            # Création du Spinner
            sites = matrice[k * n + j][1], matrice[k * n + j][2], matrice[k * n + j][3]
            angle = matrice[k * n + j][4]

            angle = angle + 1
            
            arrow_direction = 360 - (angle / 2) * 120
            arrow = FancyArrow(x, y, -largeur* np.cos(np.radians(arrow_direction)), -largeur* np.sin(np.radians(arrow_direction)), 
                               width=hauteur, head_length=0.12, head_width=0.2 , edgecolor='black', facecolor='black')
            ax.add_patch(arrow)
            
            arrow_opposite_direction = arrow_direction + 180  # Direction opposée en ajoutant 180 degrés
            opposite_x = x - largeur * np.cos(np.radians(arrow_opposite_direction))
            opposite_y = y - largeur * np.sin(np.radians(arrow_opposite_direction))
            plt.plot([x, opposite_x], [y, opposite_y], color='black', linewidth=5)  # Dessiner une ligne


            x += a
        y += a * np.cos(30 * np.pi / 180)

    # Réglages des axes pour que le rectangle soit correctement visible
    ax.set_xlim(0, L)
    ax.set_ylim(L * np.cos(30 * np.pi / 180), 0)

    # Supprimer les marges
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Afficher la figure
    plt.gca().set_aspect('equal', adjustable='box')
    plt.gca().set_axis_off()

    # Sauvegarder la figure au format JPG
    plt.savefig(direc + add + '_arrow.jpg', format='jpg', bbox_inches='tight', dpi=300)

    # plt.show()
    plt.close()


direc = 'C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\'
for i in range(2,11):
    add = 'ppm_' + str(i) + 'x' + str(i) + '_L0.025_Ehight' 
    #add = 'ppm_5x5_L0.025000_distance5_all' 

    n = 2

    print_spinner(i, i, direc, add)
