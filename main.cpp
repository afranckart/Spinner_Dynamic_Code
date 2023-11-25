#include <iostream>
#include <omp.h>
#include "spinner_ppm.h"
#include <ctime>
#include <cmath>
#include <string>
#include "spinner_CUDA.cuh"


int main() {

    clock_t debut = clock();

    double L = 0.025;
    int nx = 2;
    int ny = 2;
    char add[] = "ppm_2x2_L0.025000_allmeta.txt";
    //char add[] = "ppm_2x2_L0.025000_allmeta.txt";
    //char addspin[] = "mnt\\c\\Users\\axelf\\OneDrive - Universite de Liege\\M�moire\\simulation\\ppm_5x5_L0.025_Elow.txt";

    spinners_t spin;
    spinners_init(&spin, L, nx, ny, 1);
    
    
    double* H = H_init(L);
    H_plot(H);

    double* HB = H_B_init( 0, 0);
    H_B_plot(HB);
    
    read_spinnersall(&spin, add, 2, 2, 0.025);
    //read_spinners(&spin, add);
    //printf("%d %d %d %d\n", spin.angles[0], spin.angles[1], spin.angles[2], spin.angles[3]);
    plotall(&spin);
    remove_equale(&spin);
    plotall(&spin);
    Finalisation_simu(&spin, H, HB); 
    
    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}

