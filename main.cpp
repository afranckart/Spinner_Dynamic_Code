#include <iostream>
#include <omp.h>
#include "spinner_ppm.h"
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <cstring>
#include "spinner_CUDA.cuh"


int main() {

    clock_t debut = clock();

    #pragma omp parallel for num_threads(1)
    for(int j = 11; j < 31; j++){
        double i = j/100.;
        printf("%d\n", j);

        int num_threads = omp_get_num_threads();

        double L = 0.025;
        int nx = 1 * j;
        int ny = 1 * j;

        spinners_t spin;
        spinners_init(&spin, L, nx, ny, 1);
    
        double* H = H_init(L);
        //H_plot(H)

        double* HB = H_B_init(0, 0);
        //H_B_plot(HB);

        srand((unsigned int)time(NULL));

        std::string direc = "/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/UM_demo/";
        
        std::string add = direc + "ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_L0.025_Elow.txt";
        //std::string add = direc + "ppm_4x4_L0.025000_allmeta.txt";
        char* addspin = new char[add.length() + 1];
        std::strcpy(addspin, add.c_str());

        std::string direc2 = "/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/UM_demo/";

        std::string add0 = direc2 + "ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_T00.350000" ;
        char* spin0 = new char[add0.length() + 1];
        std::strcpy(spin0, add0.c_str());
        
        //print_Emin(&spin, spin0 , 1000);

        read_spinnersall(&spin, addspin, nx, ny, L);
        recuitN(&spin, H, HB, 0.35, 0.001, 0.95, 1000 * nx * ny, 1000, 8);
        //plot_E_mean(&spin, H, HB, i);
        //print_E_Histo(&spin, spin0, H, HB);
        print_dist(&spin, spin0, "EG_allmeta", H, HB, dist_EG);
        print_dist(&spin, spin0, "EL_allmeta", H, HB, dist_EL);
        print_dist(&spin, spin0, "H_allmeta", H, HB, dist_H);
        print_dist(&spin, spin0, "HI_allmeta", H, HB, dist_HI);
        
        //FILE* fichier = openfile_out(addspin);
        //print_spinners(&spin, fichier);
        //fclose(fichier);
        

        Finalisation_simu(&spin, H, HB);
        delete[] addspin;
        delete[] spin0;
    }
    
    
    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}

