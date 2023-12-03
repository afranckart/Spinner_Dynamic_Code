#include <iostream>
#include <omp.h>
#include "spinner_ppm.h"
#include <ctime>
#include <cmath>
#include <string>
#include <cstring>
#include "spinner_CUDA.cuh"


int main() {

    clock_t debut = clock();

    
     #pragma omp parallel for num_threads(5)
    for(int j = 50; j < 51; j++){
        double i = j/100.;

        int num_threads = omp_get_num_threads();

        double L = 0.025;
        int nx = 2;
        int ny = 2;

        spinners_t spin;
        spinners_init(&spin, L, nx, ny, 1);
    
        double* H = H_init(L);
        //H_plot(H)

        double* HB = H_B_init(0, 0);
        //H_B_plot(HB);

        srand((unsigned int)time(NULL));

        std::string direc = "/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/allmeta/";
        
        //std::string add = direc + "ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_L0.025_T0"+std::to_string(i) + "_recuit1000.txt";
        std::string add = direc + "ppm_2x2_L0.025000_allmeta.txt";
        char* addspin = new char[add.length() + 1];
        std::strcpy(addspin, add.c_str());


        std::string add0 = direc + "ppm_2x2_L0.025000_allmeta_E.txt";
        char* spin0 = new char[add0.length() + 1];
        std::strcpy(spin0, add0.c_str());
        

        read_spinnersall(&spin, addspin, nx, ny, L);
        print_E(&spin, spin0, H, HB);
        //print_dist(&spin, spin0, "EG_allmeta", H, HB, dist_EG);
        //print_dist(&spin, spin0, "EL_allmeta", H, HB, dist_EL);
        //print_dist(&spin, spin0, "H_allmeta", H, HB, dist_H);
        //print_dist(&spin, spin0, "HI_allmeta", H, HB, dist_HI);
        //recuitN(&spin, H, HB, i, 0.001, 0.95, 5 * nx * ny, 1000, 10);
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

