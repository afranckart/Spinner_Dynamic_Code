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

    //char add[] = "ppm_2x2_L0.025000_allmeta.txt";

    
    
    for(int j = 1; j < 11; j++){
        double i = j/10.;

        int num_threads = omp_get_num_threads();

        double L = 0.025;
        int nx = 10;
        int ny = 10;

        spinners_t spin;
        spinners_init(&spin, L, nx, ny, 1);
    
        double* H = H_init(L);
        //H_plot(H);

        double* HB = H_B_init(0, 0);
        //H_B_plot(HB);

        srand((unsigned int)time(NULL));

        std::string direc = "/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/";
        
        std::string add = direc + "ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_L0.025_T0"+std::to_string(i) + "_recuit1000.txt";
        char* addspin = new char[add.length() + 1];
        std::strcpy(addspin, add.c_str());


        std::string add0 = direc + "EmaxMin/ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_L0.025_Elow.txt";
        char* spin0 = new char[add0.length() + 1];
        std::strcpy(spin0, add0.c_str());
        
        read_spinners(&spin, spin0);
        recuitN(&spin, H, HB, i, 0.001, 0.95, 5 * nx * ny, 1000, 10);
        FILE* fichier = openfile_out(addspin);
        print_spinners(&spin, fichier);
        fclose(fichier);
        

        Finalisation_simu(&spin, H, HB);
        delete[] addspin;
    }
    
    
    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}

