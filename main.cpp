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

    //#pragma omp parallel for num_threads(8)
    for(int j = 1; j < 2; j++){
        double i = j/1000.;
        //printf("%d\n", j);

        int num_threads = omp_get_num_threads();

        double L = 0.025;
        int nx = 5;
        int ny = 5;

        double T0 = 0.35;

        spinners_t spin;
        spinners_init(&spin, L, nx, ny, 1);
    
        double* H = H_init(L);
        //H_plot(H)

        double* HB = H_B_init(0, 0);
        //H_B_plot(HB);

        srand((unsigned int)time(NULL));

        std::string direc = "/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/UM_demo/";
        
        std::string add = direc + "ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_L0.025_Elow.txt";
        char* addspin = new char[add.length() + 1];
        std::strcpy(addspin, add.c_str());

        std::string direc2 = "/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/";

        std::string add0 = direc2 + "ppm_"+ std::to_string(nx) +"x"+ std::to_string(ny) + "_T0"+ std::to_string(T0) ;
        char* add2 = new char[add0.length() + 1];
        std::strcpy(add2, add0.c_str());

        std::string add00 = direc2 + "ppm_30x30_L0.025000_T0" + std::to_string(i) + "_recuit1000.txt";
        char* add3 = new char[add00.length() + 1];
        std::strcpy(add3, add00.c_str());

        //plot_MinMax(, L, nx, ny);
        
        print_Emin(L, nx, ny, add3, 10, 3);

        //read_spinnersall(&spin, add3, nx, ny, L);
        //recuitN(&spin, H, HB, T0, 0.001, 0.95, 1000 * nx * ny, 2000, 8);
        //plot_E_mean(&spin, H, HB, i);
        //print_E_Histo(&spin, spin0, H, HB);
        //double DEG = print_dist(&spin, add2, "EG_allmeta", H, HB, dist_EG);
        //double DEL = print_dist(&spin, add2, "EL_allmeta", H, HB, dist_EL);
        //double DH = print_dist(&spin, add2, "H_allmeta", H, HB, dist_H);
        //double DHI = print_dist(&spin, add2, "HI_allmeta", H, HB, dist_HI);

        //printf("%d\t%f\t%f\t%f\t%f\t%d\n", spin.nx,DEG, DEL, DH, DHI, spin.Ngrid);
        
        //FILE* fichier = openfile_out(add3);
        //print_spinners(&spin, fichier);
        //fclose(fichier);
        

        Finalisation_simu(&spin, H, HB);
        delete[] addspin;
        delete[] add2;
        delete[] add3;
    }
    
    
    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}

