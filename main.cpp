#include <iostream>
#include <omp.h>
//#include "Spinner.h"
//#include "Functions.h"
//#include "Spinner_Dynamic.h"
#include "spinner_ppm.h"
#include <ctime>
#include <cmath>
#include <string>
//#include <fstream>


int main() {

    clock_t debut = clock();

    double L = 0.025;
    int nx = 5;
    int ny = 5;
    char add[] = "C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\ppm_5x5";
    char addspin[] = "C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\ppm_5x5_L0.025_Elow.txt";

    spinners_t spin;
    spinners_init(&spin, L, nx, ny);
    
    
    double* H = H_init(L);
    H_plot(H);

    double* HB = H_B_init( 0, 0);
    H_B_plot(HB);
    
  
    Finalisation_simu(&spin, H, HB); 

    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}

