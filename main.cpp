#include <iostream>
#include <omp.h>
#include "Spinner.h"
#include "Functions.h"
#include "Spinner_Dynamic.h"
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>



int main() {

    clock_t debut = clock();

    int graine = 500;
    int n = 2;

    double lam = 0.90;
    double T0 = 2;
    double Tf = 0.003;

    couplage J;
    J.J = 6;
    J.Jdouble = -1;
    J.Joppose = 0.0;
    J.Jlow = 0.0;
    J.Jlowc = 0.0;

    int Niter = 1000;
    int Nsimu = 20;

    int p = 15;

    std::string path = "C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\++-\\";


    /*for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            std::vector<int> A = { -1,-2,-5 };
            Spinner B(1, A, i);

            std::vector<int> N = { 1,1,1 };
            Spinner V(1, N, j);

            std::cout << i << "\t" << j << "\t" << interaction(B, V, 5, JJ) << std::endl;
        }
    }*/

    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}
