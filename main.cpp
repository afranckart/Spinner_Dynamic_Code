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

    int graine = 200;
    int n = 10;

    double lam = 0.90;
    double T0 = 2;
    double Tf = 0.003;  


    double J = -6.;
    double Jd = 1.;

    int Niter = 1000;
    int Nsimu = 20;

    int p = 12;


    std::vector<Spinner> spin = configuration_rand(n, n, graine);
    
    clock_t debut = clock();

    std::string add = path("simulation/recuit_multi_state/G200_10x10/", graine, lam, T0, Tf, J, Jd, Niter, Nsimu, n, n);

    //std::vector<std::vector<int>> A = Metropolis(spin, lam, T0, Tf, J, Jd, Niter, Nsimu, n, n, p);

    //print_data(A, add, n, n);

    //std::vector<std::vector<int>> A = read_data(add, n, n);

    //Recuit_allmeta(spin, A, add, lam, 0.1, Tf, J, Jd, Niter, Nsimu, 3, n, n, p);
  

    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}
