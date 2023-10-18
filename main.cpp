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

    int graine = 200;
    int n = 10;

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

    std::string path = "C:\\Users\\axelf\\OneDrive - Universite de Liege\\Mémoire\\simulation\\from_meta_1\\";

    std::vector<Spinner> spin = configuration_rand(n, n, graine); 

    std::vector<int> U= read_angle(path + "G200_10x10_meta_1", n, n);
    for (int j = 0; j < n * n; j++) { spin[j].update_orientation(U[j]); }

    
    //std::vector<std::vector<int>> A = Metropolis(spin, 0.9, 10, 0.1, J, 1000, 10, n, n, p);

    //A = find_meta(spin, A, J, 1000, n, n, p);

    //print_data(A, path + "G200_10x10_meta_1_T010_meta", n, n);
    
    //print_dist(A, path + "G200_10x10_meta_1_T010_meta", n, n);

    //print_E(A, spin, path + "G200_10x10_meta_1_t010_meta", J, n, n);
    
    //std::vector<std::vector<int>> A = read_data( path + "G200_10x10_meta_1_T010_meta", n, n);

    //Recuit_allmeta(spin, A, path + "G200_10x10_meta_1_T010_meta", 0.9, 1, 0.01, J, 1000, 100, 1, n, n, p);
    
    std::vector<std::vector<int>> A = read_data(path + "G200_10x10_meta_1_T010_meta_meta_RC1_T01.000000", n, n);
    for (int i = 0; i < A.size(); i++)
    {
        std::vector<Spinner> spin = configuration_rand(n, n, graine);
        for (int j = 0; j < n * n; j++) { spin[j].update_orientation(A[i][j]); }
        std::cout << i << "\t" << E_total(spin, J, n, n) << std::endl;
    }


    


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
