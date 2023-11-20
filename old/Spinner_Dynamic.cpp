#include <iostream>
#include <fstream>
#include <omp.h>
#include "Spinner.h"
#include "Functions.h"
#include "Spinner_Dynamic.h"
#include <ctime>
#include <cmath>
#include <vector>

/* ----------- Determine le nombre de voisins ------------*/
std::vector<int> neighbour(int kx, int ky, int nx, int ny) { //  on part de en haut à gauche et on commence à kx=ky=0


    std::vector<int> voisins(6, -1);


    if (kx != 0) {
        voisins[3] = ky * nx + kx - 1; // voisin de gauche
    }

    if (ky != 0)     //peut avoir des voisins au dessus
    {
        if (ky % 2 == 1)     //peut avoir un voisin au dessu à gauche
        {
            voisins[2] = (ky - 1) * nx + kx;
        }
        else if (ky % 2 == 0 && kx > 0) {
            voisins[2] = (ky - 1) * nx + kx - 1;
        }

        if (ky % 2 == 0)     //peut avoir un voisin au dessu à droite
        {
            voisins[1] = (ky - 1) * nx + kx;
        }
        else if (ky % 2 == 1 && kx < nx - 1) {
            voisins[1] = (ky - 1) * nx + kx + 1;
        }
    }

    if (kx != nx - 1) {
        voisins[0] = ky * nx + kx + 1;    //voisin de droite
    }


    if (ky != ny - 1)     //peut avoir des voisins en dessous
    {
        if (ky % 2 == 1)     //peut avoir un voisin en dessous à gauche
        {
            voisins[4] = (ky + 1) * nx + kx;
        }
        else if (ky % 2 == 0 && kx > 0) {
            voisins[4] = (ky + 1) * nx + kx - 1;
        }

        if (ky % 2 == 0)     //peut avoir un voisin en dessous à droite
        {
            voisins[5] = (ky + 1) * nx + kx;
        }
        else if (ky % 2 == 1 && kx < nx - 1) {
            voisins[5] = (ky + 1) * nx + kx + 1;
        }
    }
    return voisins; // droite, au dessus à droite, au dessus à gauche, gauche, en dessous à gauche, en dessous à droite

}


/* -- Renvoie l'énergie d'interaction entre 2 spinner voisin  -- */
//float interaction(Spinner& A, Spinner& voisins,  int nvoisins, couplage J) //nvoisins = 0 et s'incrémente dans le sens antihorloger
//{
//    float U = 0.; 
//    int alpha = A.orientation();
//    int theta = voisins.orientation();
//
//    int l = nvoisins * 2 ;
//
//    if (nvoisins % 2 == 0)
//    {
//        if (alpha % 2 == 0)
//        {
//            if (theta % 2 == 1) //  - -
//            {
//                U = J.J * A.charge()[(alpha + l) % 3] * voisins.charge()[(theta + l) % 3];
//            }
//            else // - =
//            {
//                U = J.Jdouble * A.charge()[(alpha + l) % 3]
//                    * (voisins.charge()[(theta + 1 + l) % 3] + voisins.charge()[(theta + 2 + l) % 3])
//                    + J.Joppose * A.charge()[(alpha + l) % 3] * voisins.charge()[(theta + l) % 3];
//            }
//        }
//        else
//        {
//            if (theta % 2 == 1) //  = -
//            {
//                U = J.Jdouble * voisins.charge()[(theta + l ) % 3]
//                    * (A.charge()[(alpha + 1 + l) % 3] + A.charge()[(alpha + 2 + l) % 3])
//                    + J.Joppose * voisins.charge()[(theta + l) % 3] * A.charge()[(alpha + l) % 3];
//            }
//            else // = =
//            {
//                U = J.Jlow * (A.charge()[(alpha + 2 + l) % 3] * voisins.charge()[(theta + 1 + l) % 3]
//                    + A.charge()[(alpha + 1 + l) % 3] * voisins.charge()[(theta + 2 + l) % 3]) + 
//                     J.Jlowc * (A.charge()[(alpha + 2 + l) % 3] * voisins.charge()[(theta + 2 + l) % 3]
//                    + A.charge()[(alpha + 1 + l) % 3] * voisins.charge()[(theta + 1 + l) % 3]);
//            }
//        }
//    }
//    else
//    {
//        if (alpha % 2 == 0)
//        {
//            if (theta % 2 == 1) //  = =
//            {
//                U = J.Jlow * (A.charge()[(alpha + 2 + l) % 3] * voisins.charge()[(theta + 1 + l) % 3]
//                    + A.charge()[(alpha + 1 + l) % 3] * voisins.charge()[(theta + 2 + l) % 3]) +
//                    J.Jlowc * (A.charge()[(alpha + 2 + l) % 3] * voisins.charge()[(theta + 2 + l) % 3]
//                        + A.charge()[(alpha + 1 + l) % 3] * voisins.charge()[(theta + 1 + l) % 3]);
//            }
//            else // = -
//            {
//                U = J.Jdouble * voisins.charge()[(theta + l) % 3]
//                    * (A.charge()[(alpha + 1 + l) % 3] + A.charge()[(alpha + 2 + l) % 3]) 
//                    + J.Joppose * A.charge()[(alpha  + l) % 3] * voisins.charge()[(theta + l) % 3];
//            }
//        }
//        else
//        {
//            if (theta % 2 == 1) //  - =
//            {
//                U = J.Jdouble * A.charge()[(alpha + l) % 3]
//                    * (voisins.charge()[(theta + 1 + l) % 3] + voisins.charge()[(theta + 2 + l) % 3]) 
//                    + J.Joppose * voisins.charge()[(theta  + l) % 3] * A.charge()[(alpha + l) % 3];
//            }
//            else // - -
//            {
//                U = J.J * A.charge()[(alpha + l) % 3] * voisins.charge()[(theta + l) % 3];
//            }
//        }
//    }
//    
//    return -U;
//}

/* -- Renvoie l'énergie d'interaction entre 2 spinner voisin  -- */
float interaction(Spinner& A, Spinner& voisins, int nvoisins, couplage J) {  //nvoisins = 0 et s'incrémente dans le sens antihorloger

    float H[6][6]; // [A][voisins]
    H[0][0] = -0.0179792;
    H[0][1] = 1.01101;
    H[0][2] = -0.0468777;
    H[0][3] = 0.988989;
    H[0][4] = -0.0753328;
    H[0][5] = -1.1096;
    H[1][0] = 0.25988;
    H[1][1] = -0.0179792;
    H[1][2] = -0.0197036;
    H[1][3] = -0.0753328;
    H[1][4] = -0.238725;
    H[1][5] = 0.0660521;
    H[2][0] = 0.0660521;
    H[2][1] = -1.1096;
    H[2][2] = 0.0221328;
    H[2][3] = -1.1096;
    H[2][4] = 0.0660521;
    H[2][5] = 1.25891;
    H[3][0] = -0.238725;
    H[3][1] = -0.0753328;
    H[3][2] = -0.0197036;
    H[3][3] = -0.0179792;
    H[3][4] = 0.25988;
    H[3][5] = 0.0660521;
    H[4][0] = -0.0753328;
    H[4][1] = 0.988989;
    H[4][2] = -0.0468777;
    H[4][3] = 1.01101;
    H[4][4] = -0.0179792;
    H[4][5] = -1.1096;
    H[5][0] = -0.0197036;
    H[5][1] = -0.0468777;
    H[5][2] = 0.150482;
    H[5][3] = -0.0468777;
    H[5][4] = -0.0197036;
    H[5][5] = 0.0221328; 

    int theta = A.orientation() - nvoisins;
    if (theta < 0) { theta += 6; }

    int beta = voisins.orientation() - nvoisins;
    if (beta < 0) { beta += 6; }

    return H[theta % 6][beta% 6];
}

/* -- Renvoie l'énergie d'interaction entre 1 spinner et ses 6 voisin  -- */
float E_local(std::vector<Spinner>& spin,  int index, couplage J, int nx, int ny) // opti : ne faire passer que les 6 voisinss
{
    int ky = (index - index % nx) / nx;
    int kx = index - ky * nx;

    std::vector<int> voisins = neighbour(kx, ky, nx, ny);

    float U = 0.;

    for (int i = 0; i < 6; i++) {
        if (voisins[i] != -1) { U += interaction(spin[index], spin[voisins[i]], i, J); }
    }
    
    return U;
}

/* ----------- Détermine l'énergie totale d'interaction ------------*/

float E_total(std::vector<Spinner>& spin, couplage J, int nx,  int ny) {

    float E = 0.;
    for (int j = 0; j < ny; j++) // colonne
    {
        for (int i = j % 2; i < nx; i += 2) // ligne
        {
            E += E_local(spin, i + j * nx, J, nx, ny);
        }
    }
    return  E; 
}


/*----Initialise une configuration de spinner------*/
std::vector<Spinner> configuration_rand(int nx, int ny, unsigned int graine) { // graine dans srand(time(NULL))

    srand(graine); // Initialistation

    int N = nx * ny;
    std::vector<Spinner> spin;

    spin.reserve(N); // réserve de la place

    for (int i = 0; i < N; i++)
    {
        std::vector<int> site(0);
        for (int j = 0; j < 3; j++)   // initialise la "charge" (+1 ou -1) de chaque site de l'hélice d'un spinneur, dans le sens trigo en commemcant par la charge sur R+
        {
            int charge = rand() % 2;
            if (charge == 0) { charge = -1; }
            site.push_back(charge);
        }

        int angle = rand() % 6;   // initialise le spinneur dans une orientation aléatoire : angle = angle*pi/3
        
        spin.emplace_back(i, site, angle);  // initialise le spinneur
    }

    return spin;
}

/*----Initialise une configuration de spinner de charge donner ------*/
std::vector<Spinner> configuration(int nx, int ny, unsigned int graine, int Q1, int Q2, int Q3) { // graine dans srand(time(NULL))

    srand(graine); // Initialistation

    int N = nx * ny;
    std::vector<Spinner> spin;

    spin.reserve(N); // réserve de la place

    std::vector<int> site = { Q1, Q2, Q3 };

    for (int i = 0; i < N; i++)
    {
        int angle = rand() % 6;   // initialise le spinneur dans une orientation aléatoire : angle = angle*pi/3

        spin.emplace_back(i, site, angle);  // initialise le spinneur
    }

    return spin;
}


/*----Recuit à T fixée------*/
void annealing_T_fixed(std::vector<Spinner>& spin, double T, couplage J, unsigned long int Niteration, int nx, int ny) {

    int N = nx * ny; 

    for (unsigned int i = 0; i < Niteration; i++) {
        
        int elem = rand() % N;  
       
        float U_old = E_local(spin, elem, J, nx, ny);
        int sign = rand() % 2; 
        if (sign == 0) { sign--; }

        int alpha_old = spin[elem].orientation();
        int alpha_new = alpha_old + sign; 
         
        if (alpha_new < 0) { alpha_new += 6; }
        if (alpha_new > 5) { alpha_new -= 6; }
       
        spin[elem].update_orientation(alpha_new); 
        float U_new = E_local(spin, elem, J, nx, ny); 
        
        if (U_new > U_old) {

            double temp = (double)rand() / (double)RAND_MAX;
            if (temp >= exp(-(U_new - U_old) / T)) { spin[elem].update_orientation(alpha_old); } // on rejette la modification
        }
    }
}


/*----Recuit Simulée, enregistre seulement les angles------*/
/* ----------- la dernière colonne est le # de fois que l'état est apparus ------------*/
std::vector<std::vector<int>> Metropolis(std::vector<Spinner>& spin, double lambda, double T0, double Tf, 
    couplage J, unsigned long int Niteration, int Nsimulation, int nx, int ny, int p){

    if (spin.size() != nx * ny) { std::cout << "ERROR in Metroplolis : spin.size() != nx * ny" << std::endl; }

    int N = nx * ny;
    std::vector<std::vector<int>> databrut(Nsimulation, std::vector<int>(N)); // l'utilisation des class est ouf niveau éfficasité
    srand(time(NULL));

    #pragma omp parallel for schedule(dynamic) num_threads(p)
    for (int j = 0; j < Nsimulation; j++)
    {
        std::vector<Spinner> spinThermalized = spin;
        for (double i = T0; i >= Tf; i *= lambda)
        {
            annealing_T_fixed(spinThermalized, i, J, Niteration, nx, ny);
        }
        for (int k = 0; k < N; k++)// extrait les angles
        {
            databrut[j][k] = spinThermalized[k].orientation();
        }
    }
    
    return  remove_equal(databrut) ;
}

/*----Recuit des états dejà recuit------*/
std::vector<std::vector<int>> Recuit(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, double lambda, double T0, double Tf,
    couplage J, unsigned long int Niteration, int Nsimulation, int nx, int ny, int p) {

    if (spin.size() != nx * ny) { std::cout << "ERROR in Recuit : spin.size() != nx * ny" << std::endl; }
    if (data[0].size() != nx * ny + 1) { std::cout << "ERROR in Recuit : data[0].size() != nx * ny + 1" << std::endl; }

    int N = nx * ny;
    std::vector<std::vector<int>> databrut(Nsimulation * data.size(), std::vector<int>(N)); // l'utilisation des class est ouf niveau éfficasité
    srand(time(NULL));

    std::vector<Spinner> spini = spin;

    for (int l = 0; l < data.size(); l++)
    {

        for (int i = 0; i < spini.size(); i++) { spini[i].update_orientation(data[l][i]); }
       
        #pragma omp parallel for schedule(dynamic) num_threads(p)
        for (int j = 0; j < Nsimulation; j++)
        {
            std::vector<Spinner> spinThermalized = spini;
            for (double i = T0; i >= Tf; i *= lambda)
            {
                annealing_T_fixed(spinThermalized, i, J, Niteration, nx, ny);
            }
            for (int k = 0; k < N; k++)// extrait les angles
            {
                databrut[j + Nsimulation * l ][k] = spinThermalized[k].orientation();
            } 
        }
    }

    return  remove_equal(databrut);
}


/*----Recuit à 0------*/
void annealing_0K(std::vector<Spinner>& spin, couplage J, unsigned long int Niteration, int nx, int ny) {

    int N = nx * ny;

    for (unsigned int i = 0; i < Niteration; i++) {

        int elem = rand() % N;    
       
        float U_old = E_local(spin, elem,  J, nx, ny);
        int sign = rand() % 2; 
        if (sign == 0) { sign--; }
        
        int alpha_old = spin[elem].orientation();
        int alpha_new = alpha_old + sign; 

        if (alpha_new < 0) { alpha_new += 6; }
        if (alpha_new > 5) { alpha_new -= 6; }
        
        spin[elem].update_orientation(alpha_new);
        float U_new = E_local(spin, elem, J, nx, ny); 
       
        if (U_new > U_old) { 
            spin[elem].update_orientation(alpha_old);
        }
    }
}

/* ----------- la dernière colonne est le # de fois que l'état est apparus ------------*/
std::vector<std::vector<int>> find_meta(std::vector<Spinner>& spin0, std::vector<std::vector<int>>& databrut,
     couplage J, int pasNiter, int nx, int ny, int p) {

    if (spin0.size() != nx * ny) { std::cout << "ERROR in find_meta : spin.size() != nx * ny" << std::endl; }
    if (databrut[0].size() != nx * ny + 1) { std::cout << "ERROR in find_meta : databrut[0].size() != nx * ny + 1" << std::endl; }

    int N = nx * ny;
    std::vector<std::vector<int>> data(databrut.size(), std::vector<int>(N)); // l'utilisation des class est ouf niveau éfficasité
    srand(time(NULL));

    #pragma omp parallel for schedule(dynamic) num_threads(p)
    for (int j = 0; j < data.size(); j++)
    {
        std::vector<Spinner> spinmeta = spin0;
        for (int i = 0; i < spinmeta.size(); i++) { spinmeta[i].update_orientation(databrut[j][i]); }
        
        while (!metastable(spinmeta, J, nx, ny))
        {
            annealing_0K(spinmeta, J, pasNiter, nx, ny);
        }
        for (int k = 0; k < N; k++)// extrait les angles
        {
            data[j][k] = spinmeta[k].orientation();
        }
    }

    return  remove_equal(data);
}

/*----Recuit des états déjà recuit------*/
void Recuit_all(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, std::string add, double lambda, double T0, double Tf,
    couplage J, unsigned long int Niteration, int Nsimulation, int RBSlevel, int nx, int ny, int p)
{
    std::vector<std::vector<int>> C = data;
    for (int i = 1; i <= RBSlevel; i++)
    {
        C = Recuit(spin, C, lambda, T0, Tf, J, Niteration, Nsimulation, nx, ny, p);

        print_E(C, spin, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), J, nx, ny);

        print_dist(C, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny, p);

        print_data(C, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny);

        print_metastable(C, spin, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), J, nx, ny, p);

    }
}

/*----Recuit des états meta------*/
void Recuit_allmeta(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, std::string add, double lambda, double T0, double Tf,
     couplage J, unsigned long int Niteration, int Nsimulation, int RBSlevel, int nx, int ny, int p)
{
    std::vector<std::vector<int>> C = data;
    for (int i = 1; i <= RBSlevel; i++)
    {
        C = Recuit(spin, C, lambda, T0, Tf, J, Niteration, Nsimulation, nx, ny, p);

        C = find_meta(spin, C, J, 5 * nx * ny, nx, ny, p);

        print_E(C, spin, add + "_meta_RC" + std::to_string(i) + "_T0" + std::to_string(T0), J, nx, ny);

        print_dist(C, add + "_meta_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny, p);

        print_data(C, add + "_meta_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny);

    }
}

