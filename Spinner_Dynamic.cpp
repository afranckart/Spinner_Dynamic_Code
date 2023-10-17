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


/* -- Renvoie l'énergie d'interaction d'un spinner avec le système -- */

double interaction( std::vector<Spinner>& spin, int i_elem, double J, double Jdouble, int nx, int ny) { // Jdouble est l'intéraction quand les charges ne sont pas alignée

    int elem = i_elem;

    // positions x et y du spinneur d'interet
    int ky = (elem - elem % nx) / nx;
    int kx = elem - ky * nx;

    std::vector<int> voisins = neighbour(kx, ky, nx, ny);
    
    double U = 0;
    int alpha = spin[elem].orientation(); // alpha apparentient à {0 1 2 3 4 5}
    
    if (alpha % 2 == 0) // on considère les voisins : droite, au dessus gauche, inférieur gauche = voisins 0 2 4
    {
        if (voisins[0] != -1) // vérifie que le voisin éxiste, droite
        {
            if (spin[voisins[0]].orientation() % 2 == 0) { // non alignée, angle->charge 0->1,2 2->0,1 4->0,2
                U -= Jdouble * spin[elem].charge()[ alpha % 3] *
                    (spin[voisins[0]].charge()[(spin[voisins[0]].orientation() + 1) % 3]
                        + spin[voisins[0]].charge()[ (spin[voisins[0]].orientation() + 2) % 3]);
            }
            else { // alignée, angle->charge 1->1 ou 3->0 ou 5->2 
                U -= J * spin[elem].charge()[ alpha % 3]
                    * spin[voisins[0]].charge()[ spin[voisins[0]].orientation() % 3];
            }
        }

        if (voisins[2] != -1) // vérifie que le voisin éxiste, au dessus gauche
        {
            if (spin[voisins[2]].orientation() % 2 == 0) { // non alignée, angle->charge 0->2,0 2->1,2 4->0,1 
                U -= Jdouble * spin[elem].charge()[( alpha + 1) % 3] *
                    (spin[voisins[2]].charge()[ spin[voisins[2]].orientation() % 3]
                        + spin[voisins[2]].charge()[(spin[voisins[2]].orientation() + 2 ) % 3]);
            }
            else { // alignée, angle->chage 1->2 3->1 5->0
                U -= J * spin[elem].charge()[( alpha + 1) % 3]
                    * spin[voisins[2]].charge()[(spin[voisins[2]].orientation() + 1) % 3];
            }
        }

        if (voisins[4] != -1) // vérifie que le voisin éxiste, inférieure gauche
        {
            if (spin[voisins[4]].orientation() % 2 == 0) { // non alignée, angle->charge 0->0,1 2->0,2 4->1,2
                U -= Jdouble * spin[elem].charge()[( alpha + 2) % 3] *
                    (spin[voisins[4]].charge()[(spin[voisins[4]].orientation() + 1) % 3]
                        + spin[voisins[4]].charge()[spin[voisins[4]].orientation() % 3]);
            }
            else { // alignée , angle->charge 1->0 3->2 5->1
                U -= J * spin[elem].charge()[( alpha + 2) % 3]
                    * spin[voisins[4]].charge()[(spin[voisins[4]].orientation() + 2) % 3];
            }
        }
    }
    else  // on considère les voisins : gauche, au dessus droite, inférieur droite = voisins 3 1 5
    {
        if (voisins[3] != -1) // vérifie que le voisin éxiste, gauche
        {
            if (spin[voisins[3]].orientation() % 2 == 1) { // non alignée, angle->charge 1->2,0 3->2,1 5->0,1
                U -= Jdouble * spin[elem].charge()[ alpha % 3] *
                    (spin[voisins[3]].charge()[(spin[voisins[3]].orientation() + 1) % 3]
                        + spin[voisins[3]].charge()[(spin[voisins[3]].orientation() + 2) % 3]);
            }
            else { // alignée, angle->charge 0->0 ou 2->2 ou 4->1 
                U -= J * spin[elem].charge()[ alpha % 3]
                    * spin[voisins[3]].charge()[spin[voisins[3]].orientation() % 3];
            }
        }

        if (voisins[1] != -1) // vérifie que le voisin éxiste, au dessus droite
        {
            if (spin[voisins[1]].orientation() % 2 == 1) { // non alignée, angle->charge 1->2,1 3->0,1 5->2,0 
                U -= Jdouble * spin[elem].charge()[ ( alpha + 2 ) % 3] *
                    (spin[voisins[1]].charge()[spin[voisins[1]].orientation() % 3]
                        + spin[voisins[1]].charge()[(spin[voisins[1]].orientation() + 1) % 3]);
            }
            else { // alignée, angle->charge 0->2 ou 2->1 ou 4->0 
                U -= J * spin[elem].charge()[( alpha + 2) % 3]
                    * spin[voisins[1]].charge()[(spin[voisins[1]].orientation() + 2) % 3];
            }
        }

        if (voisins[5] != -1) // vérifie que le voisin éxiste, inférieur droite
        {
            if (spin[voisins[5]].orientation() % 2 == 1) { // non alignée, angle->charge 1->0,1 3->0,2 5->1,2 
                U -= Jdouble * spin[elem].charge()[(alpha + 1) % 3] *
                    (spin[voisins[5]].charge()[spin[voisins[5]].orientation() % 3]
                        + spin[voisins[5]].charge()[(spin[voisins[5]].orientation() + 2) % 3]);
            }
            else { // alignée, angle->charge 0->1 ou 2->0 ou 4->2 
                U -= J * spin[elem].charge()[(alpha + 1) % 3]
                    * spin[voisins[5]].charge()[(spin[voisins[5]].orientation() + 1) % 3];
            }
        } 
    }

    return U;
}

/* ----------- Détermine l'énergie totale d'interaction ------------*/

double E_total(std::vector<Spinner>& spin, double J, double Jdouble, int nx, int ny) {

    double E = 0;
    int N = nx * ny;
    for (int i = 0; i < N; i++)
    {
        E += interaction(spin, i, J, Jdouble, nx, ny);
    }
    return  E; // Attention on compte plusieure fois les paires.
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


/*----Recuit à T fixée------*/
void annealing_T_fixed(std::vector<Spinner>& spin, double T, double J, double Jdouble, unsigned long int Niteration, int nx, int ny) {

    int N = nx * ny; 

    for (unsigned int i = 0; i < Niteration; i++) {
        
        int elem = rand() % N;     // choisis un spinneur au hasard
       
        double U_old = interaction(spin, elem, J, Jdouble, nx, ny);
        int sign = rand() % 2; // signe du changement aléatoire d'angle
        if (sign == 0) { sign--; }

        int alpha_old = spin[elem].orientation();
        int alpha_new = alpha_old + sign; // angle après changement aléatoire de +/- 1
         
        if (alpha_new < 0) { alpha_new += 6; }
        if (alpha_new > 5) { alpha_new -= 6; }
       
        spin[elem].update_orientation(alpha_new); 
        double U_new = interaction(spin, elem, J, Jdouble, nx, ny); // Energie après rotation
        
        if (U_new > U_old) {

            double temp = rand() / (double)RAND_MAX;
            if (temp >= exp(-(U_new - U_old) / T)) { spin[elem].update_orientation(alpha_old); } // on rejette la modification
        }
    }
}


/*----Recuit Simulée, enregistre seulement les angles------*/
/* ----------- la dernière colonne est le # de fois que l'état est apparus ------------*/
std::vector<std::vector<int>> Metropolis(std::vector<Spinner>& spin, double lambda, double T0, double Tf, 
    double J, double Jdouble, unsigned long int Niteration, int Nsimulation, int nx, int ny, int p){

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
            annealing_T_fixed(spinThermalized, i, J, Jdouble, Niteration, nx, ny);
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
    double J, double Jdouble, unsigned long int Niteration, int Nsimulation, int nx, int ny, int p) {

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
                annealing_T_fixed(spinThermalized, i, J, Jdouble, Niteration, nx, ny);
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
void annealing_0K(std::vector<Spinner>& spin, double J, double Jdouble, unsigned long int Niteration, int nx, int ny) {

    int N = nx * ny;

    for (unsigned int i = 0; i < Niteration; i++) {

        int elem = rand() % N;     // choisis un spinneur au hasard
       
        double U_old = E_total(spin, J, Jdouble, nx, ny);
        int sign = rand() % 2; // signe du changement aléatoire d'angle
        if (sign == 0) { sign--; }
        
        int alpha_old = spin[elem].orientation();
        int alpha_new = alpha_old + sign; // angle après changement aléatoire de +/- 1

        if (alpha_new < 0) { alpha_new += 6; }
        if (alpha_new > 5) { alpha_new -= 6; }
        
        spin[elem].update_orientation(alpha_new);
        double U_new = E_total(spin, J, Jdouble, nx, ny); // Energie après rotation
       
        if (U_new > U_old) { // a vérfier qu il y a tjr une barri§re entre 2 angles
            spin[elem].update_orientation(alpha_old);
        }
    }
}

/* ----------- la dernière colonne est le # de fois que l'état est apparus ------------*/
std::vector<std::vector<int>> find_meta(std::vector<Spinner>& spin0, std::vector<std::vector<int>>& databrut, double J, double Jdouble, int pasNiter, int nx, int ny, int p) {

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
        
        while (!metastable(spinmeta, J, Jdouble, nx, ny))
        {
            annealing_0K(spinmeta, J, Jdouble, pasNiter, nx, ny);
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
    double J, double Jdouble, unsigned long int Niteration, int Nsimulation, int RBSlevel, int nx, int ny, int p)
{
    std::vector<std::vector<int>> C = data;
    for (int i = 1; i <= RBSlevel; i++)
    {
        C = Recuit(spin, C, lambda, T0, Tf, J, Jdouble, Niteration, Nsimulation, nx, ny, p);

        print_E(C, spin, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), J, Jdouble, nx, ny);

        print_dist(C, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny);

        print_data(C, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny);

        print_metastable(C, spin, add + "_RC" + std::to_string(i) + "_T0" + std::to_string(T0), J, Jdouble, nx, ny, p);

    }
}

/*----Recuit des états meta------*/
void Recuit_allmeta(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, std::string add, double lambda, double T0, double Tf,
    double J, double Jdouble, unsigned long int Niteration, int Nsimulation, int RBSlevel, int nx, int ny, int p)
{
    std::vector<std::vector<int>> C = data;
    for (int i = 1; i <= RBSlevel; i++)
    {
        C = Recuit(spin, C, lambda, T0, Tf, J, Jdouble, Niteration, Nsimulation, nx, ny, p);

        C = find_meta(spin, C, J, Jdouble, 5 * nx * ny, nx, ny, p);

        print_E(C, spin, add + "_meta_RC" + std::to_string(i) + "_T0" + std::to_string(T0), J, Jdouble, nx, ny);

        print_dist(C, add + "_meta_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny);

        print_data(C, add + "_meta_RC" + std::to_string(i) + "_T0" + std::to_string(T0), nx, ny);

    }
}
