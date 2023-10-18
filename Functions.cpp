//
// Fonctions associée à la class Spinner
//
#include "Spinner.h"
#include "Spinner_Dynamic.h"
#include "Functions.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <map>

/* ----------- génération des adresse ------------*/
std::string path(std::string add, int graine, double lambda, double T0, double Tf,
    double J, double Jdouble, unsigned long int Niteration, int Nsimulation, int nx, int ny) {

    return add + "G" + std::to_string(graine) + "_" + std::to_string(nx) + "x" + std::to_string(nx) + "_lam" + std::to_string(lambda)
        + "_T0" + std::to_string(T0) + "_Tf" + std::to_string(Tf) + "_J" + std::to_string(J) + "_Jd"
        + std::to_string(Jdouble) + "_Nit" + std::to_string(Niteration) + "_Nsim" + std::to_string(Nsimulation);
}

/* ----------- écrit l'état de tous les spinner dans .txt ------------*/

void print_spinner(std::vector<Spinner> &spin, std::string add) {
    for (Spinner& s : spin)
    {
        s.print(add + ".txt");
    }
}

/* ----------- écrit l'état de tous les simulations de métropolis dans .txt ------------*/
/* ----------- ligne : états diff, colone : angle, dernière colonne : multiplicité ------------*/

void print_data(std::vector<std::vector<int>>& data, std::string add, int nx, int ny) {


    if (data[0].size() != nx * ny + 1) { std::cout << "ERROR in print_data : data[0].size() != nx * ny + 1" << std::endl; }

    std::ofstream fichier(add + ".txt");

    if (!fichier.is_open()) {
        std::cout << "Erreur : Impossible d'ouvrir le fichier." << add + ".txt" << std::endl;
    }
    else
    {
        for (int i = 0; i < data.size(); i++)
        {
            fichier << data[i][0];
            for (int j = 1; j < data[0].size(); j++)
            {
                fichier << "\t" << data[i][j];
            }
            if (i != data.size() -1) { fichier << std::endl; }
        }
    }
    fichier.close();
}

/* ----------- écrit l'état de tous les simulations de métropolis dans .txt ------------*/
/* ----------- ligne : états diff, colone : angle, dernière colonne : multiplicité ------------*/

std::vector<std::vector<int>> read_data(std::string add, int nx, int ny) {


    std::ifstream fichier(add + ".txt", std::ios::out);

    int N = nx * ny + 1;
    std::vector<std::vector<int>> data(0, std::vector<int>(N));

    if (!fichier.is_open()) {
        std::cout << "Erreur : Impossible d'ouvrir le fichier." << add + ".txt" << std::endl;
    }
    else
    {
        while (fichier.good())
        {
            std::vector<int> A(N);
            for (int i = 0; i < N; i++)
            {
                fichier >> A[i];
            }
            data.push_back(A); 
        }
    }
    return data;
}


/* ----------- lit l'état de tous les spinner depuis .txt ------------*/

std::vector<Spinner> read_spinner( std::string add, int nx, int ny) {

    std::ifstream fichier(add + ".txt");
    std::vector<Spinner> spin;

    if (!fichier) {
        std::cout << "Impossible d'ouvrir le fichier: " << add << std::endl;
    }
    else
    {
        spin.reserve(nx * ny);

        for (int i = 0; i < ny * nx; i++)
        {
            int index, angle;
            std::vector<int> site(3);
            fichier >> index;
            fichier >> site[0];
            fichier >> site[1];
            fichier >> site[2];
            fichier >> angle;
            spin.emplace_back(index, site, angle);
        }
    }
    fichier.close();
  
    return spin;
}

/* ----------- lit l'angle de tous les spinner depuis .txt ------------*/

std::vector<int> read_angle(std::string add, int nx, int ny) {
    std::ifstream fichier(add + ".txt");
    std::vector<int> angle(nx * ny);
    if (fichier.is_open()) {
        for (int i = 0; i < ny * nx; i++)
        {
            int a;
            fichier >> a;
            fichier >> a;
            fichier >> a;
            fichier >> a;
            fichier >> angle[i];
        }
    }
    else
    {
        std::cout << "read_print : Impossible d'ouvrir le fichier: " << add + ".txt" << std::endl;
    }
    fichier.close();

    return angle;
}

/* ----------- calcule la distance  entre 2 configuration------------*/
float dist(const std::vector<int>& A, const std::vector<int>& B) {
    size_t taille = A.size();
    float d = 0;
    for (size_t i = 0; i < taille - 1; i++) {
        int diff = (float)(B[i] - A[i]) * (B[i] - A[i]);
        switch (diff) {
        case 16:
            diff = 4;
            break;
        case 25:
            diff = 1;
            break;
        default:
            break;
        }
        d += diff;
    }
    return sqrt(d);
}

/* ----------- change ------------*/
void change(float& a, float& b)
{
    float c = a;
    a = b;
    b = c;
}

/* ----------- distance entre les lignes de la matrice des distance Dij ------------*/
float distline(std::vector<float>& A, std::vector<float>& B)
{
    float d = 0;;
    for (int i = 0; i < A.size(); i++) { d += (float)abs(A[i] - B[i]); }
    return d;
}

/* ----------- trie la matrice de distance ------------*/
void tri(std::vector<std::vector<float>>& A) 
{
    int n = A.size();
    std::vector<int> newindex(n);
    for (int i = 0; i < n - 1; i++) // parcourt les colonnes et cherche les lignes les plus proches
    {
        int k = i + 1;
        float d = distline(A[i], A[i+1]);
        for (int j = i + 2; j < n; j++)
        {
            if (distline(A[i], A[j]) < d) 
            {
                d = distline(A[i], A[j]);
                k = j;
            }
        }
        newindex[i + 1] = k;
        for (int j = 0; j < n; j++)
        {
            change(A[i + 1][j], A[newindex[i + 1]][j]);
        }
    }
    for (int i = 1; i < n - 1; i++) // change les colonnes de la même facon que les lignes
    {
        for (int j = 0; j < n; j++)
        {
            change(A[j][i], A[j][newindex[i]]);
        }
    }
}


/* ----------- trie la matrice de distance par clustering ------------*/
void clustering(std::vector<std::vector<float>>& A) // méthode de clustering
{
    int n = A.size();
    std::vector<int> newindex(n);
    for (int i = 0; i < n - 1; i++) // parcourt les colonnes et cherche les lignes les plus proches
    {
        int k = i + 1;
        int d = distline(A[i], A[i + 1]);
        for (int j = i + 2; j < n; j++)
        {
            if (distline(A[i], A[j]) < d)
            {
                d = distline(A[i], A[j]);
                k = j;
            }
        }
        newindex[i + 1] = k;
        for (int j = 0; j < n; j++)
        {
            change(A[i + 1][j], A[newindex[i + 1]][j]);
        }
    }
    for (int i = 1; i < n - 1; i++) // change les colonnes de la même facon que les lignes
    {
        for (int j = 0; j < n; j++)
        {
            change(A[j][i], A[j][newindex[i]]);
        }
    }
}

/* ----------- calcule la matrice de distance et la probabilité de trouver ces distance ------------*/

void print_dist(std::vector<std::vector<int>>& data, std::string add, int nx, int ny) {

    if (data[0].size() != nx * ny + 1) { std::cout << "ERROR in print_dist : data[0].size() != nx * ny + 1" << std::endl; }

    int N = data.size();
    int n = data[0].size() - 1; 
    std::vector<std::vector<float>> matrice(N, std::vector<float>(N));
    float dmax = 0;
    for (int i = 0; i < N ; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            float d = dist(data[i], data[j]);
            matrice[i][j] = d;
            matrice[j][i] = d;
            if (d > dmax) { dmax = d; } 
        }
    }

    std::vector<double> prob(dmax, 0);
    for (int i = 0; i < N; i++) // -1 car la dernière colonne c'est la multiplicité de l'état
    {
        for (int j = 0; j < i; j++)
        {
            if (matrice[i][j] != 0)
            {
                prob[(int)matrice[i][j] - 1] += data[i][n] * data[j][n];
            }
        }
    } 
    int normalisation = 0;
    for (int i = 0; i < prob.size(); i++) { normalisation += prob[i]; }
    
    tri(matrice);

    std::ofstream fileDist(add + "_MD" + ".txt", std::ios::out);

    if (!fileDist.is_open()) {
        std::cout << "Erreur print_dist : Impossible d'ouvrir le fichier." << add + ".txt" << std::endl;
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            
            fileDist << matrice[i][0];
            for (int j = 1; j < N; j++)
            {
                fileDist << "\t" << matrice[i][j];
            }
            if (i != N - 1) { fileDist << std::endl; }
        }
    }
    fileDist.close();

    std::ofstream fileProb(add + "_PD" + ".txt", std::ios::out);

    if ( !fileProb.is_open()) {
        std::cout << "Erreur print_dist : Impossible d'ouvrir le fichier." << add + ".txt" << std::endl;
    }
    else
    {
        for (int i = 0; i < prob.size(); i++)
        {
            if(prob[i] != 0){ fileProb << i + 1 << "\t" << prob[i] / normalisation << std::endl; }
        }
    }
    fileProb.close();
}

/* ----------- élimine les états équivalents, rend leurs multiplicité dans la dernière colonne ------------*/

std::vector<std::vector<int>> remove_equal(std::vector<std::vector<int>>& databrut)
{
    std::vector<std::vector<int>> data(0, std::vector<int>(databrut[0].size() + 1));

    for(std::vector<int>& state : databrut)
    {
        bool unique = true;
        for (std::vector<int>& stateunique : data)
        {
            for (int k = 0; k < state.size(); k++)
            {
                unique = false;
                if (state[k] != stateunique[k])
                {
                    unique = true;
                    break;
                }
            }
            if (!unique) // si il y est déjà, multiplicité++ et on passe au suivant
            {
                stateunique.back()++;
                break;
            }
        }

        if (unique)
        {
            std::vector<int> A = state;
            A.push_back(1);
            data.push_back(A);
        }
    }
    return data;
}

/* ----------- vérifie que un états est métastable : tous changement entraine une augementation d E ------------*/
bool metastable(const std::vector<int>& A, const  std::vector<Spinner>& spin0, couplage J, int nx, int ny)
{
   std::vector<Spinner> spin = spin0;
   bool stable = true;
   for (int i = 0; i < spin.size(); i++) { spin[i].update_orientation(A[i]); }
   for (int i = 0; i < spin.size(); i++) // size de spin0 et non de A, car il y a la multiplicté dans A aussi
   {
       int angle = spin[i].orientation();
       float Eref = E_local(spin, i, J, nx, ny);

       int anglep = angle + 1;
       if (anglep > 5) { anglep -= 6; }
       int anglem = angle - 1;
       if (anglem < 0) { anglem += 6; }

       spin[i].update_orientation(anglep);
       float Ep = E_local(spin, i, J, nx, ny);
       spin[i].update_orientation(anglem);
       float Em = E_local(spin, i, J, nx, ny);

       spin[i].update_orientation(angle);
       if (Em < Eref || Ep < Eref) {
           stable = false;
           break;
       }
   }
   return stable;
}

/* ----------- vérifie que un états est métastable : tous changement entraine une augementation d E ------------*/

bool metastable(const std::vector<Spinner>& spin0, couplage J, int nx, int ny) {
    std::vector<Spinner> spin = spin0;
    bool stable = true;
    const int max_angle = 5;
    const int min_angle = 0;

    for (size_t i = 0; i < spin.size(); i++) {
        int angle = spin[i].orientation();
        float Eref = E_local(spin, i, J, nx, ny);

        int anglep = (angle + 1) % 6;
        int anglem = (angle - 1 + 6) % 6;

        spin[i].update_orientation(anglep);
        float Ep = E_local(spin, i, J, nx, ny);
        spin[i].update_orientation(anglem);
        float Em = E_local(spin, i, J, nx, ny);

        spin[i].update_orientation(angle);
        if (Em < Eref || Ep < Eref) {
            stable = false;
            break;
        }
    }
    return stable;
}


/* ----------- vérifie que un états est métastable : tous changement entraine une augementation d E ------------*/
void print_metastable(std::vector<std::vector<int>>& data, std::vector<Spinner>& spin0, std::string add, couplage J, int nx, int ny, int p)
{
    std::ofstream file(add + "_meta" + ".txt", std::ios::out);

    if (!file.is_open()) {
        std::cout << "Erreur print_metastable : Impossible d'ouvrir le fichier." << add + ".txt" << std::endl;
    }
    else
    {
        std::vector<std::vector<int>> A = find_meta(spin0, data, J, 10 * nx * ny, nx, ny, p);
        for (std::vector<int>& state : A)
        {
            file << state[0];
            for (int i = 1; i < state.size() ; i++) { file << "\t" << state[i]; }
            if (state != data.back()) { file << std::endl; }
        }
    }
    file.close();
}

/*---------- - donne l histogramme de l'énergie des états ------------ */
void print_E(std::vector<std::vector<int>>&data, std::vector<Spinner>& spin0, std::string add, couplage J, int nx, int ny)
{
    if (data[0].size() != nx * ny + 1) { std::cout << "ERROR in print_E : data[0].size() != nx * ny + 1" << std::endl; }
    
    std::vector<int> E;
    E.reserve(data.size());
    std::vector<Spinner> spin = spin0;
    std::map<float, int> histogramme;
    for (std::vector<int>& state : data)
    {
        for (int i = 0; i < spin.size(); i++) { spin[i].update_orientation(state[i]); }
        histogramme[E_total(spin, J, nx, ny)]++;
    }
    std::ofstream file(add + "_E" + ".txt", std::ios::out);

    if (!file.is_open()) {
        std::cout << "Erreur print_E : Impossible d'ouvrir le fichier." << add + ".txt" << std::endl;
    }
    else
    {
        for (const auto& pair : histogramme) {
            file << pair.second << "\t" << pair.first << "\n";
        }
    }
    file.close();
}

