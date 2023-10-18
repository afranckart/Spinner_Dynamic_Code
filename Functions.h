//
// Fonctions associée à la class Spinner
//

#ifndef UNTITLED8_FUNCTIONS_H
#define UNTITLED8_FUNCTIONS_H

#endif //UNTITLED8_FUNCTIONS_H

#include <vector>
class Spinner;
struct couplage;

/**
 * @brief make a convential structure path for file
 * @param add is the diretory path
 * @return std::string
 */
std::string path(std::string add,int graine, double lambda, double T0, double Tf,
    double J, double Jdouble, unsigned long int Niteration, int Nsimulation, int nx, int ny);

/**
 * @brief print all information of a configuration of spinner : index charge charge charge angle
 * @param add is the diretory path
 * @return std::string
 */
void print_spinner(std::vector<Spinner>& psin, std::string add);


/**
 * @brief print all angle configuration for all state in data
 * @param add is the diretory path
 * @return file with a state by line, all angle abd state multiplicity by colunm
 */
void print_data(std::vector<std::vector<int>>& data, std::string add, int nx, int ny);


/**
 * @brief read from a file, all information of a configuration of spinner : index charge charge charge angle
 * @param add is the diretory path
 * @return std::vector<Spinner>
 */
std::vector<Spinner> read_spinner(std::string add, int nx, int ny);

/**
 * @brief read from a file, all the angles 
 * @param add is the diretory path
 * @return std::vector<int>
 */
std::vector<int> read_angle(std::string add, int nx, int ny);


/**
 * @brief performe hamming distance between 2 lattice of spinner
 * @param A, B are 2 angle configuration of the same spinner lattice : 
 * @return int
 */
int dist(const std::vector<int>& A, const std::vector<int>& B);

/**
 * @brief change de value between 2 int
 * @param 2 int
 * @return int
 */
void change(int& a, int& b);

/**
 * @brief performe a distance beatween 2 line or collunm of distance matrice Dij
 * @param A, B are tow  line or collunmof distance matrice Dij
 * @return int
 */
int distline(std::vector<int>& A, std::vector<int>& B);


/**
 * @brief trie of distance matrice Dij by clustering algorithme
 * @param A is the  distance matrice Dij
 * @return void
 */
void tri(std::vector<std::vector<int>>& A);

void clustering(std::vector<std::vector<int>>& A);

/**
 * @brief performe the distance matrice Dij and the probability of have a distance d, and print this in 2 file
 * @param add is the diretory path
 * @return void
 */
void print_dist(std::vector<std::vector<int>>& data, std::string add, int nx, int ny);

/**
 * @brief read data print by print_data, 
 * @param add is the diretory path
 * @return return vector<vector<int>> of file with a state by line, all angle abd state multiplicity by colunm
 */
std::vector<std::vector<int>> read_data(std::string add,  int nx, int ny);

/**
 * @brief remove equale state of the output of métropolise ct abd add the multiplicity of the state at the laste colunm
 * @param databrut is the ouput of metropolis fct and this matrice with line as state and column as angle
 * @return return vector<vector<int>> of file with a state by line, all angle abd state multiplicity by colunm
 */
std::vector<std::vector<int>> remove_equal(std::vector<std::vector<int>>& databrut);



/**
 * @brief check if any modification of a angle of a lattice of sponner, cost energy
 * @param databrut is the ouput of metropolis fct and this matrice with line as state and column as angle
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * return true if metastable
 */
bool metastable(const std::vector<int>& A, const std::vector<Spinner>& spin0,  couplage J, int nx, int ny);

bool metastable(const std::vector<Spinner>& spin0, couplage J, int nx, int ny);

/**
 * @brief find and print metatsable state
 * @param databrut is the ouput of metropolis fct and this matrice with line as state and column as angle
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * return a file with the metastable state angle conguration
 */
void print_metastable(std::vector<std::vector<int>>& data, std::vector<Spinner>& spin0, std::string add, couplage J, int nx, int ny, int p);


/**
 * @brief print the state that any modification of a angle of a lattice of sponner, cost energy
 * @param databrut is the ouput of metropolis fct and this matrice with line as state and column as angle
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * return a file with the metastable state angle conguration
 */
void print_E(std::vector<std::vector<int>>& data, std::vector<Spinner>& spin0, std::string add, couplage J, int nx, int ny);


