

#ifndef UNTITLED8_SPINNER_DYNAMIC_H
#define UNTITLED8_SPINNER_DYNAMIC_H

#endif //UNTITLED8_SPINNER_DYNAMIC_H


#include <iostream>
#include <vector>

/**
 * @brief return a vector int with the index of neighbour spinner 
 * @param kx coordinate x
 * @param ky coordinate y
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 */
std::vector<int> neighbour(int kx, int ky, int nx, int ny);


float interaction(Spinner& A, Spinner& voisins, int nvoisins, couplage J);

float E_local(std::vector<Spinner>& spin, int index, couplage J, int nx, int ny);

/**
 * @brief return total energy of the lattice of spinner
 * @param spin is the vector with all the spinner of lattice
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 * @return totale energy : double
 */
float E_total(std::vector<Spinner>& spin, couplage J , int nx,  int ny);

/**
 * @brief return lattice of spinner 
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 * @param graine seed is the int used to seed random numbers with srand(graine)
 * @return a lattice of spinner : std::vector<Spinner>
 */
std::vector<Spinner> configuration_rand(int nx, int ny, unsigned int graine);

/**
 * @brief thermalize a lattice of spinner at a fixed T
 * @param spin is the vector with all the spinner of lattice
 * @param T is the temperature of simulation
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * @param Niteration is the number of random change made on the grid
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 * @return void, just the parameter spin is edited
 */
void annealing_T_fixed(std::vector<Spinner>& spin, double T, couplage J, unsigned  long int Niteration, int nx, int ny);

/**
 * @brief simulated annealing
 * @param spin is the initial vector with all the spinner of lattice
 * @param lambda is the temperature decreasing parameter Ti+1=lambda * Ti, we muste have 0<lambda<1
 * @param T0 is the initial temperature  of simulation
 * @param Tf is the final temperature  of simulation
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * @param Niteration is the number of random change made on the grid
 * @param Nsimulation si the number of performed simulated annealing
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 * @param p is the number of threats to performe this code with openMP
 * @return a vector of Nsimulation ( vector int with the final angle configuration after annealing  and the multiplicity ) : std::vector<std::vector<int>>
 */
std::vector<std::vector<int>> Metropolis(std::vector<Spinner>& spin, double lambda, double T0, double Tf, couplage J
	, unsigned long int Niteration, int Nsimulation, int nx, int ny, int p);


/**
 * @brief annealing some state that have already annealed
 * @param spin is the initial vector, give the information on charges
 * @param angle of annealed state input
 * @param lambda is the temperature decreasing parameter Ti+1=lambda * Ti, we muste have 0<lambda<1
 * @param T0 is the initial temperature  of simulation
 * @param Tf is the final temperature  of simulation
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * @param Niteration is the number of random change made on the grid
 * @param Nsimulation si the number of performed simulated annealing
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 * @param p is the number of threats to performe this code with openMP
 * @return a vector of Nsimulation ( vector int with the final angle configuration after annealing and the multiplicity ) : std::vector<std::vector<int>>
 */
std::vector<std::vector<int>> Recuit(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, double lambda, double T0, double Tf,
	 couplage J, unsigned long int Niteration, int Nsimulation, int nx, int ny, int p);


void annealing_0K(std::vector<Spinner>& spin, couplage J, unsigned long int Niteration, int nx, int ny);

std::vector<std::vector<int>> find_meta(std::vector<Spinner>& spin0, std::vector<std::vector<int>>& databrut,  couplage J,
	int pasNiter, int nx, int ny, int p);

/**
 * @brief annealing some state that have already annealed
 * @param spin is the initial vector, give the information on charges
 * @param angle of annealed state input
 * @param path of directoy
 * @param lambda is the temperature decreasing parameter Ti+1=lambda * Ti, we muste have 0<lambda<1
 * @param T0 is the initial temperature  of simulation
 * @param Tf is the final temperature  of simulation
 * @param J is the intercation coupling constante with aligned charges
 * @param Jdouble is the intercation coupling constante with no-aligned charges
 * @param Niteration is the number of random change made on the grid
 * @param Nsimulation si the number of performed simulated annealing
 * @param RSBLevel is the number of serial recuit
 * @param nx size of lattice along x
 * @param nx size of lattice along y
 * @param p is the number of threats to performe this code with openMP
 * @return a vector of Nsimulation ( vector int with the final angle configuration after annealing and the multiplicity ) : std::vector<std::vector<int>>
 */
void Recuit_all(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, std::string add, double lambda, double T0, double Tf,
	couplage J, unsigned long int Niteration, int Nsimulation, int RBSlevel, int nx, int ny, int p);

void Recuit_allmeta(std::vector<Spinner>& spin, std::vector<std::vector<int>>& data, std::string add, double lambda, double T0, double Tf,
	couplage J, unsigned long int Niteration, int Nsimulation, int RBSlevel, int nx, int ny, int p);