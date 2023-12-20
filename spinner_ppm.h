#pragma once

#ifndef UNTITLED8_SPINNER_PPM_H
#define UNTITLED8_SPINNER_PPM_H

#endif //UNTITLED8_SPINNER_PPM_H

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#define PI 3.14159265358979323846 
#define PI3 3.14159265358979323846/3. 
#define MU0 0.0000012566370614 //[SI]
#define HREF 0.274348 // interaction d-d pour m = 1  et L =0.025 dans le cas ou l'interaction est max [SI]

#define SIZE_H 6*6
#define SIZE_NEIGHBOUR 6
#define R 0.008 //distance entre le centre du dipo  le et le centre du spinner [m]
#define STRING_MAX 256

#define GETELEMENT(arr, i, v) arr[ i * 6 + v]
#define SETELEMENT(arr, i, v, value) arr[ i * 6 + v] = value

typedef struct spinners {
	
	int nx;
	int ny;
	int Ngrid;
	double L;
	int* angles;

} spinners_t;

typedef struct matrice_line{
	double* col;
	int pos;
}matrice_line_t;

typedef struct matrice{
	int N;
	matrice_line_t* line;
}matrice_t;

typedef struct cluster{
	int size;
	int* pos;
	bool notmerged;
}cluster_t;

typedef struct tree{
	int N;
	cluster_t* cluster;
}tree_t;





/********************************************************************************/
/*                                                                              */
/*                            comput Energy dipole fct                          */
/*                                                                              */
/********************************************************************************/

/**
 * @brief compute a dipole dipole interaction 
 * 
 * @param  m1_x [IN] composante x of dipole 1
 * 
 * @param  m1_y [IN] composante y of dipole 1
 * 
 * @param  m2_x [IN] composante x of dipole 2
 * 
 * @param  m2_y [IN] composante y of dipole 2
 * 
 * @param  r1_x [IN] composante x of position of dipole 1
 * 
 * @param  r1_y [IN] composante y of position of dipole 1
 * 
 * @param  r2_x [IN] composante x of position of dipole 2
 * 
 * @param  r2_y [IN] composante y of position of dipole 2
 * 
 * @return energy of dipole dipole interaction 
 */
double Udipole(double m1_x, double m1_y, double m2_x, double m2_y, double r1_x, double r1_y, double r2_x, double r2_y);


/**
 * @brief compute a dipole B-fiedl interaction
 *
 * @param  mx [IN] composante x of dipole
 *
 * @param  my [IN] composante y of dipole 
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 *
 * @return energy of dipole B-field interaction -mu.B
 */
double UB(double mx, double my, double bx, double by);

/**
 * @brief compute a ppm spinner spinner interaction
 *
 * @param  theta [IN] angle of spinner of interest
 *
 * @param  thetav [IN]  angle of niegtboor of spinner of interest
 *
 * @param  L [IN] distance between 2 spinners
 *
 * @param  nv [IN] position of neigtboor : right 1 to Below right 5 by counterclockwise
 *
 * @return energy of ppm spinner spinner interaction normalised by HERF
 */
double Uspinner(int theta, int thetav, double L, int nv);


/**
 * @brief compute a ppm spinner B-field interaction
 *
 * @param  theta [IN] angle of spinner of interest
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 *
 * @return energy of ppm spinner B-field interaction normalised by HERF
 */
double UBspinner(int theta, double bx, double by);

/**
 * @brief compute the interaction tensor by dipole dipole interaction
 *
 * @param  L [IN] distance between 2 spinner
 *
 * @return interaction tensor of the right neghitbourg
 */
double* H_init(double L);

/**
 * @brief compute the interaction tensor  by dipole B-field interaction
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 *
 * @return interaction tensor of the right neghitbourg
 */
double* H_B_init(double bx, double by);

/**
 * @brief plot interaction tensor 
 *
 * @param  H [IN]  interaction tensor dipole dipole
 */
void H_plot(double* H);


/**
 * @brief plot interaction tensor
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void H_B_plot(double* HB);

/********************************************************************************/
/*                                                                              */
/*          Initialisation spinner lattice and lattice fct                      */
/*                                                                              */
/********************************************************************************/

/**
 * @brief intitialize a spinner lattice structure
 *
 * @param  spin [OUT]  spinner_t
 * 
 * @param  L [IN]  distance between tow spinners
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  Ngrid [IN]  number de grid
 */
void spinners_init(spinners_t* spin, double L, int nx, int ny, int Ngrid);

/**
 * @brief free interaction tensor and spinners_t
 *
 * @param  spin [INOUT]  spinner_t
 * 
 * @param  H [INOUT]  interaction tensor
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void Finalisation_simu(spinners_t* spin, double* H, double* HB);

/**
 * @brief comput the index of spinner neighbour
 *
 * @param  spin [IN]  spinner lattice
 * 
 * @param  index [IN]  index of interesst spinner
 * 
 * @return array with neighbour position
 */
int* neighbour( spinners_t* spin, int index);

/********************************************************************************/
/*                                                                              */
/*                           Energy of spinner lattice                          */
/*                                                                              */
/********************************************************************************/


/**
 * @brief comput the Energy betwenn a spinner and all neighbour
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  index [IN]  index of interesst spinner
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return local Energy
 */
double E_local(spinners_t* spin, int index, double* H, double* HB, int offset);


/**
 * @brief comput the total Energy of spinners grid
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return total Energy
 */
double E_total(spinners_t* spin, double* H, double* HB, int offset);

/**
 * @brief detremine if spinner configuration is metastable
 * 
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return true if spinner configuration is metastable, or false
 */
bool metastable(spinners_t* spin, double* H, double* HB, int offset);


/********************************************************************************/
/*                                                                              */
/*                                      utilitaire                              */
/*                                                                              */
/********************************************************************************/


/**
 * @brief interverted tow int
 *
 * @param  A [INOUT]   int that will be interverted 
 * 
 * @param  B [INOUT]   int that will be interverted 
 */
void change(int* A, int* B);

/**
 * @brief open a file in w write and erase mode
 *
 * @param  add [IN]  pathfile
 *
 * @return a FILE*
 */
FILE* openfile_out(char* add);


/**
 * @brief open a file in read mode
 *
 * @param  add [IN]  pathfile
 *
 * @return a FILE*
 */
FILE* openfile_in(char* add);

/**
 * @brief print array spin->angles in a FILE
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  fichier [IN]  file where write
 */
void print_spinners(spinners_t* spin, FILE* fichier);


/**
 * @brief print a matrice sizex * sizey, a arrayr of pointeurs
 *
 * @param  matrice [IN] array of pointeur of line of the matrice
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  sizex [IN]  horizontale size of matrice
 * 
 * @param  sizey [IN]  verticale size of matrice
 * 
 * @param  fichier [IN]  file where write
 */
void print_matrice(double** matrice, const int sizex, const int sizey, FILE* fichier);


/**
 * @brief print a matrice_t
 *
 * @param  matrice [IN] ùatrice_t that will be printed
 * 
 * @param  fichier [IN]  file where write
 */
void print_matrice(matrice_t* matrice, FILE* fichier);


/**
 * @brief read from a file : array spin->angles and put in a spinners_t, if spin->Ngrid = 1
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  add [IN]  input file path
 */
void read_spinners(spinners_t* spin, char* add);


/**
 * @brief read from a file : array spin->angles and put in a spinners_t, if spin->Ngrid != 1
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  add [IN]  input file path
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  L [IN]  distance between tow spinners
 */
void read_spinnersall(spinners_t* spin, char* add, int nx, int ny, double L);


/**
 * @brief plot all grid in un spinner_t
 *
 * @param  spin [IN]  spinner_t with Ngrid
 */
void plotall(spinners_t* spin);

/**
 * @brief annealing simualed of a spinner grid
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 *
 * @param  TO [IN]  initial temparure
 * 
 * @param  TF [IN]  initial temparure
 * 
 * @param  lamnbda [IN]  paramter T *= lambda
 * 
 * @param  Niter [IN]  number of change at a fixed temperature
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 * 
 * @param  seed [IN]  seed for rand_d threat safe
 */
void recuit(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset, unsigned int* seed);

void plot_MinMax(char* add, double L, int nx, int ny);

/**
 * @brief performe n!/nmin!
 *
 * @param  spin [IN]  spinner_t with Ngrid
 * 
 * @return n!/nmin!
 */
unsigned long factorielle(int n, int nmin);


/**
 * @brief check if 2 grid are equal
 *
 * @param  A [IN]  spinner_t with Ngrid
 * 
 * @param  size [IN]  nx * ny
 * 
 * @param  offset1 [IN]  position of a grid in spin->angles
 *
 * @param  offset2 [IN]  position of a grid in spin->angles
 * 
 * @return returne true if equale, false else
 */
bool isequale(spinners_t* A, const int size, int offset1, int offset2);


/**
 * @brief remove all duplicate grid of spin->angles
 *
 * @param  spin [INOUT]  spinner_t with Ngrid
 */
void remove_equale(spinners_t* spin);

/**
 * @brief remove all duplicate and no metastable grid of spin->angles
 *
 * @param  spin [INOUT]  spinner_t with Ngrid
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void remove_equale_allmeta(spinners_t* spin, double* H, double* HB);


/**
 * @brief print the energy of all grid of a spinner_t
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  output file path
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void print_E(spinners_t* spin, char* add, double* H, double* HB);


/**
 * @brief plot the mean energy of all grid of a spinner_t and the standard deviation
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  track [IN]  double that will plot with the mean energy, that can be use like an abscisse
 */
void plot_E_mean(spinners_t* spin, double* H, double* HB, double track);


/**
 * @brief compute the energy and compute the histograme of all grid of a spinner_t
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  output file path
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void print_E_Histo(spinners_t* spin, char* add, double* H, double* HB);

/********************************************************************************/
/*                                                                              */
/*                                distance                                      */
/*                                                                              */
/********************************************************************************/

/**
 * @brief compute the global energy difference :  |E_i-E_j|/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return |E_i-E_j|/N
 */
double dist_EG(spinners_t* spin, int i, int j, int N, double* H, double * HB);

/**
 * @brief compute the local energy difference : \sqrt{\sum_k (E_ik-E_jk)^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k (E_ik-E_jk)^2}/N
 */
double dist_EL(spinners_t* spin, int i, int j, int N, double* H, double * HB);

/**
 * @brief compute the Hamming distance between two states : \sqrt{\sum_k (theta_ik-theta_jk)^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k (theta_ik-theta_jk)^2}/N
 */
double dist_H(spinners_t* spin, int i, int j, int N, double* H, double * HB);


/**
 * @brief compute a rotation invariant Hamming distance between tow stats : \sqrt{\sum_k \sum_neigbour (theta_ik - theta_ineigbourg[k]-theta_jk + theta_jneigbourg[k])^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k \sum_neigbour (theta_ik - theta_ineigbourg[k]-theta_jk + theta_jneigbourg[k])^2}/N
 */
double dist_HI(spinners_t* spin, int i, int j, int N, double* H, double * HB);

double fromUM(matrice_t* matrice_dist, matrice_t* matrice_ultra);

/**
 * @brief print the matrice of distance bteween two grid of all grid of a spinner_t. The matrice is "trier" and the ultrametric distance
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  output file path
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param dist [in] a fonction that performe the distance betxeen two state
 * 
 * @return the distance betwenn metric and ultrmetric distance
 */
double print_dist(spinners_t* spin, char* add, char* distchar, double*H, double *HB, double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB));

/********************************************************************************/
/*                                                                              */
/*                                clustering                                    */
/*                                                                              */
/********************************************************************************/

/**
 * @brief compute the distacne between two line of the distance matrice d = sum_l | matrice_il - matrice_jl|
 *
 * @param  B [IN]  a pointeur of the a line of the matrice
 * 
 * @param  A [IN]  a pointeur of the a line of the matrice
 *
 * @param  N [IN]  the size of a line of the matrice 
 * 
 * @return d = sum_l | matrice_il - matrice_jl|
 */
double distline(double* A, double* B, int N);

/**
 * @brief "trie" the distance matrice, by hamming distance on the ultrametric distance
 *
 * @param  matrice [INOUT]  the matrice_t of distance between two states
 * 
 * @param  ultra [INOUT]  the matrice_t of ultrametric distance between two states
 */
void tri(matrice_t* matrice, matrice_t* ultra);

/**
 * @brief performe the ultrametric distance with an average link clustering
 *
 * @param  tree [IN]  inital partition tree_t of cluster of size 1
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance between two states
 * 
 * @param  matrice_ultra [OUT]  the matrice_t of umtrametric distance between two states
 */
void cluster_fusion(tree_t* tree, matrice_t* matrice_dist, matrice_t* matrice_ultra) ;


/**
 * @brief performe the ultrametric distance matrice
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance between two states
 * 
 * @param  matrice_ultra [OUT]  the matrice_t of umtrametric distance between two states
 */
void matrice_ultra(matrice_t* matrice_dist,  matrice_t* matrice_ultra);

/********************************************************************************/
/*                                                                              */
/*                                experimente                                   */
/*                                                                              */
/********************************************************************************/

/**
 * @brief print all mestastable state of a size of grid with 0 fields
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  L [IN]  lenght between tow spinner
 * 
 * @param  add [IN]  pathfile
 */
void print_allmeta(spinners_t* spin, double L, char *add);


/**
 * @brief print all mestastable state for a B-field
 *
 * @param  spin [IN]  spin grid
 *
 * @param  L [IN]  lenght between tow spinner
 *
 * @param  add [IN]  pathfile
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 */
void print_allmeta_B(spinners_t* spin, double L, char* add, double bx, double by);

/**
 * @brief print numbert of mestastable state of a size of grid of a range of L
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  pathfile
 * 
 * @param  L0 [IN]  minimal value of L intervale
 * 
 * @param  LF [IN]  maximale value of L intervale
 * 
 * @param  Lpas [IN]  pas of discretisation of L intervale
 */
void print_allmetaofL(spinners_t* spin, char* add, double L0, double LF, double Lpas);


/**
 * @brief print numbert of mestastable state of a size of grid of a range of B
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  pathfile
 * 
 * @param  B0 [IN] B minimal
 *
 * @param  BF [IN] m maximale
 * 
 * @param  Bpas [IN] pas of range
 */
void print_allmetaofBX(spinners_t* spin, double L, char* add, double B0, double BF, double Bpas);

void print_Emin(spinners_t* spin, char* add, int Niters);

void print_Emax(spinners_t* spin, char* add, int Niters);

void print_neighbours_state_rand(spinners_t* spin, char* add, int Niters, int distance);

void print_neighbours_state_all_for(spinners_t* spin, int distancemax, int distance, int* index, int *track, FILE* fichier);

void print_neighbours_state_all(spinners_t* spin, char* add, int distance);


/**
 * @brief compute Nsimu annealing of a initial sate
 *
 * @param  spin [INOUT]  spin grid initial and all final grid
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 *
 * @param  TO [IN]  initial temparure
 * 
 * @param  TF [IN]  initial temparure
 * 
 * @param  lamnbda [IN]  paramter T *= lambda
 * 
 * @param  Niter [IN]  number of change at a fixed temperature
 * 
 * @param Nismu [IN] number of annelaing performed
 * 
 * @param p [IN] number of threats for OpenMP
 */
void recuitN(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu, int p);

/********************************************************************************/
/*                                                                              */
/*                                avalange                                      */
/*                                                                              */
/********************************************************************************/
