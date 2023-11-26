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
 * @brief print a matrice sizex * sizey
 *
 * @param  matrice [IN] array to print 
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  sizex [IN]  horizontale size of matrice
 * 
 * @param  sizey [IN]  verticale size of matrice
 * 
 * @param  fichier [IN]  file where write
 */
void print_matrice(double* matrice, const int sizex, const int sizey, FILE* fichier);


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
 */
void recuit(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset);

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

/********************************************************************************/
/*                                                                              */
/*                                distance                                      */
/*                                                                              */
/********************************************************************************/

void print_dist(spinners_t* spin, char* add, char* distchar, double (*dist)(spinners_t*, int i, int j, int N));

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

void recuitN(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu);
