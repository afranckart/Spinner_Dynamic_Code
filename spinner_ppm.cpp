#include "spinner_ppm.h"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <vector>
#include <omp.h>
#include <float.h>


/********************************************************************************/
/*                                                                              */
/*                            comput Energy dip�le fct                          */
/*                                                                              */
/********************************************************************************/

double Udipole(double m1_x, double m1_y, double m2_x, double m2_y, double r1_x, double r1_y, double r2_x, double r2_y) {

	double dx = r1_x - r2_x; 
	double dy = r1_y - r2_y; 
	double r = sqrt(pow(dx, 2) + pow(dy, 2)); 
	return MU0 * 0.25 / PI * (m1_x * m2_x + m1_y * m2_y - 3. * (m1_x * dx + m1_y * dy) * (m2_x * dx + m2_y * dy) / pow(r, 2))
		/ pow(r, 3) ;
}

double UB(double mx, double my, double bx, double by) {
	return -mx * bx - my * by;
}

double UBspinner(int theta, double bx, double by) {
	double U = 0.;
	for (int i = 0; i < 6; i +=2) {
		double ang = PI3 * (theta + i);
		double Udd = UB(cos(ang), sin(ang), bx, by);
		if (i == 4) { Udd *= -1.; }
		U += Udd;
	}
	return U / HREF;
}

double Uspinner(int theta, int thetav, double L, int nv) {

	double U = 0.;
	for (int i = 0; i < 6; i += 2) {
		for (int j = 0; j < 6; j += 2) {

			double angv = PI3 * (thetav + j); 
			double ang = PI3 * (theta + i); 

			double Udd = Udipole(cos(ang), sin(ang), cos(angv), sin(angv),
				R * cos(ang), R * sin(ang), L * cos(PI3 * nv) + R * cos(angv), L * sin(PI3 * nv) + R * sin(angv));

			if ((j == 4 && i != 4) || (j != 4 && i == 4)) { Udd *= -1.; }
			U += Udd; 
		}
	}
	return U / HREF;
}


double*  H_init(double L) {

	double* H = (double*)malloc(SIZE_H * sizeof(double));

	if (H == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 6; i++) {
		for (int v = 0; v < 6; v++) {
			SETELEMENT(H, i, v, Uspinner(i, v, L, 0));
		}
	}
	return H;
}

double* H_B_init( double bx, double by) {

	double* HB = (double*)malloc(6 * sizeof(double));

	if (HB == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}


	for (int i = 0; i < 6; i++) { HB[i] = UBspinner(i, bx, by); }
	return HB;
}

void H_plot(double *H) {

	for (int i = 0; i < 6; i++) {
		printf("%f", GETELEMENT(H, i, 0));
		for (int v = 1; v < 6; v++) {
			printf("\t %f", GETELEMENT(H, i, v));
		}
		printf("\n");
	}
	printf("\n\n");
}

void H_B_plot(double* HB) {

	printf("%f", HB[0]);
	for (int i = 1; i < 6; i++) {
		printf("\t %f", HB[i]);
	}
	printf("\n\n");
}

/********************************************************************************/
/*                                                                              */
/*          Initialisation spinner lattice and lattice fct                      */
/*                                                                              */
/********************************************************************************/

void spinners_init(spinners_t *spin, double L, int nx, int ny, int Ngrid){

	spin->L = L;
	spin->nx = nx;
	spin->ny = ny;
	spin->Ngrid = Ngrid;

	spin->angles = (int*)malloc(nx * ny * Ngrid * sizeof(int));

	if (spin->angles == NULL) {
		fprintf(stderr, "spinner_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}
}

void Finalisation_simu(spinners_t* spin, double *H, double* HB) {
	free(spin->angles);
	free(H);
	free(HB);
}

int* neighbour( spinners_t* spin, int index) { //  on part de en haut � gauche et on commence � kx=ky=0

    int ky = (index - index % spin->nx) / spin->nx;
	int kx = index - ky * spin->nx; //printf("kx %d ky %d index %d\n", kx, ky, index);

	int* voisins = (int*)malloc(SIZE_NEIGHBOUR * sizeof(int));

	if (voisins == NULL) {
		fprintf(stderr, "neighbour : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}


    if (kx != 0) {
        voisins[3] = ky * spin->nx + kx - 1; // voisin de gauche
    }
	else { voisins[3]  = - 1; }

    if (ky != 0)     //peut avoir des voisins au dessus
    {
        if (ky % 2 == 1)     //peut avoir un voisin au dessu � gauche
        {
            voisins[2] = (ky - 1) * spin->nx + kx;
        }
        else if (ky % 2 == 0 && kx > 0) {
            voisins[2] = (ky - 1) * spin->nx + kx - 1;
        }
		else {
			voisins[2] = -1;
		}

        if (ky % 2 == 0)     //peut avoir un voisin au dessu � droite
        {
            voisins[1] = (ky - 1) * spin->nx + kx;
        }
        else if (ky % 2 == 1 && kx < spin->nx - 1) {
            voisins[1] = (ky - 1) * spin->nx + kx + 1;
		}
		else { voisins[1] = -1; }
    }
	else {
		voisins[2] = -1;
		voisins[1] = -1;
	}

    if (kx != spin->nx - 1) {
        voisins[0] = ky * spin->nx + kx + 1;    //voisin de droite
	}
	else { voisins[0] = -1; }

    if (ky != spin->ny - 1)     //peut avoir des voisins en dessous
    {
        if (ky % 2 == 1)     //peut avoir un voisin en dessous � gauche
        {
			voisins[4] = (ky + 1) * spin->nx + kx; 
        }
        else if (ky % 2 == 0 && kx > 0) {
            voisins[4] = (ky + 1) * spin->nx + kx - 1; 
		}
		else { voisins[4] = -1; }

        if (ky % 2 == 0)     //peut avoir un voisin en dessous � droite
        {
            voisins[5] = (ky + 1) * spin->nx + kx;
        }
        else if (ky % 2 == 1 && kx < spin->nx - 1) {
            voisins[5] = (ky + 1) * spin->nx + kx + 1;
		}
		else { voisins[5] = -1; }
	} else {// droite, au dessus � droite, au dessus � gauche, gauche, en dessous � gauche, en dessous � droite
		voisins[4] = -1;
		voisins[5] = -1;
	}
	return voisins;
}


/********************************************************************************/
/*                                                                              */
/*                           Energy of spinner lattice                          */
/*                                                                              */
/********************************************************************************/

double E_local(spinners_t* spin, int index, double* H, double* HB, int offset) {

	int* voisins = neighbour(spin, index);
	double U = 0.;
	for (int v = 0; v < SIZE_NEIGHBOUR; v++) {
		if (voisins[v] != -1) {
			U += GETELEMENT(H, (spin->angles[index + offset] - v + 6) % 6, (spin->angles[voisins[v] + offset] - v + 6) % 6); 
		}
	}
	free(voisins); 
	return U + HB[( spin->angles[index] + 6) % 6];
}

double E_total(spinners_t* spin, double*H, double* HB, int offset) {
	double U = 0.;
	for (int j = 0; j < spin->ny; j++) // colonne
	{
		for (int i = 0; i < spin->nx; i ++) // ligne
		{
			U += E_local(spin, i + j * spin->nx, H, HB, offset);
		}
	}
	return U / 2.;
}

bool metastable(spinners_t* spin, double* H, double* HB, int offset) {
	int N = spin->ny * spin->nx;
	for (int i = 0; i < N; i++) {

		spin->angles[i]++;
		double EP = E_local(spin, i, H, HB, offset);
		spin->angles[i] -= 2;
		double EM = E_local(spin, i, H, HB, offset);
		spin->angles[i]++;
		double ER = E_local(spin, i, H, HB, offset); 
		if (EM < ER || EP < ER) { return false; }
	}
	return true;
}


/********************************************************************************/
/*                                                                              */
/*                                      utilitaire                              */
/*                                                                              */
/********************************************************************************/

void change(int* A, int* B){
	int c = *A;
	*A = *B;
	*B = c;
}

FILE* openfile_out(char* add) {

	 FILE* fichier = fopen(add, "w");

    if (fichier == NULL) {
        fprintf(stderr, "openfile_out : Impossible d'ouvrir le fichier %s pour l'�criture.\n", add);
        exit(EXIT_FAILURE);
    }
	return fichier;
}

FILE* openfile_in(char* add) {

	FILE* fichier = fopen( add, "r");

	if (fichier == NULL) {
		fprintf(stderr, "openfile_in : Impossible d'ouvrir le fichier %s pour lecture.\n", add);
		exit(EXIT_FAILURE);
	}
	return fichier;
}


void print_spinners(spinners_t* spin, FILE* fichier) {

	int N = spin->nx * spin->ny;
	for(int j = 0; j < spin->Ngrid; j++){
		fprintf(fichier, "%d", spin->angles[0]);
		for (int i = 1; i < N; i++) {
			fprintf(fichier, "\t%d", (spin->angles[i + N * j] + 6) % 6);
		}
		fprintf(fichier, "\n");
	}
}

void print_matrice(double** matrice, const int sizex, const int sizey, FILE* fichier) {

	for(int i = 0; i < sizey ; i++){
		fprintf(fichier, "%f", matrice[i][0]);
		for(int j = 1; j < sizex ; j++){
			fprintf(fichier, "\t%f", matrice[i][j]);
		}
		if(i < sizey - 1){fprintf(fichier, "\n");}
	}
}


void print_matrice(matrice_t* matrice, FILE* fichier) {

	for(int i = 0; i < matrice->N ; i++){
		fprintf(fichier, "%f", matrice->line[i].col[0]);
		for(int j = 1; j < matrice->N ; j++){
			fprintf(fichier, "\t%f", matrice->line[i].col[j]);
		}
		if(i < matrice->N - 1){fprintf(fichier, "\n");}
	}
}


void read_spinners(spinners_t* spin, char* add) {

	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}
	FILE* fichier = openfile_in(add);
	int N = spin->nx * spin->ny;
	for (int i = 0; i < N; i++) {fscanf(fichier, "%d", &spin->angles[i]); }
	fclose(fichier);
}

void read_spinnersall(spinners_t* spin, char* add, int nx, int ny, double L) {

	FILE* fichier = openfile_in(add);
	int N = nx * ny;
	int value;
	std::vector<int> input(0);
	while (fscanf(fichier, "%d", &value) == 1) {input.push_back(value);}
	spin->Ngrid = (int)input.size() / nx / ny; 
	free(spin->angles);
	spin->angles = (int*)malloc(input.size() * sizeof(int));
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			spin->angles[N * i + j] = input[N * i + j];
		}
	}
	fclose(fichier);
}

void plotall(spinners_t* spin){
	int N = spin->nx * spin->ny;
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			printf("%d\t", spin->angles[i * N + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void recuit(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset, unsigned int* seed) {
	int N = spin->ny * spin->nx;
	for (double T = T0; T >= TF; T *= lambda) {
		for (int i = 0; i < Niter; i++) {
			int index = rand_r(seed) % N;
			double E0 = E_local(spin, index, H, HB, offset);
			int signe = rand_r(seed) % 2 == 0 ? 1 : -1;
			spin->angles[index + offset] += signe;
			double EF = E_local(spin, index, H, HB, offset); 
			if (E0 < EF) { 
				double temp = (double)rand_r(seed) / (double)RAND_MAX;
				if (temp > exp((E0 - EF) / T)) { spin->angles[index + offset] -= signe; }
			}
			spin->angles[index + offset] = (spin->angles[index + offset] + 6) % 6;
		}
	}
}

void plot_MinMax(char* add, double L, int nx, int ny){
	spinners_t spin;
	spinners_init(&spin,  L, nx, ny, 1);
	spinners_t spinmin;
	spinners_init(&spinmin,  L, nx, ny, 1);
	spinners_t spinmax;
	spinners_init(&spinmax,  L, nx, ny, 1);
	double* H = H_init(L);
	double* HB = H_B_init(0, 0);
	const int N = nx * ny;

	read_spinnersall(&spin, add, nx, ny, L);

	double Emax = -DBL_MAX;
	double Emin = DBL_MAX;

	for(int i = 0; i < spin.Ngrid; i++){
		double E = E_total(&spin, H, HB, N * i);
		if(E < Emin){
			Emin = E;
			for(int j = 0; j < N; j++){spinmin.angles[j] = spin.angles[N * i + j];}
		}
		if(E > Emax){
			Emax = E;
			for(int j = 0; j < N; j++){spinmax.angles[j] = spin.angles[N * i + j];}
		}
	}
	printf("Emin %f\n\n", Emin);
	for(int j = 0; j < N; j++){printf("%d\t", spinmin.angles[j]);}
	printf("\n\nEmax %f\n\n", Emax);
	for(int j = 0; j < N; j++){printf("%d\t", spinmax.angles[j]);}
	printf("\n\n");
	Finalisation_simu(&spin, H, HB);
	free(spinmax.angles);
	free(spinmin.angles);
}

unsigned long factorielle(int n, int nmin) {
	if (n == 0 || n == 1) {return 1; }
	if (nmin == n) { return nmin; }
	else {
		return n * factorielle(n - 1, nmin);
	}
}

bool isequale(spinners_t* A, const int size, int offset1, int offset2){ // passer en cuda/openMP
	for(int i = 0; i < size; i++){
		if (A->angles[offset1 + i] != A->angles[offset2 + i]){
			return false;
		}
	}
	return true;
}

void remove_equale(spinners_t* spin){
	int N = spin->nx * spin->ny;
	bool* unicity = (bool*)malloc(spin->Ngrid * sizeof(bool));
	if (unicity == NULL) {
		fprintf(stderr, "remove_equal, unicity : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < spin->Ngrid; i++){unicity[i] = true;}
	int sizeout = spin->Ngrid;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int  j = i + 1; j < spin->Ngrid; j++){ // accelerable sur GPU
				if(isequale(spin, N, i * N,  j * N) && unicity[i]){
					unicity[i] = false;
					sizeout--;
				}
			}
		}
	}

	int* out = (int*)malloc(sizeout * N * sizeof(int));
	if (out == NULL) {
		fprintf(stderr, "remove_equal, out : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	int k = 0;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int j = 0; j < N ; j++){
				out[N * k + j] = spin->angles[N * i + j];
			}
			k++;
		}
	}

	spin->angles = (int*)realloc(spin->angles, sizeout * N * sizeof(int));
	if (spin->angles == NULL) {
        fprintf(stderr, "remove_equal, new_angles :Réallocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	
	memcpy(spin->angles, out, sizeout * N * sizeof(int));
	free(out);
	free(unicity);
	printf("Ngird = %d, after remove_equale %d\n", spin->Ngrid, sizeout );
	spin->Ngrid = sizeout;
}

void remove_equale_allmeta(spinners_t* spin, double* H, double* HB){
	int N = spin->nx * spin->ny;
	bool* unicity = (bool*)malloc(spin->Ngrid * sizeof(bool));
	if (unicity == NULL) {
		fprintf(stderr, "remove_equal_allmeta, unicity : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < spin->Ngrid; i++){unicity[i] = true;}
	int sizeout = spin->Ngrid;
	for(int i =0; i < spin->Ngrid; i++){
		if(!metastable(spin, H, HB, N*i)){
			unicity[i] = false;
			sizeout--;
		};
	}
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int  j = i + 1; j < spin->Ngrid; j++){ // accelerable sur GPU
				if(isequale(spin, N, i * N,  j * N) && unicity[i]){
					unicity[i] = false;
					sizeout--;
				}
			}
		}
	}

	int* out = (int*)malloc(sizeout * N * sizeof(int));
	if (out == NULL) {
		fprintf(stderr, "remove_equal_allmeta, out : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	int k = 0;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int j = 0; j < N ; j++){
				out[N * k + j] = spin->angles[N * i + j];
			}
			k++;
		}
	}

	spin->angles = (int*)realloc(spin->angles, sizeout * N * sizeof(int));
	if (spin->angles == NULL) {
        fprintf(stderr, "remove_equal_allmeta, new_angles :Réallocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	
	memcpy(spin->angles, out, sizeout * N * sizeof(int));
	free(out);
	free(unicity);
	printf("Ngird = %d, after remove_equale_allmeta %d\n", spin->Ngrid, sizeout );
	spin->Ngrid = sizeout;
}

void print_E(spinners_t* spin, char* add, double* H, double* HB){
	FILE* fichier = openfile_out(add);
	int N = spin->nx * spin->ny;
	double* E = (double*)malloc(spin->Ngrid * sizeof(double));
	if (E == NULL) {
        fprintf(stderr, "print_E, E : allocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	for(int i = 0; i < spin->Ngrid; i++){ E[i] = E_total(spin, H, HB, N * i); }
	for(int i = 0; i < spin->Ngrid - 1; i++){ fprintf(fichier, "%f\n", E[i]); }
	fprintf(fichier, "%f", E[spin->Ngrid - 1]);
	free(E);
	fclose(fichier);
}

void plot_E_mean(spinners_t* spin, double* H, double* HB, double track){
	int N = spin->nx * spin->ny;
	double* E = (double*)malloc(spin->Ngrid * sizeof(double));
	if (E == NULL) {
        fprintf(stderr, "print_E, E : allocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	double Emin = DBL_MAX;
	double Emax = -DBL_MAX;
	for(int i = 0; i < spin->Ngrid; i++){ 
		E[i] = E_total(spin, H, HB, N * i);
		if(E[i] > Emax) Emax = E[i];
		if(E[i] < Emin) Emin = E[i];
	}
	double moyenne = 0.0;
    for (int i = 0; i < spin->Ngrid; i++) { moyenne += E[i]; }
    moyenne /= spin->Ngrid;
    double ecart_type = 0.0;
    for (int i = 0; i < spin->Ngrid; i++) { ecart_type += pow(E[i] - moyenne, 2);}
    ecart_type = sqrt(ecart_type / spin->Ngrid);
	printf("%f %f %f %f %f\n", track , moyenne, ecart_type, Emin, Emax);
	free(E);
}


void print_E_Histo(spinners_t* spin, char* add, double* H, double* HB){
	FILE* fichier = openfile_out(add);
	int N = spin->nx * spin->ny;
	std::map<double, int> histogram;
	for(int i = 0; i < spin->Ngrid; i++){ histogram[E_total(spin, H, HB, N * i)]++; }
	for (const auto& entry : histogram){ fprintf(fichier, "%f %d\n", entry.first, entry.second); }
	fclose(fichier);
}

/********************************************************************************/
/*                                                                              */
/*                                distance                                      */
/*                                                                              */
/********************************************************************************/

double dist_EG(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	return abs(E_total(spin, H, HB, N * i) - E_total(spin, H, HB, N * j)) / (double)N;
}

double dist_EL(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
    double E = E_local(spin, k, H, HB, i * N) - E_local(spin, k, H, HB, j * N);
    d += E *E ;
  }
	return sqrt(d / (double)N);
}

double dist_H(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
   int delta = (spin->angles[i * N + k] - spin->angles[j * N + k]) * (spin->angles[i * N + k] - spin->angles[j * N + k]);
    switch( delta ){
      case 25:
        d += 1.;
        break;
      case 16 :
        d += 4.;
        break;
      default:
        d += delta;
        break;
    }
  }
	return sqrt(d / (double)N) ;
}

double dist_HI(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
		double dd = 0;
		int* voisins = neighbour(spin, k);
		for(int l = 0; l < SIZE_NEIGHBOUR; l++){
			if(voisins[l] != -1){
        
       int delta = spin->angles[i * N + k] - spin->angles[i * N + voisins[l]] - spin->angles[j * N + k] + spin->angles[j * N + voisins[l]];
       delta *= delta;
       switch( delta ){
          case 25:
          dd += 1.;
          break;
        case 16 :
          dd += 4.;
          break;
        default:
          dd += delta;
          break;
        }
      }
		}
		d += dd;
		free(voisins);
	}
	return sqrt(d / (double)N) ;
}

double fromUM(matrice_t* matrice_dist, matrice_t* matrice_ultra){
	double d = 0;
	for(int i = 0; i < matrice_ultra->N ; i++){
		for(int j = i + 1; j < matrice_ultra->N ; j++){
			d += (matrice_dist->line[i].col[j] - matrice_ultra->line[i].col[j]) * (matrice_dist->line[i].col[j] - matrice_ultra->line[i].col[j]);
		}
	}
	return sqrt(d / (double)matrice_ultra->N);
}

double print_dist(spinners_t* spin, char* add, char* distchar, double*H, double *HB, double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB)){
	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L"); 
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", spin->L);
	strcat(path, "_dist"); 
	strcat(path, distchar); 
	strcat(path, ".txt");


	char path2[STRING_MAX];
	strcpy(path2, add);
	strcat(path2, "_L"); 
	snprintf(path2 + strlen(path2), sizeof(path2) - strlen(path2), "%f", spin->L);
	strcat(path2, "_dist"); 
	strcat(path2, distchar); 
	strcat(path2, "_ultra.txt");

	FILE* fichier = openfile_out(path);
	FILE* fichier_ultra = openfile_out(path2);

	matrice_t matrice;
	matrice_t ultra;
	matrice.N = spin->Ngrid;
	ultra.N = spin->Ngrid;
	matrice.line = (matrice_line_t*)malloc( matrice.N * sizeof(matrice_line_t));
	ultra.line = (matrice_line_t*)malloc( ultra.N * sizeof(matrice_line_t));
	if (matrice.line == NULL || ultra.line == NULL) {
		fprintf(stderr, "print_matrice, matrice.line and ultra.line : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < matrice.N; i++){
		matrice.line[i].col = (double*)malloc( matrice.N * sizeof(double));
		ultra.line[i].col = (double*)malloc( ultra.N * sizeof(double));
		if (matrice.line[i].col == NULL || ultra.line[i].col == NULL) {
			fprintf(stderr, "print_matrice, matrice.line[%d].col and ultra.line[%d].col : Allocation de memoire echouee.\n", i, i);
			exit(EXIT_FAILURE);
		}
		matrice.line[i].pos = i;
		ultra.line[i].pos = i;
	}

	const int N = spin->nx * spin->ny;
	for(int i = 0; i < spin->Ngrid ; i ++){
		for(int j = i; j < spin->Ngrid ; j ++){
			double d = dist(spin, i, j, N, H, HB);
			matrice.line[i].col[j] = d;
			matrice.line[j].col[i] = d;
		}
	}
	
	matrice_ultra(&matrice, &ultra);
	tri(&matrice, &ultra);
	print_matrice(&matrice, fichier);
	print_matrice(&ultra, fichier_ultra);
	double deltaUM =fromUM(&matrice, &ultra);
	fclose(fichier);
	fclose(fichier_ultra);
	for(int i = 0; i < matrice.N; i++){free(matrice.line[i].col);}
	free(matrice.line);
	for(int i = 0; i < ultra.N; i++){free(ultra.line[i].col);}
	free(ultra.line);
	return deltaUM;
}

/********************************************************************************/
/*                                                                              */
/*                                clustering                                    */
/*                                                                              */
/********************************************************************************/

double distline(double* A, double* B, int N){
	double d = 0;
	for(int i = 0; i < N; i++){ d += fabs(A[i] - B[i]);}
	return d;
}

void tri(matrice_t* matrice, matrice_t* ultra){
	
	for (int i = 0; i < ultra->N - 1; i++) {
        double dmin = distline(ultra->line[i].col, ultra->line[i + 1].col, ultra->N);
        for (int j = i + 2; j < ultra->N; j++) {
			double d = distline(ultra->line[i].col, ultra->line[j].col, ultra->N);
            if (d < dmin ) {
				double* change = ultra->line[j].col;
				ultra->line[j].col = ultra->line[i + 1].col;
				ultra->line[i + 1].col = change; 

				dmin = d;

				int temp = ultra->line[i + 1].pos;
				ultra->line[i + 1].pos = ultra->line[j].pos;
				ultra->line[j].pos = temp;
            }
        }
    }

	for (int i = 0; i < ultra->N - 1; i++) {
        double dmin = ultra->line[i].col[0];
        for (int j = i + 1; j < ultra->N; j++) {
			double d = ultra->line[j].col[0];
            if (d < dmin ) {

				double* change = ultra->line[j].col;
				ultra->line[j].col = ultra->line[i].col;
				ultra->line[i].col = change; 

				dmin = d;

				int temp = ultra->line[i].pos;
				ultra->line[i].pos = ultra->line[j].pos;
				ultra->line[j].pos = temp;
            }
        }
    }
	
	double* change1 = (double*)malloc(matrice->N * sizeof(double));
	double* change2 = (double*)malloc(matrice->N * sizeof(double));
	double** change3 = (double**)malloc(matrice->N * sizeof(double*));
	if (change1 == NULL || change2 == NULL || change3 == NULL) {
		fprintf(stderr, "tri, change : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < matrice->N; i++){change3[i] = matrice->line[i].col;}
	for(int i = 0; i < matrice->N; i++){matrice->line[i].col = change3[ultra->line[i].pos];}
	free(change3);

	for(int i = 0; i < ultra->N; i++){
		memcpy(change1, matrice->line[i].col, matrice->N * sizeof(double));
		memcpy(change2, ultra->line[i].col, ultra->N * sizeof(double));
		for(int j = 0; j < ultra->N; j++){ultra->line[i].col[j] = change2[ultra->line[j].pos];}
		for(int j = 0; j < matrice->N; j++){matrice->line[i].col[j] = change1[ultra->line[j].pos];}
	}
	free(change1);
	free(change2);
}


void cluster_fusion(tree_t* tree, matrice_t* matrice_dist, matrice_t* matrice_ultra) {
	
	int Ncluster = tree->N;
	int* newClusterPos = (int*)malloc(tree->N * sizeof(int));
    if (newClusterPos == NULL) {
        fprintf(stderr, "cluster_fusion, newClusterPos : Allocation de memoire echouee.\n");
        exit(EXIT_FAILURE);
    }

    while (Ncluster > 1) {
        
        int minI = -1, minJ = -1;
        double minSimilarite = DBL_MAX;
		double similarity = 0.;

        for (int i = 0; i < tree->N; i++) {
            if (tree->cluster[i].notmerged) {
                for (int j = i + 1; j < tree->N; j++) {
                    if (tree->cluster[j].notmerged) {
                        similarity = 0.;
						for(int k = 0; k < tree->cluster[i].size;  k++){
							for(int l = 0; l < tree->cluster[j].size;  l++){
								similarity += matrice_dist->line[tree->cluster[i].pos[k]].col[tree->cluster[j].pos[l]];
							}
						}
						similarity /= tree->cluster[i].size * tree->cluster[j].size; 
                        if ( similarity < minSimilarite) {
                            minSimilarite = similarity;
                            minI = i;
                            minJ = j;
                        }
                    }
                }
            }
        }

		for (int k = 0; k < tree->cluster[minI].size; k++) {
            for (int l = 0; l < tree->cluster[minJ].size; l++) {
				matrice_ultra->line[tree->cluster[minI].pos[k]].col[tree->cluster[minJ].pos[l]] = minSimilarite;
				matrice_ultra->line[tree->cluster[minJ].pos[l]].col[tree->cluster[minI].pos[k]] = minSimilarite;
        	}
        }

        int newClusterSize = tree->cluster[minI].size + tree->cluster[minJ].size;

        memcpy(newClusterPos, tree->cluster[minI].pos, tree->cluster[minI].size * sizeof(int));
        memcpy(newClusterPos + tree->cluster[minI].size, tree->cluster[minJ].pos, tree->cluster[minJ].size * sizeof(int));

        tree->cluster[minI].size = newClusterSize;
        tree->cluster[minI].pos = (int*)realloc(tree->cluster[minI].pos, newClusterSize * sizeof(int));
		if (tree->cluster[minI].pos == NULL) {
            fprintf(stderr, "cluster_fusion, tree->cluster[%d].pos : Allocation de memoire echouee.\n", minI);
            exit(EXIT_FAILURE);
        }
		memcpy(tree->cluster[minI].pos, newClusterPos, newClusterSize * sizeof(int));

        tree->cluster[minJ].notmerged = false;

		Ncluster--;
    }
	free(newClusterPos);
}

void matrice_ultra(matrice_t* matrice_dist, matrice_t* matrice_ultra){
	tree_t tree;
	tree.N = matrice_dist->N;
	tree.cluster = (cluster_t*)malloc( tree.N * sizeof(cluster_t));
	if (tree.cluster == NULL) {
		fprintf(stderr, "matrice_ultra, tree.cluster : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < tree.N; i++){
		tree.cluster[i].pos = (int*)malloc(sizeof(int));
		if (tree.cluster[i].pos == NULL) {
			fprintf(stderr, "matrice_ultra, tree.cluster[%d].pos : Allocation de memoire echouee.\n", i);
			exit(EXIT_FAILURE);
		}
		tree.cluster[i].size = 1;
		tree.cluster[i].pos[0] = i;
		tree.cluster[i].notmerged = true;
	}

	for(int i = 0; i < matrice_ultra->N; i++){matrice_ultra->line[i].col[i] = 0;}

	cluster_fusion(&tree, matrice_dist, matrice_ultra);

	for(int i = 0; i < tree.N; i++){free(tree.cluster[i].pos);}
	free(tree.cluster);
}


/********************************************************************************/
/*                                                                              */
/*                                experimente                                   */
/*                                                                              */
/********************************************************************************/

void print_all4x4(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;
	
  #pragma omp parallel
  {
  	spinners spin;
	spinners_init(&spin, L, 4, 4, 1);
	#pragma omp for collapse(6)
	for(int i0 = 0 ; i0 < 6; i0++){
		for(int i1 = 0 ; i1 < 6; i1++){
			for(int i2 = 0 ; i2 < 6; i2++){
				for(int i3 = 0 ; i3 < 6; i3++){
					for(int i4 = 0 ; i4 < 6; i4++){
						for(int i5 = 0 ; i5 < 6; i5++){
							spin.angles[0] = i0;
							spin.angles[1] = i1;
							spin.angles[2] = i2;
							spin.angles[3] = i3;
							spin.angles[4] = i4;
							spin.angles[5] = i5;
							for(int i6 = 0 ; i6 < 6; i6++){
								spin.angles[6] = i6;
								for(int i7 = 0 ; i7 < 6; i7++){
									spin.angles[7] = i7;
									for(int i8 = 0 ; i8 < 6; i8++){
										spin.angles[8] = i8;
										for(int i9 = 0 ; i9 < 6; i9++){
											spin.angles[9] = i9;
											for(int i10 = 0 ; i10 < 6; i10++){
												spin.angles[10] = i10;
												for(int i11 = 0 ; i11 < 6; i11++){
													spin.angles[11] = i11;
													for(int i12 = 0 ; i12 < 6; i12++){
														spin.angles[12] = i12;
														for(int i13 = 0 ; i13 < 6; i13++){
															spin.angles[13] = i13;
															for(int i14 = 0 ; i14 < 6; i14++){
																spin.angles[14] = i14;
																for(int i15 = 0 ; i15 < 6; i15++){
																	spin.angles[15] = i15;
																	if(metastable(&spin, H, HB, 0)) {
																		print_spinners(&spin, fichier);
																		n++;
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
  free(spin.angles);
  }
	printf("%f\t%d\n",L, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_all3x3(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;
	
  #pragma omp parallel
  {
  	spinners spin;
	spinners_init(&spin, L, 3, 3, 1);
	#pragma omp for collapse(6)
	for(int i0 = 0 ; i0 < 6; i0++){
		for(int i1 = 0 ; i1 < 6; i1++){
			for(int i2 = 0 ; i2 < 6; i2++){
				for(int i3 = 0 ; i3 < 6; i3++){
					for(int i4 = 0 ; i4 < 6; i4++){
						for(int i5 = 0 ; i5 < 6; i5++){
							spin.angles[0] = i0;
							spin.angles[1] = i1;
							spin.angles[2] = i2;
							spin.angles[3] = i3;
							spin.angles[4] = i4;
							spin.angles[5] = i5;
							for(int i6 = 0 ; i6 < 6; i6++){
								spin.angles[6] = i6;
								for(int i7 = 0 ; i7 < 6; i7++){
									spin.angles[7] = i7;
									for(int i8 = 0 ; i8 < 6; i8++){
										spin.angles[8] = i8;
										if(metastable(&spin, H, HB, 0)) {
											print_spinners(&spin, fichier);
											n++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
  free(spin.angles);
  }
	printf("%f\t%d\n",L, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_all2x2(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;
	
 
  	spinners spin;
	spinners_init(&spin, L, 2, 2, 1);
	for(int i0 = 0 ; i0 < 6; i0++){
		spin.angles[0] = i0;
		for(int i1 = 0 ; i1 < 6; i1++){
			spin.angles[1] = i1;
			for(int i2 = 0 ; i2 < 6; i2++){
				spin.angles[2] = i2;
				for(int i3 = 0 ; i3 < 6; i3++){
					spin.angles[3] = i3;
					if(metastable(&spin, H, HB, 0)) {
						print_spinners(&spin, fichier);
						n++;
					}
				}
			}
		}
	}
	free(spin.angles);
  
	printf("%f\t%d\n",L, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_Emin( spinners_t* spin,  char* add, int Niters)
{
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	FILE* fichier = openfile_out(add);
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);
	int N = spin->nx * spin->ny;
	double Emin = DBL_MAX;
	int* spinmin = (int*)malloc(N * sizeof(int));

	unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

	for (int i = 0; i < Niters; i++) {
		for (int j = 0; j < N; j++) { spin->angles[j] = rand() % 6; }
		recuit(spin, H, HB, 1, 0.001, 0.95, 100 * N, 0, &seed);
		double E = E_total(spin, H, HB, 0);
		if (E < Emin && metastable(spin, H, HB, 0)) { 
			Emin = E;
			for (int l = 0; l < N; l++) { spinmin[l] = spin->angles[l]; }
		}

	}
	for (int l = 0; l < N; l++) { spin->angles[l] = spinmin[l]; }
	printf("nx %d ny %d EF = %f metastable : %d\n",spin->nx, spin->ny, Emin, metastable(spin, H, HB, 0));
	print_spinners(spin, fichier);
	free(H);
	free(HB);
	free(spinmin);
	fclose(fichier);
}

void print_Emax( spinners_t* spin,  char* add, int Niters){
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	FILE* fichier = openfile_out(add);
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);
	int N = spin->nx * spin->ny;
	double Emax = -DBL_MAX;
	int* spinmin = (int*)malloc(N * sizeof(int));

	unsigned int seed = (unsigned int)time(NULL);

	for (int i = 0; i < Niters; i++) {
		for (int j = 0; j < N; j++) { spin->angles[j] = rand() % 6; }
		recuit(spin, H, HB, 0.0001, 0.00001, 0.5, 5 * N, 0, &seed);
		double E = E_total(spin, H, HB, 0);
		if (E > Emax && metastable(spin, H, HB, 0)) { 
			Emax = E;
			for (int l = 0; l < N; l++) { spinmin[l] = spin->angles[l]; }
		}
	}
	for (int l = 0; l < N; l++) { spin->angles[l] = spinmin[l]; }
	printf("nx = %d ny = %d, EF = %f metastable : %d\n",spin->nx, spin->ny, Emax, metastable(spin, H, HB, 0));
	print_spinners(spin, fichier);
	free(H);
	free(HB);
	free(spinmin);
	fclose(fichier);
}

void print_neighbours_state_rand(spinners_t* spin, char* add, int Niters, int distance){
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}
	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", spin->L);
	strcat(path, "_distance");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%d", distance);
	strcat(path, "_iter");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%d", Niters);
	strcat(path, ".txt");

	FILE* fichier = openfile_out(path);
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);
	int N = spin->nx * spin->ny;

	srand((unsigned int)time(NULL));
	int* indexes = (int*)malloc(distance * sizeof(int));
	int* changes = (int*)malloc(distance * sizeof(int));

	int track = 0;

	for (int i = 0; i < Niters; i++) {
		bool go = true;
		for (int j = 0; j < distance; j++) { indexes[j] = rand() % N;}
		for (int j = 0; j < distance; j++) { 
			for (int l = j + 1; l < distance; l++) { 
				if ( indexes[j] == indexes[l]) {
					go = false;
					break;
				};
			}
		}
		if (go) {
			for (int j = 0; j < distance; j++) { changes[j] = rand() % 2 == 0 ? 1 : -1; }
			for (int j = 0; j < distance; j++) { spin->angles[indexes[j]] += changes[j]; }
			if (metastable(spin, H, HB, 0)) { 
				print_spinners(spin, fichier);
				track++;
			}
			for (int j = 0; j < distance; j++) { spin->angles[indexes[j]] -= changes[j]; }
		}
	}
	printf("nombre d'etats trouver = %d\n", track);
	free(H);
	free(indexes);
	free(changes);
	free(HB);
	fclose(fichier);
}


/*modifier car des �tats sont compter plusieur fois*/
void print_neighbours_state_all_for(spinners_t* spin, int distancemax, int distance, int* index, int* track, FILE* fichier){
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);
	int N = spin->nx * spin->ny;

	if (distance == 1) {
		for (int i = 0; i < N; i++) {
			bool go = true;
			for (int j = 0; j < distancemax; j++) {
				if (i == index[j]) {
					go = false;
					break;
				}
			}
			if (go) {
				spin->angles[i] += 1;
				if (metastable(spin, H, HB, 0)) {
					print_spinners(spin, fichier);
					*track += 1; 
				}
				spin->angles[i] -= 2;
				if (metastable(spin, H, HB, 0)) {
					print_spinners(spin, fichier);
					*track += 1;
				}
				spin->angles[i] += 1;
			}
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			index[distance - 1] = i;
			spin->angles[i] += 1;
			print_neighbours_state_all_for(spin, distancemax, distance - 1, index, track, fichier);
			spin->angles[i] -= 2;
			print_neighbours_state_all_for(spin, distancemax,  distance - 1, index, track, fichier);
			spin->angles[i] += 1;
		}
	}
	free(H);
	free(HB);
}

void print_neighbours_state_all(spinners_t* spin, char* add, int distance)
{

	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}
	
	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", spin->L);
	strcat(path, "_distance");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%d", distance);
	strcat(path, "_all.txt");
	int* index = (int*)malloc(distance * sizeof(int));
	index[0] = -1;
	FILE* fichier = openfile_out(path);
	int N = spin->nx * spin->ny;
	int track = 0;
	print_neighbours_state_all_for(spin, distance, distance, index, &track, fichier);
	printf("distance %d , nombre d'etats trouver = %d sur %lu\n", distance, track, factorielle(N, N - distance));
	free(index);
	fclose(fichier);
}

void recuitN(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu, int p){

	int N = spin->nx * spin->ny;
	int* angles0 = (int*)malloc(N * sizeof(int));
	if (angles0  == NULL) {
        fprintf(stderr, "recuitN, angles0 failled to malloc\n");
        exit(EXIT_FAILURE);
	}

	memcpy(angles0, spin->angles, N * sizeof(int));
	spin->angles = (int*)realloc(spin->angles, Nsimu * N * sizeof(int));
	if (spin->angles == NULL) {
    	fprintf(stderr, "recuitN, échec de l'allocation mémoire.\n");
    	exit(EXIT_FAILURE);	
	}

	spin->Ngrid = Nsimu;
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			spin->angles[i * N + j] = angles0[j];
		}
	}
	free(angles0);

	#pragma omp parallel num_threads(p)
	{	
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < spin->Ngrid ; i++){
			recuit(spin, H, HB, T0, TF, lambda, Niter, N * i, &seed);
		}
	}

	remove_equale_allmeta(spin, H, HB);
}



void recuit_GPU(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu, int p){

	int N = spin->nx * spin->ny;
	int* angles0 = (int*)malloc(N * sizeof(int));
	if (angles0  == NULL) {
        fprintf(stderr, "recuitN, angles0 failled to malloc\n");
        exit(EXIT_FAILURE);
	}

	memcpy(angles0, spin->angles, N * sizeof(int));
	spin->angles = (int*)realloc(spin->angles, Nsimu * N * sizeof(int));
	if (spin->angles == NULL) {
    	fprintf(stderr, "recuitN, échec de l'allocation mémoire.\n");
    	exit(EXIT_FAILURE);	
	}

	spin->Ngrid = Nsimu;
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			spin->angles[i * N + j] = angles0[j];
		}
	}
	free(angles0);

	#pragma omp parallel num_threads(p)
	{	
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < spin->Ngrid ; i++){
			for (double T = T0; T >= TF; T *= lambda) {
				for (int j = 0; j < Niter; j++) {
					int index = rand_r(&seed) % N;
					double E0 = E_local(spin, index, H, HB, i * N);
					int signe = rand_r(&seed) % 2 == 0 ? 1 : -1;
					spin->angles[index + i * N] += signe;
					double EF = E_local(spin, index, H, HB, i * N); 
					if (E0 < EF) { 
						double temp = (double)rand_r(&seed) / (double)RAND_MAX;
						if (temp > exp((E0 - EF) / T)) { spin->angles[index + i * N] -= signe; }
					}
				spin->angles[index + i * N] = (spin->angles[index + i * N] + 6) % 6;
				}
			}
		}
	}

	remove_equale_allmeta(spin, H, HB);
}


void print_Emin( double L,  int nx, int ny, char* add, int Niters, int p)
{
	spinners_t spin;
	spinners_init(&spin, L, nx, ny, Niters);
	spinners_t spinmin;
	spinners_init(&spinmin, L, nx, ny, 1);

	double* H = H_init(spin.L);
	double* HB = H_B_init(0, 0);
	int N = spin.nx * spin.ny;
	double Emin = DBL_MAX;

	#pragma omp parallel num_threads(p)
	{	
		unsigned int seed = (unsigned int)(time(NULL) + omp_get_thread_num());
		#pragma omp for
		for(int i = 0; i < N ; i++){
			for(int j = 0; j < N ; j++){spin.angles[i * N + j] = rand_r(&seed) % 6;}
		}

		#pragma omp for
		for(int i = 0; i < spin.Ngrid ; i++){
			recuit(&spin, H, HB, 1, 0.001, 0.95, 100 * N, N * i, &seed);
		}
	}
	
	for (int i = 0; i < spin.Ngrid; i++) {
		double E = E_total(&spin, H, HB, i * N);
		if (E < Emin && metastable(&spin, H, HB, i * N)) { 
			Emin = E;
			for (int l = 0; l < N; l++) { spinmin.angles[l] = spin.angles[l + i * N]; }
		}

	}
	
	printf("nx %d ny %d EF = %f metastable : %d\n",spin.nx, spin.ny, Emin, metastable(&spinmin, H, HB, 0));
	if(metastable(&spinmin, H, HB, 0)){
		FILE* fichier = openfile_out(add);
		print_spinners(&spinmin, fichier);
		fclose(fichier);
	}
	Finalisation_simu(&spin, H, HB);
	free(spinmin.angles);
}

/********************************************************************************/
/*                                                                              */
/*                                avalange                                      */
/*                                                                              */
/********************************************************************************/