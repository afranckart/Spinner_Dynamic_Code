#include "spinner_ppm.h"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <time.h>
//#include <map>
#include <vector>

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

	spin->angles = (int*)malloc(nx * ny * sizeof(int));

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
		fprintf(stderr, "openfile_in : Impossible d'ouvrir le fichier %s pour l'�criture.\n", add);
		exit(EXIT_FAILURE);
	}
	return fichier;
}


void print_spinners(spinners_t* spin, FILE* fichier) {

	int N = spin->nx * spin->ny;
	for(int j = 0; j < spin->Ngrid; j++){
		fprintf(fichier, "%d", spin->angles[0]);
		for (int i = 1; i < N; i++) {
			fprintf(fichier, "\t%d", (spin->angles[i] + 6) % 6);
		}
	}
	fprintf(fichier, "\n");
}

void print_matrice(double* matrice, const int sizex, const int sizey, FILE* fichier) {

	for(int i = 0; i < sizey ; i++){
		fprintf(fichier, "%f", matrice[i * sizex]);
		for(int j = 1; j < sizex ; j++){
			fprintf(fichier, "\t%f", matrice[i * sizex + j]);
		}
		if(i < sizey - 1){fprintf(fichier, "\n");}
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

void recuit(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset) {
	int N = spin->ny * spin->nx;
	for (double T = T0; T >= TF; T *= lambda) {
		for (int i = 0; i < Niter; i++) {
			int index = rand() % N;
			double E0 = E_local(spin, index, H, HB, offset);
			int signe = rand() % 2 == 0 ? 1 : -1;
			spin->angles[index] += signe;
			double EF = E_local(spin, index, H, HB, offset); 
			if (E0 < EF) { 
				double temp = (double)rand() / (double)RAND_MAX;
				if (temp > exp((E0 - EF) / T)) { spin->angles[index] -= signe; }
			}
			spin->angles[index] = (spin->angles[index] + 6) % 6;
		}
	}
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

/********************************************************************************/
/*                                                                              */
/*                                distance                                      */
/*                                                                              */
/********************************************************************************/


void print_dist(spinners_t* spin, char* add, char* distchar, double (*dist)(spinners_t*, int i, int j, int N)){
	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L"); 
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", spin->L); 
	strcat(path, "_dist.txt"); 
	strcat(path, distchar); 
	strcat(path, ".txt"); 

	FILE* fichier = openfile_out(path);

	const int N = spin->nx * spin->ny;
	double* matrice = (double*)malloc( spin->Ngrid * spin->Ngrid * sizeof(double)); 
	for(int i = 0; i < spin->Ngrid ; i ++){
		for(int j = i; j < spin->Ngrid ; j ++){
			matrice[i * N + j] = dist(spin, i, j, N);
			matrice[j * N + i] = matrice[i * N + j];
		}
	}
	print_matrice(matrice, spin->Ngrid, spin->Ngrid, fichier);
	fclose(fichier);
	free(matrice);
}

/********************************************************************************/
/*                                                                              */
/*                                experimente                                   */
/*                                                                              */
/********************************************************************************/

void print_allmeta(spinners_t* spin, double L, char* add) {
	
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L"); 
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", L); 
	strcat(path, "_allmeta.txt"); 

	FILE* fichier = openfile_out(path);

	int N = spin->nx * spin->ny;

	for (int l = 0; l < N; l++) { spin->angles[l] = 0; }
	double* H = H_init(L);
	double* HB = H_B_init(0, 0);
	while (true) {

		// Incr�menter la combinaison
		int j = N - 1;
		while (j >= 0 && spin->angles[j] == 5) {
			spin->angles[j] = 0;
			j--;
		}

		// Si tous les chiffres sont 5, nous avons termin�.
		if (j < 0) {
			break;
		}

		// Incr�menter le dernier chiffre non 5
		spin->angles[j]++;

		if (metastable(spin, H, HB, 0)) { print_spinners(spin, fichier); }
	}
	free(H);
	free(HB);
	fclose(fichier);
}

void print_allmeta_B(spinners_t* spin, double L, char* add, double bx, double by) {

	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d\n", spin->Ngrid );
	}

	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", L);
	strcat(path, "_BX");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", bx);
	strcat(path, "_BY");
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", by);
	strcat(path, "_allmeta.txt");

	FILE* fichier = openfile_out(path);

	int N = spin->nx * spin->ny;

	for (int l = 0; l < N; l++) { spin->angles[l] = 0; }
	double* H = H_init(L);
	double* HB = H_B_init(bx, by);
	while (true) {

		// Incr�menter la combinaison
		int j = N - 1;
		while (j >= 0 && spin->angles[j] == 5) {
			spin->angles[j] = 0;
			j--;
		}

		// Si tous les chiffres sont 5, nous avons termin�.
		if (j < 0) {
			break;
		}

		// Incr�menter le dernier chiffre non 5
		spin->angles[j]++; 
		
		if (metastable(spin, H, HB, 0)) { print_spinners(spin, fichier); }
	}
	free(H);
	free(HB);
	fclose(fichier);
}

void print_allmetaofL(spinners_t* spin, char* add, double L0, double LF, double Lpas)
{
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}
	FILE* fichier = openfile_out(add);

	int N = spin->nx * spin->ny;

	for (int l = 0; l < N; l++) { spin->angles[l] = 0; }
	double* HB = H_B_init(0, 0);
	for (double l = L0; l <= LF; l += Lpas) {
		double* H = H_init(l);
		int Nstable = 0;
		while (true) {

			// Incr�menter la combinaison
			int j = N - 1;
			while (j >= 0 && spin->angles[j] == 5) {
				spin->angles[j] = 0;
				j--;
			}

			// Si tous les chiffres sont 5, nous avons termin�.
			if (j < 0) {
				break;
			}

			// Incr�menter le dernier chiffre non 5
			spin->angles[j]++;

			if (metastable(spin, H, HB, 0)) { Nstable++; }
		}
		free(H);
		fprintf(fichier, "%f\t%d ", l, Nstable);
		if (l + Lpas < LF) { fprintf(fichier, "\n"); }
	}
	free(HB);
	fclose(fichier);
}

void print_allmetaofBX(spinners_t* spin, double L, char* add, double B0, double BF, double Bpas)
{
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	FILE* fichier = openfile_out(add);

	int N = spin->nx * spin->ny;

	for (int l = 0; l < N; l++) { spin->angles[l] = 0; }
	for (double b = B0; b <= BF; b += Bpas) {
		double* H = H_init(L);
		double* HB = H_B_init( b, 0);
		int Nstable = 0;
		while (true) {

			// Incr�menter la combinaison
			int j = N - 1;
			while (j >= 0 && spin->angles[j] == 5) {
				spin->angles[j] = 0;
				j--;
			}

			// Si tous les chiffres sont 5, nous avons termin�.
			if (j < 0) {
				break;
			}

			// Incr�menter le dernier chiffre non 5
			spin->angles[j]++;

			if (metastable(spin, H, HB, 0)) { Nstable++; }
		}
		free(H);
		free(HB);
		fprintf(fichier, "%f\t%d ", b, Nstable);
		if (b + Bpas < BF) { fprintf(fichier, "\n"); }
	}

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
	double Emin = 1000000 * N;
	int* spinmin = (int*)malloc(N * sizeof(int));

	srand((unsigned int)time(NULL));

	for (int i = 0; i < Niters; i++) {
		for (int j = 0; j < N; j++) { spin->angles[j] = rand() % 6; }
		recuit(spin, H, HB, 10, 0.0001, 0.95, 5 * N, 0);
		double E = E_total(spin, H, HB, 0);
		if (E < Emin && metastable(spin, H, HB, 0)) { 
			Emin = E;
			for (int l = 0; l < N; l++) { spinmin[l] = spin->angles[l]; }
		}
	}
	for (int l = 0; l < N; l++) { spin->angles[l] = spinmin[l]; }
	printf("EF = %f metastable : %d\n", Emin, metastable(spin, H, HB, 0));
	print_spinners(spin, fichier);
	free(H);
	free(HB);
	free(spinmin);
	fclose(fichier);
}

void print_Emax( spinners_t* spin,  char* add, int Niters)
{
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	FILE* fichier = openfile_out(add);
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);
	int N = spin->nx * spin->ny;
	double Emax = -1000000. * N;
	int* spinmin = (int*)malloc(N * sizeof(int));

	srand((unsigned int)time(NULL));

	for (int i = 0; i < Niters; i++) {
		for (int j = 0; j < N; j++) { spin->angles[j] = rand() % 6; }
		recuit(spin, H, HB, 0.0001, 0.00001, 0.5, 5 * N, 0);
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

void print_neighbours_state_rand(spinners_t* spin, char* add, int Niters, int distance)
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
void print_neighbours_state_all_for(spinners_t* spin, int distancemax, int distance, int* index, int* track, FILE* fichier)
{
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

void recuitN(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu){

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
	for(int i = 0; i < spin->Ngrid ; i++){
		recuit(spin, H, HB, T0, TF, lambda, Niter, N * i);
	}
	remove_equale(spin);
}