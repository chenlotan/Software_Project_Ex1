#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct linkedlist linkedlist;
typedef linkedlist* link;
int dimension, k, N;
double** read_file(char fileName[]);
int find_dimension(char line[]);
void initialize(double** vectors_list, double* mu[]);
double compute_distance(double vec1[],double vec2[]);
void reset_clusters(double** vectors_list, double* mu[], double* new_sum[]);
double calculating_epsilon(double *mu[], double *new_mu[]);
void create_output(double *mu[], char op_filename[]);
int get_len_linked_list(link head);
int check_allocation(double* p);
void free_memory(double** array, int len);

/*
 * missions:
 * 2. check code on Nova
 * 6. warning in function - read_file - using atof
 * 7. warning in main - using atoi
 * 9. check if memory got released 
 * die
 */

int main(int argc, char* argv[]) {
    char* input_file = argv[2];
    char* output_file = argv[3];
    int max_iter = 200;
    int i,j,q;
    double** vectors_list = read_file(input_file);
    k = (int )strtol(argv[1], '\0', 10);
    double *mu[k], eps, *new_mu[k];
    if (argc == 5){
        max_iter = (int )strtol(argv[4], '\0', 10);
    }
    if (k<1 || k>N){
        printf("Invalid Input!");
        exit(1);
    }
    initialize(vectors_list, mu);
    for (i = 0; i < max_iter; ++i) {
        double *new_mu[k];
        reset_clusters(vectors_list, mu, new_mu);
        eps = calculating_epsilon(mu, new_mu);
        for (j = 0; j < k; ++j) {
            mu[j] = new_mu[j];
        }
        free_memory(new_mu, k);
        if (eps < 0.000001){
            break;
        }
    }
    create_output(mu,output_file);
    free_memory(vectors_list, N);
    free_memory(mu, k);
    return 0;
}




// Returns a linked list with all the vectors from a given file name
double** read_file(char fileName[]){
    FILE *file = fopen(fileName,"r");
    char buff[255], copy_buff[255], *ptr;
    int ch, max_size = 100, n = 0; //--------max_size is the current size of all_vectors, n is the real amount of vectors we entered---------
    double* vector, *place;
    double** all_vactors;
    if (file) {
        all_vactors = (double **)malloc(max_size * sizeof(double*));
        ch = fscanf(file, "%s", buff);
        strcpy(copy_buff, buff);
        dimension = find_dimension(copy_buff);
        while ((ch != '\n') && (ch != EOF)) {
            vector = (double *)(malloc(dimension * sizeof (double)));
            check_allocation(vector);
            place = vector;
            ptr = strtok(buff, ",");
            while (ptr != NULL) {
                *(vector) = atof(ptr);
                ptr = strtok(NULL, ",");
                vector++;
            }
            if (n >= max_size){
                max_size *= 1.5;
                all_vactors = (double **) realloc(all_vactors, max_size * sizeof(double*));
            }
            all_vactors[n] = place;
            n++;
            ch = fscanf(file, "%s", buff);
        }
        fclose(file);
        if (n<max_size){
            all_vactors = (double **)realloc(all_vactors, n*sizeof(double*)); // ---------if we entered less veectors than cuurent size of all_vectors - realloc to array in size n---------------
        }
        N = n;
        return all_vactors;
    }
    else{
        printf("Invalid Input!");
        exit(1);
    }

}

// Returns the distance between two given vectors
double compute_distance(double vec1[], double vec2[]){
    double sum = 0, dist;
    int i;
    for (i = 0; i < dimension; i++){
        sum += (double)pow(vec1[i] - vec2[i], 2);
    }
    dist = (double)pow(sum, 0.5);
    return dist;
}

// Compute the dimension of the vectors in the file
int find_dimension(char line[]){
    int i = 0;
    char *ptr = strtok(line, ",");
    while (ptr != NULL) {
        ptr = strtok(NULL, ",");
        i++;
    }
    return i;
}

// Initializing the first k centroids from the list of all the vectors
void initialize(double** vectors_list, double* mu[]){
    int i;
    double* vec;
    for (i = 0; i < k; ++i) {
        mu[i] = (double *)calloc(dimension, sizeof (double));
        check_allocation(mu[i]);
        vec = vectors_list[i];
        mu[i] = vec;
    }
}

// Returns the best centroid's index for a given vector
int calc_argmin(double *mu[], double *vector){
    double min_val = INFINITY, sum_p;
    int min_mu = 0, i;
    for (i = 0; i < k; ++i) {
        sum_p = compute_distance(mu[i], vector);
        if (sum_p < min_val){
            min_val = sum_p;
            min_mu = i;
        }
    }
    return min_mu;
}

// Compute and return the new centroids
void reset_clusters(double** verctors_list, double* mu[], double* new_sum[]) {
    int count[k];
    double *vec;
    int i, t, j, r, s, q, min_mu;

    for (i = 0; i < k; ++i) {
        new_sum[i] = (double *) calloc(dimension, sizeof(double));
        check_allocation(new_sum[i]);
        count[i] = 0;
    }
    for (t = 0; t < N; t++) {
        min_mu = calc_argmin(mu, verctors_list[t]);
        count[min_mu]++;
        vec = verctors_list[t];
        for (j = 0; j < dimension; ++j) {
            new_sum[min_mu][j] += *vec;
            vec++;
        }
    }
    for (r = 0; r < k; ++r) {
        if (count[r] == 0) {
            for (s = 0; s < dimension; ++s) {
                new_sum[r][s] = mu[r][s];
            }
        } else {
            for (q = 0; q < dimension; ++q) {
                new_sum[r][q] = new_sum[r][q] / (double) count[r];
            }
        }
    }
    free(count);
}

// Calculating the new epsilon between the new centroids list to the old ones
double calculating_epsilon(double *mu[], double *new_mu[]){
    double eps = 0, dist;
    int i;
    for (i = 0; i < k; ++i) {
        dist = compute_distance(mu[i], new_mu[i]);
        if (eps < dist){
            eps = dist;
        }
    }
    return eps;
}

// Writing to a given file name the centroids
void create_output(double *mu[], char op_filename[]){
    FILE *f;
    int i,j;
    f = fopen(op_filename, "w");
    for (i = 0; i < k; ++i) {
        for (j = 0; j < dimension; ++j) {
            fprintf(f, "%0.4f", mu[i][j]);
            if (j != dimension - 1){
                fprintf(f, ",");
            }
        }
        fprintf(f,"\n");
    }
    fclose(f);
}


// Check if allocation of memory worked for double pointer
int check_allocation(double* p){
    if (p == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    return 0;
}


void free_memory(double** array, int len){
    int q;
    for (q=0; q<len; q++){      
        free(array[q]); 
    }
    free(array);
}