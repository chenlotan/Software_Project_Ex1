#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int dimension, k, N;


double **read_file(char fileName[]);

void number_of_vectors(char fileName[]);

void find_dimension(char line[]);

void initialize(double **vectors_list, double *mu[], double *new_mu[]);

double compute_distance(double vec1[], double vec2[]);

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[]);

double calculating_epsilon(double *mu[], double *new_mu[]);

void create_output(double *mu[], char op_filename[]);

int check_allocation(const double *p);

void free_memory(double **array, int len);

int main(int argc, char *argv[]) {
    char *input_file = argv[argc - 2];
    char *output_file = argv[argc - 1];
    int max_iter = 200;
    int i, j, q;
    double eps;
    double **new_mu, **mu;
    double **vectors_list;
    number_of_vectors(input_file);
    vectors_list = read_file(input_file);
    k = (int) strtol(argv[1], NULL, 10);
    if (argc == 5) {
        max_iter = (int) strtol(argv[2], NULL, 10);
    }
    if (k < 1 || k > N) {
        printf("Invalid Input!");
        exit(1);
    }
    mu = (double **) malloc(k * sizeof(double *));
    new_mu = (double **) malloc(k * sizeof(double *));
    initialize(vectors_list, mu, new_mu);
    for (i = 0; i < max_iter; ++i) {
        reset_clusters(vectors_list, mu, new_mu);
        eps = calculating_epsilon(mu, new_mu);
        for (j = 0; j < k; ++j) {
            for (q=0; q<dimension; ++q){
                mu[j][q] = new_mu[j][q];
                new_mu[j][q] = 0;
            }
        }
        if (eps < 0.000001) {
            break;
        }
    }
    create_output(mu, output_file);
    free_memory(vectors_list, N);
    free_memory(mu, k);
    free_memory(new_mu, k);
    return 0;
}


double **read_file(char fileName[]) {
    FILE *file = fopen(fileName, "r");
    char buff[1024], *ptr;
    int ch, i = 0;
    double *vector, *place;
    double **all_vectors;
    if (file) {
        all_vectors = (double **) malloc(N * sizeof(double *));
        ch = fscanf(file, "%s", buff);
        while ((ch != '\n') && (ch != EOF)) {
            vector = (double *) (malloc(dimension * sizeof(double)));
            check_allocation(vector);
            place = vector;
            ptr = strtok(buff, ",");
            while (ptr != NULL) {
                *(vector) = strtod(ptr, NULL);
                ptr = strtok(NULL, ",");
                vector++;
            }
            all_vectors[i] = place;
            i++;
            ch = fscanf(file, "%s", buff);
        }
        fclose(file);
        return all_vectors;
    } else {
        printf("Invalid Input!");
        exit(1);
    }

}

double compute_distance(double vec1[], double vec2[]) {
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        sum += (double) pow(vec1[i] - vec2[i], 2);
    }
    return sum;
}

void find_dimension(char line[]) {
    int i = 0;
    char *ptr = strtok(line, ",");
    while (ptr != NULL) {
        ptr = strtok(NULL, ",");
        i++;
    }
    dimension = i;
}


void number_of_vectors(char fileName[]){
    FILE *file = fopen(fileName, "r");
    char buff[1024];
    int ch, n = 0;
    ch = fscanf(file, "%s", buff);
    while((ch != 'n') && (ch != EOF)){
        if (n == 0){
            find_dimension(buff);
        }
        n++;
        ch = fscanf(file, "%s", buff);
    }
    N = n;
    fclose(file);
}

void initialize(double **vectors_list, double *mu[], double *new_mu[]) {
    int i, j;
    double *vec;
    for (i = 0; i < k; ++i) {
        mu[i] = (double *) calloc(dimension, sizeof(double));
        new_mu[i] = (double *) calloc(dimension, sizeof(double));
        check_allocation(mu[i]);
        check_allocation(new_mu[i]);
        vec = vectors_list[i];
        for (j = 0; j < dimension; j++)
        {
            mu[i][j] = vec[j];
        }
    }
}

int calc_argmin(double *mu[], double *vector) {
    double min_val = HUGE_VAL, sum_p;
    int min_mu = 0, i;
    for (i = 0; i < k; ++i) {
        sum_p = compute_distance(mu[i], vector);
        if (sum_p < min_val) {
            min_val = sum_p;
            min_mu = i;
        }
    }
    return min_mu;
}

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[]) {
    int *count;
    double *vec;
    int i, t, j, r, s, q, z, min_mu;
    count = calloc(k, sizeof(int));

    for (t = 0; t < N; t++) {
        min_mu = calc_argmin(mu, vectors_list[t]);
        count[min_mu]++;
        vec = vectors_list[t];
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

double calculating_epsilon(double *mu[], double *new_mu[]) {
    double eps = 0, dist;
    int i;
    for (i = 0; i < k; ++i) {
        dist = compute_distance(mu[i], new_mu[i]);
        if (eps < dist) {
            eps = dist;
        }
    }
    return eps;
}

void create_output(double *mu[], char op_filename[]) {
    FILE *f;
    int i, j;
    f = fopen(op_filename, "w");
    for (i = 0; i < k; ++i) {
        for (j = 0; j < dimension; ++j) {
            fprintf(f, "%0.4f", mu[i][j]);
            if (j != dimension - 1) {
                fprintf(f, ",");
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
}


int check_allocation(const double *p) {
    if (p == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    return 0;
}


void free_memory(double **array, int len) {
    int q;
    for (q = 0; q < len; q++) {
        free(array[q]);
    }
    free(array);
}