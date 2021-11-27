#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct linkedlist linkedlist;
typedef linkedlist* link;
int dimension, k;

link read_file(char fileName[]);
void add_to_linked_list(link *head, double vec[]);
link create_linked_list();
int find_dimension(char line[]);
void initialize(link head, double* mu[]);
double compute_distance(double vec1[],double vec2[]);
void reset_clusters(link head, double* mu[], double* new_sum[]);
double calculating_epsilon(double *mu[], double *new_mu[]);
void create_output(double *mu[], char op_filename[]);

/*
 * missions:
 * 1. check if k is valid value
 * 2. check code on Nova
 * 3. change Linked List - add elements at start of the list
 * 4. there is note in the function - add_to_linked_list
 * 5. warning in function - create_linked_list
 * 6. warning in function - read_file - using atof
 * 7. warning in main - using atoi
 * 8. handle situation where we get name of input file
 * 9. die
 */
int main(int argc, char* argv[]) {
    // assert(argc>=4);
    char* input_file = argv[2];
    char* output_file = argv[3];
    int max_iter = 200;
    double *mu[k], eps, *new_mu[k];
    int i,j;
    link vectors_list = read_file("input_1.txt");
    k = atoi(argv[1]);
    if (argc == 5){
        max_iter = atoi(argv[4]);
    }
    initialize(vectors_list, mu);
    for (i = 0; i < max_iter; ++i) {
        double *new_mu[k];
        reset_clusters(vectors_list, mu, new_mu);
        eps = calculating_epsilon(mu, new_mu);
        for (j = 0; j < k; ++j) {
            mu[j] = new_mu[j];
        }
        if (eps < 0.000001){
            break;
        }
    }
    create_output(mu,output_file);
}

// Defining the data structure Linked List
struct linkedlist{
    double *data;
    struct linkedlist *next;
};

// Creating new linked list and returning the pointer to the head of the list
link create_linked_list(){
    link  head = NULL;
    return head;
}

// Add given vector to a given linked list
void add_to_linked_list(link* head, double vec[]){
    link node = (link) malloc(sizeof (linkedlist));
    link p;
    // if(head->data == NULL){
    //     head->data = vec;
    // }
    node->next = NULL;
    node->data = vec;
    if (*head == NULL){
        *head = node;
    }
    else {
        p = *head;
        while (p->next != NULL)
        {
            p = p->next;
        }
        p->next = node;
    }


}

// Returns a linked list with all the vectors from a given file name
link read_file(char fileName[]){
    FILE *file = fopen(fileName,"r");
    char buff[255], copy_buff[255], *ptr;
    link list_of_vectors;
    int ch;
    double* vector, *place;
    if (file) {
        list_of_vectors = create_linked_list();
        ch = fscanf(file, "%s", buff);
        strcpy(copy_buff, buff);
        dimension = find_dimension(copy_buff);
        while ((ch != '\n') && (ch != EOF)) {
            vector = (double *)(malloc(dimension * sizeof (double)));
            place = vector;
            ptr = strtok(buff, ",");
            while (ptr != NULL) {
                *(vector) = atof(ptr);
                ptr = strtok(NULL, ",");
                vector++;
            }
            add_to_linked_list(&list_of_vectors, place);
            ch = fscanf(file, "%s", buff);
        }
        fclose(file);
        return list_of_vectors;
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
void initialize(link head, double* mu[]){
    int i;
    double* vec;
    for (i = 0; i < k; ++i) {
        mu[i] = (double *)calloc(dimension, sizeof (double));
        vec = head->data;
        mu[i] = vec;
        head = head->next;
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
void reset_clusters(link head, double* mu[], double* new_sum[]) {
    int count[k];
    double *vec;
    int i, j, r, s, q, min_mu;

    for (i = 0; i < k; ++i) {
        new_sum[i] = (double *) calloc(dimension, sizeof(double));
        count[i] = 0;
    }
    while (head != NULL) {
        min_mu = calc_argmin(mu, head->data);
        count[min_mu]++;
        vec = head->data;
        for (j = 0; j < dimension; ++j) {
            new_sum[min_mu][j] += *vec;
            vec++;
        }
        head = head->next;
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

