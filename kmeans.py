import math
import sys


def validate(condition):
    if not condition:
        print('Invalid Input!')
        exit(1)


# Open the file and create array of vectors
def read_file(file_name):
    f = open(file_name, "r")
    vectors = []
    for line in f:
        if line != "\n":
            numbers = line.split(',')
            for i in range(len(numbers)):
                numbers[i] = float("%0.4f" % float(numbers[i]))
            vectors.append(numbers)
    return vectors


# Creating K clusters
def initialize(vectors, k):
    return vectors[0:k]


# Compute the best centroid for certain point
def calc_argmin(mu, point):
    min_val = math.inf
    min_mu = 0
    for i in range(len(mu)):
        sum_p = compute_distance(point, mu[i])
        if sum_p < min_val:
            min_val = sum_p
            min_mu = i
    return min_mu


# Compute new centroids
def calc_new_centroids(vectors, mu, k, d):
    sum_by_mu = [[0 for j in range(d)] for i in range(k)]
    count_by_mu = [0 for i in range(k)]
    for i in range(0, len(vectors)):
        min_mu = calc_argmin(mu, vectors[i])
        for j in range(d):
            sum_by_mu[min_mu][j] += vectors[i][j]
        count_by_mu[min_mu] += 1
    for i in range(k):
        for j in range(d):
            if count_by_mu[i] != 0:
                sum_by_mu[i][j] = sum_by_mu[i][j] / count_by_mu[i]
            else:
                sum_by_mu[i][j] = mu[i][j]
    return sum_by_mu


# Compute epsilon
def calc_eps(old_mu, new_mu):
    epsilon = 0
    for i in range(len(old_mu)):
        dist = compute_distance(old_mu[i], new_mu[i])
        if epsilon < dist:
            epsilon = dist
    return epsilon


# Compute distance between two given vectors
def compute_distance(vec1, vec2):
    assert len(vec1) == len(vec2)
    dist = "%.4f" % ((sum([(vec1[i] - vec2[i]) ** 2 for i in range(len(vec1))])) ** 0.5)
    return float(dist)


# Writing to the given file the final centroids
def create_output(mu, op_file_name):
    output = open(op_file_name, "w")
    for mu_i in mu:
        line = ""
        for i in range(len(mu_i) - 1):
            line += "%0.4f" % mu_i[i] + ","
        line += "%0.4f" % mu_i[len(mu_i) - 1]
        output.write(line)
        output.write("\n")


args = sys.argv[1:]
validate(3 <= len(args) <= 4)
k, max_iter, input_filename, output_filename = args[0], args[1] if len(args) == 4 else 200, args[-2], args[-1]
validate(args[0].isdigit())
if len(args) == 4:
    validate(max_iter.isdigit())
    max_iter = int(max_iter)
k = int(k)
vectors = read_file(input_filename)
d = len(vectors[0])
mu = initialize(vectors, k)
validate(1 < k < len(vectors))
for i in range(max_iter):
    new_mu = calc_new_centroids(vectors, mu, k, d)
    epsilon = calc_eps(mu, new_mu)
    if epsilon < 0.001:
        break
    mu = new_mu.copy()
create_output(mu, output_filename)
