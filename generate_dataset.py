import numpy as np
from rpy2 import robjects
from rpy2.robjects import pandas2ri


# Output sequence is samples_A, precision_matrix_A, covariance_matrix_A,
# samples_B,  precision_matrix_B, covariance_matrix_B
def generate_reference_models(number_of_nodes, number_of_samples, number_of_changes, type="ScaleFree",
                              density_of_graph=0.2, power=1, mult=1):
    pandas2ri.activate()
    robjects.r.source('R_codes/Libraray.R')
    r_func = robjects.globalenv['generate_reference_models']
    output = r_func(number_of_nodes, number_of_samples, type, density_of_graph, power, number_of_changes, mult)
    samples_A, precision_mat_A, cov_mat_A = np.array(output[0]), np.array(output[1]), np.array(output[2])
    samples_B, precision_mat_B, cov_mat_B = np.array(output[3]), np.array(output[4]), np.array(output[5])
    return samples_A, precision_mat_A, cov_mat_A, samples_B, precision_mat_B, cov_mat_B


def base_differential_network(number_of_nodes, number_of_samples=1, decay_ratio = 0.3, alternative_value = 0.3,
                              numberOfChanges = 2):
    pandas2ri.activate()
    robjects.r.source('R_codes/Libraray.R')
    r_func = robjects.globalenv['base_differential_network']
    output = r_func(number_of_nodes, number_of_samples, decay_ratio, alternative_value, alternate_diag_index)
    samples_A, precision_mat_A, cov_mat_A = np.array(output[0]), np.array(output[1]), np.array(output[2])
    samples_B, precision_mat_B, cov_mat_B = np.array(output[3]), np.array(output[4]), np.array(output[5])
    return samples_A, precision_mat_A, cov_mat_A, samples_B, precision_mat_B, cov_mat_B
