from generate_dataset import generate_reference_models
import numpy as np
import random

np.random.seed(17)
random.seed(17)


def get_kronecker_value(mat1, mat2, i, j):
    d = mat1.shape[0]
    return mat1[i//d, j//d] * mat2[i%d, j%d]


def get_kronecker_col(mat1, mat2, col):
    d = mat1.shape[0]
    return [get_kronecker_value(mat1, mat2, i, col) for i in range(d**2)]


def get_kronecker_multi_col(mat1, mat2, indices):
    col_list = []
    for idx in indices:
        col_list.append(get_kronecker_col(mat1, mat2, idx))
    col_sub_kron = np.array(col_list).T
    return col_sub_kron


def get_kronecker_sub_mat(mat1, mat2, indices):
    all_rows_kronecker = get_kronecker_multi_col(mat1, mat2, indices)
    return all_rows_kronecker[np.ix_(indices)]


def format_float(x, d):
    return float(f"{x:.{d}f}")


def round_matrix(mat, d):
    vectorized_format_float = np.vectorize(format_float)
    formatted_matrix = vectorized_format_float(mat, d)
    return formatted_matrix


# def get_va_sa(v,  indices):
#     va = v[np.ix_(indices)]
#     sa = np.sign(va)
#     return va, sa


# def calc_hit_value(idx, cols_kron, val1, val2, v):
    # print(idx, cols_kron.shape, val1.shape, val2.shape)
    # numerator = cols_kron[idx] @ val1 - v[idx]
    # return max(numerator/(1 - cols_kron[idx] @ val2), numerator/(-1 - cols_kron[idx] @ val2))

# def get_lambda_hit(cols_kron, val1, val2, v, active_complement):
#     print('aaaa', cols_kron.shape, val1.shape, val2.shape)
#     active_complement = np.array(active_complement)
#     vectorized_function = np.vectorize(calc_hit_value)
#     results = np.array([calc_hit_value(i, cols_kron, val1, val2, v) for i in active_complement])
#     return np.max(results)


def get_lambda_hit(indices, cols_kron, val1, val2, v, lambda_t):
    final_result = np.zeros((cols_kron.shape[0], 2))
    cols_kron_idx = cols_kron[indices]  # Select rows corresponding to indices
    numerator = np.dot(cols_kron_idx, val1) - v[indices]

    denom1 = (1 - np.dot(cols_kron_idx, val2)).reshape(-1, 1)
    denom2 = (-1 - np.dot(cols_kron_idx, val2)).reshape(-1, 1)
    result1 = numerator / denom1
    result2 = numerator / denom2
    cat_results = np.concatenate((result2, result1), axis=1)
    final_result[indices] = cat_results
    cat_results = round_matrix(final_result, 5)
    cat_results[cat_results >= lambda_t] = 0

    flatten_cat = cat_results.flatten()
    # print('jjkkkkk', sorted(flatten_cat, reverse=True)[:5])
    # print(cat_results[np.ix_([13, 20])])
    # print(result1.shape, result2.shape, cat_results.shape)
    max_value = np.max(cat_results)
    max_indices = np.where(cat_results == max_value)
    # max_indices = list(zip(max_indices[0], max_indices[1]))
    indices_sign = np.sign(max_indices[1] - .5)
    max_indices = list(zip(max_indices[0], indices_sign))
    sorted_indices = sorted(max_indices, key=lambda x: x[0])
    # print(sorted_indices)
    return max_value, sorted_indices


def get_lambda_cross(val1, val2, lambda_t):
    cross_vals = - val1/val2
    cross_vals = round_matrix(cross_vals, 5)
    cross_vals[cross_vals >= lambda_t] = 0
    max_value = np.max(cross_vals)
    max_indices = np.where(cross_vals == max_value)[0]
    return max_value, list(max_indices)


def convert_line_index_to_matrix(line_index, d):
    matrix_index = [(i%d, i//d) for i in line_index]
    return matrix_index


def alg_1(sigma, sigma_hat, active_set_pack, lambda_t):
    active_set = [i[0] for i in active_set_pack]
    sa = [i[1] for i in active_set_pack]
    d = sigma.shape[0]
    active_complement = [i for i in range(d**2) if i not in active_set]
    v = (sigma - sigma_hat).flatten('F').reshape(-1, 1)
    if len(active_set) == 0:
        v_val_cat = np.concatenate((v, -v))
        v_val_cat = round_matrix(v_val_cat, 5)
        next_lambda = np.max(v_val_cat)
        max_indices = np.where(v_val_cat == next_lambda)
        indices_sign = np.sign(max_indices[1] - .5)
        max_indices = list(zip(max_indices[0], indices_sign))
        sorted_indices = sorted(max_indices, key=lambda x: x[0])
        # new_active_set = [i[0] for i in sorted_indices]
        # new_sa = [i[1] for i in sorted_indices]
        return next_lambda, sorted_indices
    else:
        va = v[np.ix_(active_set)]
        cols_kron_1 = get_kronecker_multi_col(sigma, sigma_hat, active_set)
        cols_kron_2 = get_kronecker_multi_col(sigma_hat, sigma, active_set)


        sub_kron_1 = get_kronecker_sub_mat(sigma, sigma_hat, active_set)
        sub_kron_2 = get_kronecker_sub_mat(sigma_hat, sigma, active_set)

        cols_kron = (cols_kron_1+cols_kron_2)/2
        sub_kron = (sub_kron_1 + sub_kron_2)/2

        sub_kron_inverse = np.linalg.inv(sub_kron)
        val1 = sub_kron_inverse @ va
        val2 = (sub_kron_inverse @ sa).reshape(-1, 1)
        # if np.any(val2 == 0):
        #     val2 += 1e-8
        lambda_hit_val, lambda_hit_idx = get_lambda_hit(active_complement, cols_kron, val1, val2, v, lambda_t)
        lambda_cross_val, lambda_cross_idx = get_lambda_cross(val1, val2, lambda_t)
        lambda_cross_idx = [active_set[i] for i in lambda_cross_idx]

        # lambda_hit_val = min(lambda_hit_val, lambda_t)
        # lambda_cross_val = min(lambda_cross_val, lambda_t)
        next_lambda = max(lambda_hit_val, lambda_cross_val)
        print('Calc Values: ', lambda_t, next_lambda, lambda_hit_val, lambda_cross_val)
    if next_lambda <= 0:
        return lambda_t, active_set_pack
    # print(lambda_cross_val, lambda_cross_idx)
    if lambda_hit_val > lambda_cross_val:
        # print('hhiiiitttt', active_set_pack + lambda_hit_idx)
        new_active_pack = sorted(active_set_pack + lambda_hit_idx, key=lambda x: x[0])
        print(f'hittt, lambda hit: {lambda_hit_val}, lambda cross: {lambda_cross_val}, lambda t: {lambda_t}')
    else:
        new_active_set = [i for i in active_set if i not in lambda_cross_idx]
        new_active_pack = [i for i in active_set_pack if i[0] in new_active_set]
        print(f'cross, lambda hit: {lambda_hit_val}, lambda cross: {lambda_cross_val}, lambda t: {lambda_t}')

    return next_lambda, new_active_pack


dataset_1, precision_1, cov_1, dataset_2, precision_2, cov_2 = generate_reference_models(5, 5, 2, mult=1)
delta_star = precision_2 - precision_1
delta_star_min = np.abs(delta_star).min()


real_cov1 = np.linalg.inv(precision_1)
real_cov2 = np.linalg.inv(precision_2)

print((precision_1 - precision_2))
print((precision_1 - precision_2).flatten('F').reshape(-1, 1))

a_set_sa, lam, stop = [], 1000000, False
for i in range(20):
    new_lam, new_a_set_sa = alg_1(cov_1, cov_2, a_set_sa, lam)
    print('new_lam', new_lam)
    if new_lam <= 0:
        break
    a_set_sa, lam = new_a_set_sa, new_lam
    print(lam, a_set_sa)
    # print(convert_line_index_to_matrix([a[0] for a in a_set_sa], 5))
