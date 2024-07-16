from generate_dataset import generate_reference_models
import numpy as np
import random
import plotly.express as px

np.random.seed(0)
random.seed(0)


def get_s_set(precision_mat1,precision_mat2):
    delta_star = precision_mat1 - precision_mat2
    s_set_edges = np.argwhere(delta_star != 0)
    s_set_edges = [[i0, i1] for i0, i1 in s_set_edges.tolist() if i0 != i1]
    s_set_nodes = [x for xs in s_set_edges for x in xs]
    s_set_nodes = list(np.unique(s_set_nodes))
    return s_set_nodes, s_set_edges


def get_s_complement_set(precision_mat1,precision_mat2):
    delta_star = precision_mat1 - precision_mat2
    np.fill_diagonal(delta_star, 1)
    s_complement_edges = np.argwhere(delta_star == 0)
    all_nodes = list(range(precision_mat1.shape[0]))
    s_set_nodes, _ = get_s_set(precision_mat1, precision_mat2)
    s_complement_nodes = [i for i in all_nodes if i not in s_set_nodes]
    return s_complement_nodes, s_complement_edges.tolist()


def get_s_complement_samples(s_complement_set, sample_size):
    max_number = min(len(s_complement_set), sample_size)
    random_indexes = random.sample(range(len(s_complement_set)), max_number)
    return [s_complement_set[i] for i in random_indexes]
    # seen_indexes = []
    # for i in s_complement_set:
    #     if list(i[::-1]) in seen_indexes or i[0] == i[1]:
    #         continue
    #     seen_indexes.append(list(i))
    # max_number = min(len(seen_indexes), sample_size)
    # random_lines = random.sample(range(len(seen_indexes)), max_number)
    # return [seen_indexes[i] for i in random_lines]


def get_s_prime_set(s_set, s_complement_set, sample_size):
    s_complement_subset = get_s_complement_samples(s_complement_set, sample_size)
    s_prime_set = s_complement_subset + s_set
    # s_set = [list(i) for i in s_set]
    # s_prime_set = []
    # for i in s_set + s_complement_subset:
    #     if i[::-1] in s_prime_set:
    #         continue
    #     s_prime_set.append(i)
    return s_prime_set


def split_dataset(samples_a, samples_b, split_ratio=.5):
    ols_samples = random.sample(range(samples_a.shape[0]), round(split_ratio*samples_a.shape[0]))
    init_lasso_samples = np.array(range(samples_a.shape[0]))
    lasso_samples = np.delete(init_lasso_samples, ols_samples, axis=0)
    # lasso_samples = [i for i in range(samples_a.shape[0]) if i not in ols_samples]
    ols_dataset_a, ols_dataset_b = samples_a[ols_samples, :], samples_b[ols_samples, :]
    lasso_dataset_a, lasso_dataset_b = samples_a[lasso_samples, :], samples_b[lasso_samples, :]
    return ols_dataset_a, ols_dataset_b, lasso_dataset_a, lasso_dataset_b


def split_dataset_with_features(samples_a, samples_b, features, split_ratio=0.5):
    ols_samples = random.sample(range(samples_a.shape[0]), round(split_ratio*samples_a.shape[0]))
    init_lasso_samples = np.array(range(samples_a.shape[0]))
    lasso_samples = np.delete(init_lasso_samples, ols_samples, axis=0)
    # lasso_samples = [i for i in range(samples_a.shape[0]) if i not in ols_samples]
    ols_dataset_a, ols_dataset_b = samples_a[np.ix_(ols_samples, features)], samples_b[np.ix_(ols_samples, features)]
    lasso_dataset_a, lasso_dataset_b = samples_a[np.ix_(lasso_samples, features)], samples_b[np.ix_(lasso_samples, features)]
    return ols_dataset_a, ols_dataset_b, lasso_dataset_a, lasso_dataset_b


def get_delta_hat(dataset_a, dataset_b):
    cov_mat_a = np.cov(dataset_a.T)
    precision_mat_a = np.linalg.inv(cov_mat_a)
    cov_mat_b = np.cov(dataset_b.T)
    precision_mat_b = np.linalg.inv(cov_mat_b)
    return precision_mat_a - precision_mat_b


def generate_mirror_statistics(delta_hat_ols, delta_hat_lasso):
    squarer = lambda t: t ** 2
    exp_f = lambda t: np.e ** t
    vfunc = np.vectorize(squarer)
    ols_delta_hat_squere = vfunc(delta_hat_ols)
    lasso_delta_hat_squere = vfunc(delta_hat_lasso)
    sum_mat = ols_delta_hat_squere + lasso_delta_hat_squere
    # exp_func = np.vectorize(exp_f)
    # exp_mat = exp_func(sum_mat)

    # get sign
    ols_diag = np.diag(delta_hat_ols)
    lasso_diag = np.diag(delta_hat_lasso)
    empty_mat = np.zeros(delta_hat_lasso.shape)
    for i in range(delta_hat_ols.shape[0]):
        for i1 in range(delta_hat_ols.shape[1]):
            empty_mat[i, i1] = 2*delta_hat_ols[i, i1] + ols_diag[i] + lasso_diag[i1]
    sign_mat = np.sign(delta_hat_lasso * empty_mat)
    #calculate mirror mat
    mirror_mat = sign_mat * sum_mat
    return mirror_mat


def get_cutoff_value(mirror_statistics, fdr_level):
    mirror_mat = mirror_statistics.copy()
    # for i in range(mirror_mat.shape[0]):
    #     mirror_mat[i, i] = 0
    np.fill_diagonal(mirror_mat, 0)
    # check_points1 = [-(i - 1e-6) for i in list(np.unique(mirror_mat)) if i < 0]
    check_points = [i - 1e-6 for i in list(np.unique(mirror_mat)) if i > 0]
    # check_points = sorted(check_points1 + check_points2)
    check_points = sorted(check_points, reverse=True)
    count = 0
    # print('++++++++++++++++', len(check_points))
    for point in check_points:
        numerator = (mirror_mat < -point).sum()
        denumerator = max((mirror_mat > point).sum(), 1)
        fdr = numerator/denumerator
        # print(point, numerator, denumerator, fdr)
        if fdr <= fdr_level:
            # print(point, numerator, denumerator, fdr)
            return point
    return None


def get_features_mapping(features):
    mapping_dict = {}
    for key in range(len(features)):
        mapping_dict[key] = features[key]
    return mapping_dict


def intersect_lists_of_lists(list1, list2):
    set1 = set(map(tuple, list1))
    set2 = set(map(tuple, list2))
    intersected_set = set1.intersection(set2)
    return [list(item) for item in intersected_set]


def get_fpr_tpr(H0_edges, H1_edges, s_prime_set, mirror_statistics, cutoff_value):
    mirror_mat = mirror_statistics.copy()
    np.fill_diagonal(mirror_mat, 0)
    feature_map = get_features_mapping(s_prime_set)
    # differential_nodes_indexes = [x for xs in np.argwhere(mirror_statistics > cutoff_value).tolist() for x in xs]
    differential_elements = np.argwhere(mirror_statistics > cutoff_value)
    s_hat_cutoff = [[feature_map[i0], feature_map[i1]] for i0, i1 in differential_elements.tolist()]

    # differential_nodes_indexes = np.unique(differential_nodes_indexes).tolist()
    # differential_nodes = [feature_map[i] for i in differential_nodes_indexes]
    # tp = list(set(differential_nodes) & set(s_set))
    # fp = list(set(differential_nodes) & set(s_complement_set))
    # print('----------')
    # print(s_hat_cutoff)
    # print(H0_edges)
    # print('---------')
    fp = intersect_lists_of_lists(s_hat_cutoff, H0_edges)
    tp = intersect_lists_of_lists(s_hat_cutoff, H1_edges)
    denum = max(1, len(s_hat_cutoff))
    print('tp len', len(tp))
    tpr, fpr = len(tp)/len(H1_edges), len(fp)/len(H0_edges)
    return fpr, tpr


def plot_fpr_tpr_fit(H0_edges, H1_edges, s_prime_set, mirror_statistics, sample_number=100):
    q_list = [i / sample_number for i in range(1, sample_number+1)]
    fpr_list, tpr_list = [], []
    for i in range(len(q_list)):
        q = q_list[i]
        if i % 10 == 0:
            print(i, q)
        # print('cutoff value')
        cutoff_value = get_cutoff_value(mirror_statistics, q)

        # print('---------')
        fpr, tpr = get_fpr_tpr(H0_edges, H1_edges, s_prime_set, mirror_statistics, cutoff_value)
        print('fpr', fpr, 'tpr', tpr)
        fpr_list.append(fpr)
        tpr_list.append(tpr)

    print(len(fpr_list))
    fig = px.line(x=fpr_list, y=tpr_list, text=q_list)
    fig.update_layout(
        title='FPR vs TPR with different q',
        xaxis_title='FPR',
        yaxis_title='TPR')
    fig.show()


def get_e_values(features_size, mirror_statistics, cutoff_value):
    mirror_mat = mirror_statistics.copy()
    np.fill_diagonal(mirror_mat, 0)
    numerator = (mirror_mat >= cutoff_value).astype(int)
    denumerator = (mirror_mat <= cutoff_value).sum() + 1
    e_vals = features_size*numerator/denumerator
    return e_vals


def get_k_hat(alpha, e_vals):
    vec_evals = e_vals.ravel()
    vec_evals[::-1].sort()
    for i in range(len(vec_evals)):
        if vec_evals[i] < len(vec_evals)/(alpha*(i+1)):
            return i
    return i + 1


def get_e_values_avg(e_values_list):
    return np.mean(e_values_list, axis=0)


def get_s_hat_mds(e_values_avg, k_hat, q_mds, features_size):
    avg_mat = e_values_avg.copy()
    np.fill_diagonal(avg_mat, 0)
    thr = features_size/(k_hat * q_mds)
    s_hat_mds = np.argwhere(avg_mat >= thr)
    return s_hat_mds.tolist()


def run_different_number_of_changes():
    change_numbers = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for number in change_numbers:
        print(number, '---------------------------------------')
        datasets = generate_reference_models(500, 500, number)
        s_set, H1 = get_s_set(datasets[1], datasets[4])
        s_complement, H0 = get_s_complement_set(datasets[1], datasets[4])
        s_prime_set = get_s_prime_set(s_set, s_complement, int(number/5))
        ols_a, ols_b, lasso_a, lasso_b = split_dataset_with_features(datasets[0], datasets[3], s_prime_set, .5)
        ols_delta_hat = get_delta_hat(ols_a, ols_b)
        lasso_delta_hat = get_delta_hat(lasso_a, lasso_b)
        mirror_stats = generate_mirror_statistics(ols_delta_hat, lasso_delta_hat)
        cutoff = get_cutoff_value(mirror_stats, .1)
        t1, t2 = get_fpr_tpr(H0, H1, s_prime_set, mirror_stats, cutoff)
        print('fpr: ', t1, 'tpr: ', t2)


# run_different_number_of_changes()
# exit()

test = generate_reference_models(500, 500, 50, mult=1)
print('data is generated')
s_set, H1 = get_s_set(test[1], test[4])
s_complement, H0 = get_s_complement_set(test[1], test[4])
print('s_set')
print(len(s_set), len(H1))
print('s_complement_set')
print(len(s_complement), len(H0))
# print('s_complement_samples')
# s_comp_subset = get_s_complement_samples(s_complement, 3)
# print(s_comp_subset)
print('s_prime_set')
s_prime_set = get_s_prime_set(s_set, s_complement, 20)
# s_prime_set = get_s_prime_set(s_set, s_complement, 1)
print(len(s_prime_set))
print('dataset spilit features')
ols_a, ols_b, lasso_a, lasso_b = split_dataset_with_features(test[0], test[3], s_prime_set, .5)
print('delta hat')
ols_delta_hat = get_delta_hat(ols_a, ols_b)
lasso_delta_hat = get_delta_hat(lasso_a, lasso_b)
print(ols_delta_hat.max(), ols_delta_hat.min())
print(lasso_delta_hat.max(), lasso_delta_hat.min())
print('mirror stats')
mirror_stats = generate_mirror_statistics(ols_delta_hat, lasso_delta_hat)
# print(mirror_stats)

# print('cutoff value')
cutoff = get_cutoff_value(mirror_stats, .1)
# print(cutoff)

print('---------')
# t1, t2 = get_fpr_tpr(H0, H1, s_prime_set, mirror_stats, cutoff)
# print('fpr: ', t1, 'tpr: ', t2)
# plot_fpr_tpr_fit(H0, H1, s_prime_set, mirror_stats, 100)
# test = get_e_values(test[1].shape[1], mirror_stats, cutoff)
# print(test.shape)