from __future__ import division
import sys, copy, bz2
import numpy as np


depth_filename = sys.argv[1]
snp_filename = sys.argv[2]
indel_filename = sys.argv[3]
likelihood_filename = sys.argv[4]

# where to read the input junctions from
depth_file = bz2.open(depth_filename, "rt")
snp_file = bz2.open(snp_filename, "rt")
indel_file = bz2.open(indel_filename, "rt")
likelihood_file = bz2.open(likelihood_filename, "wt")
#likelihood_file.open()

np.random.seed(123456789)

n_iter = 10000

# 1/2 bc we're using two metrics
P_star = (0.05) ** (1/2)
trajectories_count = 0
#lte =less than or equal
#gte = grater than or equal

times_to_exclude = {}






def get_likelihood(t_star, t, d_pm):
    t_star_index = np.where(t==t_star)[0][0]
    d_pm_gt = d_pm[t_star_index+1:]
    n_gt = len(d_pm_gt)
    d_pm_lte = d_pm[:t_star_index+1]
    n_lte = len(d_pm_lte)

    lambda_0_t_star = (1/n_lte) * sum(d_pm_lte)
    r_t_star = (1/lambda_0_t_star) * (1/n_gt) *  sum(d_pm_gt)
    sigma_lt_t_star = np.sqrt( (1/n_lte)* sum((d_pm_lte - lambda_0_t_star) **2) )
    sigma_gt_t_star = np.sqrt( (1/n_gt)* sum((d_pm_gt - (lambda_0_t_star*r_t_star)) **2) )
    l_t_star = -1*n_lte * np.log(sigma_lt_t_star) - n_gt*np.log(sigma_gt_t_star)
    return (t_star, n_gt, sigma_gt_t_star, n_lte, sigma_lt_t_star, r_t_star, l_t_star)


def get_test_statistics(t_pm, A_pm, D_pm):
    # autocorrelation
    f_pm = A_pm/D_pm
    f_bar = min((sum(A_pm) / sum(D_pm)), 0.5)
    autocorr_filter = np.where(abs(A_pm- f_bar*D_pm) < 1)[0]
    t_pm_auto = np.delete(t_pm, autocorr_filter)
    A_pm_auto = np.delete(A_pm, autocorr_filter)
    D_pm_auto = np.delete(D_pm, autocorr_filter)
    f_pm_auto = np.delete(f_pm, autocorr_filter)
    C_num_star = 0
    C_denom_star = 0
    for t_i in list(range(len(t_pm_auto)-1)):
        C_num_star += np.sqrt(D_pm_auto[t_i +1] * D_pm_auto[t_i]) * \
        (f_pm_auto[t_i +1] - f_bar) * (f_pm_auto[t_i] - f_bar) * \
        np.heaviside(abs(A_pm_auto[t_i +1] - f_bar*D_pm_auto[t_i +1])-1, 1/2) * \
        np.heaviside(abs(A_pm_auto[t_i] - f_bar*D_pm_auto[t_i])-1, 1/2)
        C_denom_star += np.sqrt(D_pm_auto[t_i +1] * D_pm_auto[t_i])

    C_star = C_num_star / ( (f_bar**2) * C_denom_star)

    # derived allele sojourn weight
    n_0 = len(np.where(A_pm == 0)[0])
    n = len(t_pm)
    if n_0 > 0.3*n:
        f_star = f_bar / (1+ np.exp( (n_0 - 0.3*n)/5 ))
    else:
        f_star = f_bar / (1+ np.exp(-1* (0.3*n - n_0)/5 ))

    t_pm_f_star = np.delete(t_pm, np.where(f_pm <= f_star)[0])
    f_pm_f_star = np.delete(f_pm, np.where(f_pm <= f_star)[0])
    if len(t_pm_f_star) < 2:
        return np.nan,np.nan
    if (len(t_pm_f_star) == 2) and ((t_pm_f_star[1] - t_pm_f_star[0]) > 500):
        return np.nan,np.nan

    I_list = []
    for win_size in list(range(2, len(t_pm_f_star)+1)):
        for win_start in  range(len(t_pm_f_star) - win_size +1):
            t1 = t_pm_f_star[win_start]
            t2 = t_pm_f_star[win_start + win_size-1]
            f_pm_window = np.asarray(f_pm_f_star[win_start: win_start + win_size])
            I_list.append( sum(f_pm_window - f_star) )
    I = max(I_list)
    # average frequency relaxation time
    #mean_freq = sum(A_pm) / sum(D_pm)
    #if n_0 > 0.3*n:
    #    T = 0
    #else:
    #    T = 0
    #    #A_pm
    #    for t_dash in range(2, len(t_pm)):
    #        partial_mean_freq = sum(A_pm[:t_dash+1]) / sum(D_pm[:t_dash+1])
    #        if (partial_mean_freq <= 0.6*mean_freq) and (t_dash > T):
    #            T = t_dash
    # dont return T
    return C_star, I


for depth in depth_file:
    depth_split = [x.strip() for x in depth.split(',')]
    D_pt_median = depth_split[-1].split(' ')
    D_pt_median = np.asarray([float(x) for x in D_pt_median])


def get_composite_p_value(t_pm, A_pm, D_pm):
    # moderate coverage timepoints
    #C_star_mut, I_mut, T_mut = get_test_statistics(t_pm, A_pm, D_pm)
    C_star_mut, I_mut = get_test_statistics(t_pm, A_pm, D_pm)
    if np.isnan(np.sum([C_star_mut, I_mut])):
        return None, None
    # error model for moderate-coverage timepoints
    f_mean = (sum(A_pm) / sum(D_pm))
    c = max(1, (1/2)*f_mean)
    f_pm = A_pm/D_pm
    C_star_mut_null = []
    I_mut_null = []
    #T_mut_null = []
    count = copy.copy(n_iter)
    while count != 0:
        p_hat_pm = np.random.permutation(f_pm) * c
        alpha_pm = np.random.poisson(D_pm*p_hat_pm)
        beta_pm = np.random.poisson(D_pm*(1-p_hat_pm))
        A_hat_pm = np.round( (alpha_pm/ (alpha_pm+beta_pm)) * D_pm )
        if np.count_nonzero(A_hat_pm) == 0:
            continue
        C_star_mut_i, I_mut_i = get_test_statistics(t_pm, A_hat_pm, D_pm)
        if np.isnan(np.sum([C_star_mut_i, I_mut_i])):
            continue
        C_star_mut_null.append(C_star_mut_i)
        I_mut_null.append(I_mut_i)
        #T_mut_null.append(T_mut_i)
        count -= 1

    P_C_star = len([x for x in C_star_mut_null if x > C_star_mut]) / n_iter
    P_I = len([x for x in I_mut_null if x > I_mut]) / n_iter
    #P_T = len([x for x in T_mut_null if x > T_mut]) / n_iter
    #if P_T == 0:
    #    P_T = 0.5
    #P_list = [P_C_star, P_I, P_T]
    P_list = [P_C_star, P_I]
    comp_stat = sum([np.heaviside(P_star-P_k, 1/2) * np.log(1/P_k) for P_k in P_list])
    comp_stat_null = []
    count2 = copy.copy(n_iter)
    while count2 != 0:
        p_hat_pm = np.random.permutation(f_pm) * c
        alpha_pm = np.random.poisson(D_pm*p_hat_pm)
        beta_pm = np.random.poisson(D_pm*(1-p_hat_pm))
        A_hat_pm = np.round( (alpha_pm/ (alpha_pm+beta_pm)) * D_pm )
        A_hat_pm[np.isnan(A_hat_pm)] = 0

        C_star_mut_i, I_mut_i = get_test_statistics(t_pm, A_hat_pm, D_pm)
        if np.isnan(np.sum([C_star_mut_i, I_mut_i])):
            continue
        # add pseudocount of 1
        P_C_star_i = (len([x for x in C_star_mut_null if x > C_star_mut_i]) +1) / (n_iter +1)
        P_I_i = (len([x for x in I_mut_null if x > I_mut_i]) +1) / (n_iter +1)
        #P_T_i = len([x for x in T_mut_null if x > T_mut_i]) / n_iter
        #if P_C_star_i == 0:
        #    P_C_star_i = 0.5
        #if P_T_i == 0:
        #    P_T_i = 0.5
        #if P_I_i == 0:
        #    P_I_i = 0.5
        #P_list_i = [P_C_star_i, P_I_i, P_T_i]
        P_list_i = [P_C_star_i, P_I_i]
        comp_stat_i = sum([np.heaviside(P_star-P_k, 1/2) * np.log(1/P_k) for P_k in P_list_i])
        comp_stat_null.append(comp_stat_i)
        count2 -= 1

    P_comp_stat = len([x for x in comp_stat_null if x > comp_stat]) / len(comp_stat_null)

    return C_star_mut, P_comp_stat




for snp in indel_file:
    trajectories_count += 1
    #if trajectories_count != 573:
    #    continue
    if (trajectories_count % 1000) == 0:
        print("Indel trajectory " + str(trajectories_count))
    snp_split = [x.strip() for x in snp.split(',')]
    t_pm = snp_split[3].split(' ')
    t_pm = np.asarray([int(x) for x in t_pm])
    A_pm = snp_split[4].split(' ')
    A_pm = np.asarray([int(x) for x in A_pm])
    D_pm = snp_split[5].split(' ')
    D_pm = np.asarray([int(x) for x in D_pm])
    count_At_Atplus1 = 0
    #print(A_pm)
    #print(D_pm)
    for i in range(0,len(A_pm)-1):
        if (A_pm[i] != 0) and (A_pm[i+1] != 0):
            count_At_Atplus1 += 1
    if count_At_Atplus1 == 0:
        continue
    if len(t_pm) < 4:
        continue
    A_D_stack = np.stack((A_pm, D_pm, A_pm/D_pm), axis=-1 )
    # criterion 1
    A_D_stack_At_2 = [x for x in A_D_stack if (x[0] >= 2) ]
    if len(A_D_stack_At_2) < 2:
        continue
    if len([x for x in A_D_stack_At_2 if (x[1] >=10) and (x[2] >=0.05)]) <1:
        continue
    # criterion 2
    if (D_pm[0] < 10) and (A_pm[0]/D_pm[0] >0.1):
        continue
    # criteria 3
    if len([x for x in A_D_stack if (x[0] >= 2) and (x[1]>=5) ]) < 3:
        continue
    # criterion 4
    if len([x for x in A_D_stack if (x[0] >= 3) and (x[1]>=5) and (x[2]>=0.1) ]) < 1:
        continue
    # criterion 5
    if max(A_pm/D_pm) - min((sum(A_pm) / sum(D_pm)), 0.5) < 0.1:
        continue

    # remove D_pmt < 5
    #remove_D_pmt = np.asarray([x for x,y in enumerate(D_pm) if y < 5])
    #D_pt_median_copy = np.empty_like(D_pt_median)
    #D_pt_median_copy[:] = D_pt_median
    #if len(remove_D_pmt) > 0:
    #    t_pm = np.delete(t_pm, remove_D_pmt)
    #    A_pm = np.delete(A_pm, remove_D_pmt)
    #    D_pm = np.delete(D_pm, remove_D_pmt)
    #    D_pt_median_copy = np.delete(D_pt_median_copy, remove_D_pmt)
    # remove low coverage timepoints
    #d_pm = D_pm / D_pt_median_copy
    d_pm = D_pm / D_pt_median
    l_list = [get_likelihood(t_, t=t_pm, d_pm=d_pm) for t_ in t_pm[1:-2] ]
    #l_list = sorted(l_list, key=lambda x: x[-1])
    l_list_candidates = [ l_list_i[-1] for l_list_i in l_list if (l_list_i[-2] < 0.5)  and (np.isinf(l_list_i[-1])==False)  ]
    if len(l_list_candidates) != 0:
        delta_l = max(l_list_candidates)
        delta_l_null_list = []
        count_low_cov = copy.copy(n_iter)
        while count_low_cov != 0:
            d_pm_null = np.copy(d_pm)
            np.random.shuffle(d_pm_null)
            l_list_null = [get_likelihood(t_, t=t_pm, d_pm=d_pm_null) for t_ in t_pm[1:-2] ]
            delta_l_null = [ l_list_null_i[-1] for l_list_null_i in l_list_null if (l_list_null_i[-2] < 0.5)  and (np.isinf(l_list_null_i[-1])==False)  ]
            if len(delta_l_null) == 0:
                continue
                #delta_l_null_list.append(0)
            else:
                delta_l_null_list.append(max(delta_l_null))
                count_low_cov -= 1
        #max_l = l_list[-1]
        # r threshold of 0.5 too conservative
        P_delta_l = (len([x for x in delta_l_null_list if x > delta_l]) +1) / (n_iter +1)
        if P_delta_l < 0.05:
            print("Deletion detected: excluding site from rest of analysis")
            C_star_mut = P_comp_stat = None
        else:
            C_star_mut, P_comp_stat = get_composite_p_value(t_pm, A_pm, D_pm)
            print(C_star_mut)
            #if C_star_mut == None:
            #    continue

    else:
        C_star_mut, P_comp_stat = get_composite_p_value(t_pm, A_pm, D_pm)
        #if C_star_mut == None:
        #    continue

    snp_split.append( " ".join([str(C_star_mut), str(P_comp_stat)]) )
    snp_split_merge = ", ".join(snp_split)

    likelihood_file.write(snp_split_merge)
    likelihood_file.write("\n")





for snp in snp_file:
    trajectories_count += 1
    #if trajectories_count != 573:
    #    continue
    if (trajectories_count % 1000) == 0:
        print("SNP trajectory " + str(trajectories_count))
    snp_split = [x.strip() for x in snp.split(',')]
    t_pm = snp_split[3].split(' ')
    t_pm = np.asarray([int(x) for x in t_pm])
    A_pm = snp_split[4].split(' ')
    A_pm = np.asarray([int(x) for x in A_pm])
    D_pm = snp_split[5].split(' ')
    D_pm = np.asarray([int(x) for x in D_pm])
    count_At_Atplus1 = 0
    #print(A_pm)
    #print(D_pm)
    for i in range(0,len(A_pm)-1):
        if (A_pm[i] != 0) and (A_pm[i+1] != 0):
            count_At_Atplus1 += 1
    if count_At_Atplus1 == 0:
        continue
    if len(t_pm) < 4:
        continue
    A_D_stack = np.stack((A_pm, D_pm, A_pm/D_pm), axis=-1 )
    # criterion 1
    A_D_stack_At_2 = [x for x in A_D_stack if (x[0] >= 2) ]
    if len(A_D_stack_At_2) < 2:
        continue
    if len([x for x in A_D_stack_At_2 if (x[1] >=10) and (x[2] >=0.05)]) <1:
        continue
    # criterion 2
    if (D_pm[0] < 10) and (A_pm[0]/D_pm[0] >0.1):
        continue
    # criteria 3
    if len([x for x in A_D_stack if (x[0] >= 2) and (x[1]>=5) ]) < 3:
        continue
    # criterion 4
    if len([x for x in A_D_stack if (x[0] >= 3) and (x[1]>=5) and (x[2]>=0.1) ]) < 1:
        continue
    # criterion 5
    if max(A_pm/D_pm) - min((sum(A_pm) / sum(D_pm)), 0.5) < 0.1:
        continue

    C_star_mut, P_comp_stat = get_composite_p_value(t_pm, A_pm, D_pm)

    snp_split.append( " ".join([str(C_star_mut), str(P_comp_stat)]) )
    snp_split_merge = ", ".join(snp_split)

    likelihood_file.write(snp_split_merge)
    likelihood_file.write("\n")




depth_file.close()
snp_file.close()
indel_file.close()
likelihood_file.close()
