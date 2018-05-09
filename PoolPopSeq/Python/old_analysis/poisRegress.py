from __future__ import division
import cleanData as cd
import pandas as pd
import os
import scipy.stats as stats
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.cov_struct import (Exchangeable,
    Independence,Autoregressive)
from statsmodels.genmod.families import Poisson
import numpy as np
from statsmodels.formula.api import ols
from patsy.contrasts import ContrastMatrix
from statsmodels.stats import anova
import statsmodels.api as sm

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/data/")

def poisRegress():
    IN_file_mut = mydir + 'gene_by_sample/D/D100/sample_by_gene.txt'
    IN_file_genes = mydir + 'reference_assemblies_task2_table/D.txt'
    IN_mut = pd.read_csv(IN_file_mut, sep = '\t')
    # remove columns with all zeros
    #print IN_mut.loc[:, (IN_mut != 0).any(axis=0)]
    IN_genes = pd.read_csv(IN_file_genes, sep = ' ')
    IN_mut = IN_mut.T.reset_index()
    IN_mut.columns = IN_mut.iloc[0]
    IN_mut = IN_mut.reindex(IN_mut.index.drop(0))
    IN_mut.rename(columns={'Unnamed: 0':'LocusTag'}, inplace=True)
    IN_merged = pd.merge(IN_mut, IN_genes, how='inner',
        on=['LocusTag'])
    IN_merged['SizeLog10'] = np.log10(IN_merged['Size'])
    # turn the wide form data into long form
    value_vars_strain = [x for x in IN_merged.columns if 'frequency_' in x]
    IN_merged_long = pd.melt(IN_merged, id_vars=['LocusTag', 'SizeLog10', \
        'GC', 'Sequence'], value_vars=value_vars_strain, \
        var_name = 'Line', value_name = 'Muts')
    IN_merged_long['TransferTime'] = IN_merged_long['Line'].apply(getTransferTime)
    IN_merged_long['Strain'] = IN_merged_long['Line'].apply(lambda x: x[-2])
    IN_merged_long['Replicate'] = IN_merged_long['Line'].apply(lambda x: x[-1])
    # make sure mutations are read as integers
    IN_merged_long.Muts = IN_merged_long.Muts.astype(int)
    # add info for day .....
    #IN_merged_long['Day'] = IN_merged_long['Line'].apply(lambda x: x[-3:])
    #print [IN_merged_long.iloc[:,i].apply(type).value_counts() for i in range(IN_merged_long.shape[1])]
    IN_merged_long = IN_merged_long.sort_values(['TransferTime', 'Replicate'], \
        ascending=[True, True])
    cov_matrix = IN_merged_long[['TransferTime', 'Replicate']]
    cov_matrix['TransferTime'] = cov_matrix['TransferTime'].apply(int)
    cov_matrix['Replicate'] = cov_matrix['Replicate'].apply(int)
    cov_matrix['TransferTime'] = np.log10(cov_matrix['TransferTime'])
    cov_matrix['Replicate'] = cov_matrix['Replicate'] - 1
    cov_matrix['TransferTime'] = cov_matrix['TransferTime'].apply(int)
    print IN_merged_long
    IN_merged_long.to_csv(mydir + 'gene_by_sample_long/gene_by_sample_long_D.txt', sep = '\t')
    #P = sm.families.Poisson()
    #NBv = sm.families.NegativeBinomial()
    #ex = sm.cov_struct.Exchangeable()
    # initial model considers only the effect of
    #model1 = sm.GEE.from_formula("Muts ~ GC + SizeLog10", "Line",
    #                   data=IN_merged_long, family=family, cov_struct=ex)
    #result1 = model1.fit()
    #print result1.summary()
    #ne = sm.cov_struct.Nested()
    #model = sm.GEE.from_formula("Muts ~ GC * SizeLog10", groups = "Line",
    #                        data=IN_merged_long, family=P, cov_struct=ne,
    #                        dep_data=cov_matrix)

    #model2 = sm.GEE.from_formula("Muts ~ GC * SizeLog10", groups = "Line",
    #                        data=IN_merged_long, family=P2, cov_struct=ne2,
    #                        dep_data=cov_matrix)
    #result = model.fit() #maxiter=10)
    #result2 = model2.fit()
    #print IN_merged_long
    #print cov_matrix
    #print result.summary()
    #print dir(result)
    #print anova.anova_lm(args = [model, model2])
    #print ne.summary()
    #print result.centered_resid



def test_multilevel2():
    data = pd.read_csv(mydir + "a-level-chemistry.txt", header=None, sep=' ')
    data.columns = ["Board", "A-Score", "Gtot", "Gnum",
                    "Gender", "Age", "Inst_Type", "LEA",
                    "Institute", "Student"]
    data["GCSE"] = data["Gtot"] / data["Gnum"]
    data["School"] = [str(x) + str(y) for x,y in
                  zip(data["LEA"], data["Institute"])]
    us = set(data["School"])
    us = {x: k for k,x in enumerate(list(us))}
    data["School"] = [us[x] for x in data["School"]]
    qdata = data[["GCSE", "Gender", "Age", "LEA", "Institute", "School"]]
    qdata = np.asarray(qdata)
    family = sm.families.Gaussian()
    ex = sm.cov_struct.Exchangeable()
    model1 = sm.GEE.from_formula("GCSE ~ Age + Gender", "LEA",
                           data=data, family=family, cov_struct=ex)
    result1 = model1.fit()
    ne = sm.cov_struct.Nested()
    model2 = sm.GEE.from_formula("GCSE ~ Age + Gender", "LEA",
                            data=data, family=family, cov_struct=ne,
                            dep_data=data["Institute"])
    result2 = model2.fit() #maxiter=10)
    print result2.fittedvalues()

def test_multilevel3():
    # Standard_deviation of errors that are independent among individuals
    obs_sd = 3
    # The standard deviation of errors that are shared within subgroups
    subgroup_sd = 2
    # The standard deviation of errors that are shared within (top-level) groups
    group_sd = 1
    # The number of groups
    n_group = 100
    # The number of subgroups in each group
    n_subgroup = 4
    # The number of observations in each subgroup
    n_obs = 5
    def generate_data():
        # exog data is iid standard normal
        exog = [np.random.normal(size=(n_obs, 2)) for i in range(n_group * n_subgroup)]
        # Storage
        endog = []
        groups = []
        nests = []
        # Generate the endog values
        ii = 0
        for i in range(n_group):
            # Group effect, shared by all observations within a group
            group_e = group_sd * np.random.normal()
            for j in range(n_subgroup):
                # Subgroup effect, shared by all observations within a subgroup
                subgroup_e = subgroup_sd * np.random.normal()
                # All regression slopes are equal to 1.
                expval = exog[ii].sum(1)
                # The total random error for one observation
                obs_e = obs_sd * np.random.normal(size=n_obs)
                # The endog value for one observation
                y = expval + group_e + subgroup_e + obs_e
                endog.append(y)
                groups.append(i*np.ones(n_obs))
                nests.append(j*np.ones(n_obs))
                ii += 1

        exog = np.concatenate(exog)
        endog = np.concatenate(endog)
        groups = np.concatenate(groups)
        nests = np.concatenate(nests)

        return endog, exog, groups, nests
    endog, exog, groups, nests = generate_data()
    #print endog, exog, groups,nests
    family = sm.families.Gaussian()
    ne = sm.cov_struct.Nested()

    mod1 = sm.GEE(endog, exog, groups=groups, family=family, cov_struct=ne, dep_data=nests)
    rslt1 = mod1.fit()
    print dir(rslt1)
    #print groups
    #print rslt1.summary()
    #params.append(rslt1.params)
    #dep_params.append(ne.dep_params)
    #scales.append(rslt1.scale)
