#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:39:50 2022

@author: evanqu
"""
import os
import time
import numpy as np
import pandas as pd
import pickle
import gzip
import warnings
from scipy import stats
from scipy.optimize import minimize 
import subprocess
from statsmodels.base.model import GenericLikelihoodModel
# import sys
# sys.set_int_max_str_digits(0) 

#%%
# os.chdir('/Users/evanqu/Dropbox (MIT)/Lieberman Lab/Personal lab notebooks/Evan/1-Projects/phlame_project/results/2022_06_skin_global/Cacnes_2023')

# path_to_cts_file='5-counts/FDT01175cheek_ref_Pacnes_C1.counts.pickle.gz'
# path_to_classifier='Cacnes_classifiers/Cacnes_ALL_phylogrouplevel.classifier'
# level_input='Cacnes_classifiers/Cacnes_ALL_phylogroup_IDs.txt'

# path_to_output_frequencies='test_frequencies.csv'
# path_to_output_data='test_fitinfo.data'

# import time

# start = time.time()
# results = Classify(path_to_cts_file,
#                     path_to_classifier,
#                     level_input,
#                     path_to_output_frequencies,
#                     path_to_output_data,
#                     nchains = 1,
#                     nparams = 4)

# results.main(model_behavior=True)

# end = time.time()

# print(end-start)
# # ###############################################################################
# # cts = CountsMat(path_to_cts_file)

# cfr = read_phlame_classifier(path_to_classifier)

# mylevel = PhyloLevel(level_input, 
#                     cfr.clades, 
#                     cfr.clade_names)

# level_cfr = cfr.grab_level(mylevel)

#%%
def reclassify_from_data(path_to_data_file,
                         path_to_frequencies_file,
                         max_pi=0.3,
                         min_snps=10,
                         min_prob=0.75):

    if min_snps!=10:
        raise Exception('Changes in the min_snps threshold from default (10) will have to be re-MCMCed')      
    
    # Load in sample
    sample_frequencies = Frequencies(path_to_frequencies_file)
    data = FrequenciesData(path_to_data_file)
    
    # What clades have enough SNPs to be modeled in the original run
    model_bool = np.logical_or.reduce((data.counts_MLE != -1),1)

    # Make new_frequencies data structure to fill
    new_frequencies = np.full_like(sample_frequencies.freqs,0)
    
    # new_probs = np.full_like(data.prob,-1)
    
    # Recalc frequencies with new thresholds
    for i in range(len(new_frequencies)):
        
        # Ignore the clades not to model
        if not model_bool[i]:
            continue
        
        # Calc new prob; move on if below threshold
        pi_chain = data.chain[i,:,0]
        if np.sum(pi_chain < max_pi)/len(pi_chain) < min_prob:
            continue
        
        # Only consider the parts of the chain below max_pi
        chain_bool = pi_chain < max_pi
        
        p_mean = np.mean(data.chain[i,chain_bool,1])
        n_mean = np.mean(data.chain[i,chain_bool,2])
        
        counts_mean = n_mean*(1-p_mean)/p_mean
        
        new_frequencies[i] = counts_mean/data.total_MLE[i][0]
            
    new_frequencies_df = pd.DataFrame(new_frequencies, 
                                      index=sample_frequencies.freqs.index,
                                      dtype=float)
    
    return new_frequencies_df
    
class Frequencies():
    '''
    Holds clade frequency information from a given sample
    '''
    
    def __init__(self, path_to_frequencies_file):
        
        self.freqs = pd.read_csv(path_to_frequencies_file,
                                       index_col=0)
                
class FrequenciesData():
    '''
    Holds csSNP counts and modeling information from a given sample.
    '''

    def __init__(self, path_to_data_file):
        
        with gzip.open(path_to_data_file, 'rb') as f:
            
            data_dct, fit_info_dct = pickle.load(f)
            
            # clade_counts structured as follows
            
            self.clade_counts = data_dct['clade_counts']
            self.clade_counts_pos = data_dct['clade_counts_pos']
            
            self.counts_MLE = fit_info_dct['counts_MLE']
            self.total_MLE = fit_info_dct['total_MLE']
            self.counts_MAP = fit_info_dct['counts_MAP']
            self.chain = fit_info_dct['chain']
            self.prob = fit_info_dct['prob']

    
class Classify:
    '''
    Main controller of the classify step.
    
    Args:
    path_to_cts_file (str): Path to PhLAMe counts file.\n
    
    path_to_classifier (str): Path to PhLAMe classifier file.
    level_input (str): Path to file defining the phylogenetic level.
    to type strains at. This file can either list clade names as 
    they appear in the PhLAMe classifier OR custom group genomes 
    into clades in a 2 column .tsv {genomeID, cladeID}.\n
    
    path_to_output_frequencies (str): Path to output frequencies file.\n
    
    path_to_output_data (str, optional): If True, outputs a data file 
    with counts and modeling information at the defined level.
    Defaults to False.\n
    
    max_perc_diff (TYPE, optional): Maximum . Defaults to 0.3.\n
    
    max_snp_diff (TYPE, optional): If True, . Defaults to False.\n
    
    min_snps (int, optional): Minimum number of csSNPs with >0 counts
    to call a lineage as present. Defaults to 10.\n
    
    '''
    
    def __init__(self, path_to_cts_file, path_to_classifier, level_input,
                 path_to_output_frequencies, path_to_output_data=False,
                 max_pi=0.3, max_snp_diff=False, min_snps=10, min_prob=0.5,
                 niter=100000, nburn=2000, nchains = 1, nparams=4):

        self.__input_cts_file = path_to_cts_file
        self.__classifier_file = path_to_classifier
        self.__levels_file = level_input
        self.__output_freqs_file = path_to_output_frequencies
        self.__output_data_file = path_to_output_data

        self.max_pi = max_pi
        self.min_snps = min_snps
        self.min_prob = min_prob
        self.max_snp_diff = max_snp_diff
        self.niter = niter
        self.nburn = nburn
        self.nchains = nchains
        self.nparams = nparams
        
        self.tmpdir = path_to_output_frequencies.rstrip('_frequencies.csv')
                
    def main(self, model_behavior):
        
        # =====================================================================
        #  Load Data
        # =====================================================================
        print("Reading in file(s)...")
        
        self.countsmat = CountsMat(self.__input_cts_file)
        
        self.classifier = read_phlame_classifier(self.__classifier_file)
        
        self.mylevel = PhyloLevel(self.__levels_file, 
                                  self.classifier.clades, 
                                  self.classifier.clade_names)
        
        # Grab just information for the specific level
        self.level_cfr = self.classifier.grab_level(self.mylevel)

        # =====================================================================
        #  Sort through counts data
        # =====================================================================
        print("Sorting counts information...")
        
        self.index_counts()
        
        self.get_allele_counts()

        # =====================================================================
        #  Calculate & Save frequencies
        # =====================================================================
        print("Modeling counts...")

        # Old dynamic bayesian
        # mle_threshold = 6
        # nsnps_threshold = 2000
        
        self.calc_frequencies(model_behavior)
            
        self.save_frequencies()
        
    def index_counts(self):
        '''
        Make structures to index through counts array.
        '''
        
        # Informative positions for this level
        # Note this will be > true pos because of pos with multiple alleles
        informative_bool = np.nonzero(self.level_cfr.csSNPs)[0]
        self.informative_pos = self.level_cfr.csSNP_pos[informative_bool]
        
        # Corresponding index on counts mat
        self.informative_pos_idx = np.arange(0,len(self.countsmat.pos))\
            [np.in1d(self.countsmat.pos,self.informative_pos)]\
                [informative_bool]

        # informative_pos_idx = np.arange(0,len(pos))[np.sum(csSNPs,1)>0]
        counts_CSS_bool = np.in1d(self.countsmat.pos,self.informative_pos)

        # Check that all positions are accounted for
        if np.sum(counts_CSS_bool) != len(np.unique(self.informative_pos)):
            raise Exception('Counts and classifier positions do not match. ',
                            'Check that reference genomes are same.')
    
    def get_allele_counts(self):
        
        # =====================================================================
        #     Grab just the allele info from counts
        # =====================================================================
        
        np.seterr(divide='ignore', invalid='ignore') #ignore divide by 0 warnings
        
        counts = self.countsmat.counts
        counts_idx = self.informative_pos_idx
        alleles = self.level_cfr.alleles
        
        #Initialize data structures
        cts_f = np.zeros(len(counts_idx))
        cts_r = np.zeros(len(counts_idx))
        tot_f = np.zeros(len(counts_idx))
        tot_r = np.zeros(len(counts_idx))
        
        for i, p in enumerate(self.informative_pos):
            
            # -1 +3 is to convert 01234 NATCG to index 
            cts_f[i] = counts[counts_idx[i],
                              int(alleles[i]-1)]
            
            cts_r[i] = counts[counts_idx[i],
                              int(alleles[i]+3)]
            
            tot_f[i] = np.sum(counts[counts_idx[i],:4])
            tot_r[i] = np.sum(counts[counts_idx[i],4:])
            
        self.allele_counts = tuple([cts_f,cts_r,tot_f,tot_r])
    
    def calc_frequencies(self, model_behavior):

        clade_names = self.level_cfr.clade_names
        
        frequencies = np.zeros(len(clade_names))
    
        save_cts = []
        save_cts_pos = []
        
        save_chain = []

        save_cts_mle=np.full((len(clade_names),2),-1,dtype=np.float32)
        save_cts_map=np.full((len(clade_names),2),-1,dtype=np.float32)
        save_total_mle=np.full((len(clade_names),2),-1,dtype=np.float32)
        save_prob=np.full((len(clade_names),1),-1,dtype=np.float32)
        
        # Reshape data to group by clade
        for c in range(len(clade_names)):
            
            byclade_cts = self.reshape_byclade(self.allele_counts,
                                               self.level_cfr.allele_cidx, c)
            
            byclade_cts_pos = self.informative_pos[np.where(self.level_cfr.allele_cidx == c)]
    
            cts2model = byclade_cts[0] + byclade_cts[1]
            total2model = byclade_cts[2] + byclade_cts[3]

            # To have a frequency, a call must have > min_SNP pos with 
            # nonzero counts in either fwd OR rev reads
            if np.count_nonzero(cts2model) >= self.min_snps:
                                
                print(f"Fit results for clade: {clade_names[c]}")
                
                cts_fit = countsCSS_JAGS(self.tmpdir,
                                        cts2model,
                                        total2model)
                
                # cts_fit = countsCSS(cts2model)
                                
                # total_fit = countsCSS(total2model)
                
                print(f"MLE fit: cts lambda={cts_fit.counts_mle[0]:.2f} cts pi={cts_fit.counts_mle[1]:.2f}")
                print(f"total lambda={cts_fit.total_mle[0]:.2f} total pi={cts_fit.total_mle[1]:.2f}")

                # Bayesian modeling if lambda is below 6 OR npts > 1000
                # if (cts_fit.mle[0] <= mle_threshold and len(cts2model) < nsnps_threshold) \
                #     or model_behavior:

                    
                # Bayesian modeling
                prob = cts_fit.fit(nchains = self.nchains,
                                   niter = self.niter,
                                   nburn = self.nburn,
                                   max_pi = self.max_pi)
                
                # save_cts_map[c] = cts_fit.map
                save_chain.append(cts_fit.chain)
                
                save_prob[c] = prob

                #     print("Bayesian modeling...")
                #     prob = cts_fit.bayesian_fit(max_pi=self.max_pi,
                #                                 niter=self.niter,
                #                                 nburn=self.nburn)
                    
                #     print(f"MAP fit: lambda={cts_fit.map[0]:.2f} pi={cts_fit.map[1]:.2f}")
                #     print(f"Prob={prob}")

                #     save_cts_map[c] = cts_fit.map
                #     save_chain[c] = cts_fit.chain
                #     save_prob[c] = prob
                
                # else:
                #     # Else just take the mle pi fit
                #     prob = 1 if cts_fit.mle[1] < self.max_pi else 0
                
                # Only count clades with good probability
                if prob > self.min_prob:
                    frequencies[c] = cts_fit.counts_mle[0]/cts_fit.total_mle[0]

                # Save fit information
                save_cts_mle[c] = cts_fit.counts_mle
                save_total_mle[c] = cts_fit.total_mle
                
            # Otherwise is zero
            else:
                print(f"clade: {clade_names[c]} does not have not enough csSNPs to model")
                save_chain.append({})

            # Save counts data
            save_cts.append( byclade_cts )
            save_cts_pos.append( byclade_cts_pos )
            
            frequencies_df = pd.DataFrame(frequencies, 
                                          index=self.level_cfr.clade_names)

        self.frequencies = frequencies_df
        self.data = {'clade_counts':save_cts,
                     'clade_counts_pos':save_cts_pos}
        self.fit_info = {'counts_MLE': save_cts_mle,
                         'total_MLE':save_total_mle,
                         'counts_MAP':save_cts_map,
                         'chain':save_chain,
                         'prob':save_prob}
    
    def save_frequencies(self):
        '''
        Write frequencies and data to files.
        '''        
        
        # Save frequencies
        self.frequencies.to_csv(self.__output_freqs_file, sep=',')
        
        # Save data and fit info
        if self.__output_data_file:
            
            with gzip.open(self.__output_data_file,'wb') as f:
                
                pickle.dump([self.data, self.fit_info],f)

    @staticmethod
    def reshape_byclade(allele_counts, clade_idxs, i):
        '''
        Grab only allele counts belonging to a certain clade (i)

        '''
        
        clade_cts_f = allele_counts[0][np.where(clade_idxs==i)]
        clade_cts_r = allele_counts[1][np.where(clade_idxs==i)]
        clade_tot_f = allele_counts[2][np.where(clade_idxs==i)]
        clade_tot_r = allele_counts[3][np.where(clade_idxs==i)]

        return tuple([clade_cts_f,
                      clade_cts_r,
                      clade_tot_f,
                      clade_tot_r])

class CountsMat():
    '''
    Hold data and methods for a counts matrix
    '''
    def __init__(self, path_to_cts_file):
        
        with gzip.open(path_to_cts_file,'rb') as f:
            counts, pos = pickle.load(f)
                
        self.counts = counts
        self.pos = pos

def read_phlame_classifier(path_to_classifier):

    with gzip.open(path_to_classifier, 'rb') as f:
        cssnp_dct = pickle.load(f)
        
        csSNPs = cssnp_dct['cssnps']
        csSNP_pos = cssnp_dct['cssnp_pos']
        clades = cssnp_dct['clades']
        clade_names = cssnp_dct['clade_names']
        
    return PhlameClassifier(csSNPs, csSNP_pos, clades, clade_names)
        
class PhlameClassifier():
    '''
    Holds data and methods for a single Phlame Classifier object
    '''
    def __init__(self,
                 csSNPs, csSNP_pos,
                 clades, clade_names):

        self.csSNPs = csSNPs
        self.csSNP_pos = csSNP_pos
        self.clades = clades
        self.clade_names = clade_names
        
        # Get allele information
        self.get_alleles()
            
    def grab_level(self, PhyloLevel):
        '''
        Grab just information for a specific level.
        '''
        
        idx=[]
        for clade in PhyloLevel.clade_names:
            idx.append(np.where(self.clade_names==clade)[0][0])
        
        level_csSNPs = self.csSNPs[:,idx]

        level_csSNP_pos = self.csSNP_pos[~np.all(level_csSNPs == 0, axis=1)]
        
        return PhlameClassifier(level_csSNPs[~np.all(level_csSNPs == 0, axis=1)],
                                level_csSNP_pos,
                                PhyloLevel.clades,
                                PhyloLevel.names)
    
    def get_alleles(self):
        '''
        Get 1D list of every allele and corresponding clade.
        '''
        self.alleles = self.csSNPs[np.nonzero(self.csSNPs)]
        # corresponding clade index
        self.allele_cidx = np.nonzero(self.csSNPs)[1]

# Outdated as of 8/11/23

# class countsCSS_NEW:
#     '''
#     Hold and model counts data covering a single set of cluster-specific SNPs.
#     '''
    
#     def __init__(self, path_to_outdir,
#                  counts, total_counts):
        
#         self.outdir = path_to_outdir

#         self.counts = counts
#         self.total_counts = total_counts
        
#         # Names of temp files for JAGS
#         self.model_filename = "model.bug"
#         self.data_filename = "data.txt"
#         self.run_filename = "jags.bash"

#         # Maximum Likelihood fit
#         self.counts_mle = self.zip_fit_mle(self.counts)
#         self.total_mle = self.zip_fit_mle(self.total_counts)
        
#         # Load in models
#         # self.models()

#     def fit(self, 
#             nchains = 1,
#             niter = 100000,
#             nburn = 5000,
#             max_pi = 0.3):
#         '''
#         Fit counts data to model using JAGS.
#         '''

#         # =====================================================================
#         #  Calc posterior over p for total_counts
#         # =====================================================================
#         scaleby=10

#         self.p_update = self._calc_p_posterior(self.total_counts, scaleby)
#         print(self.p_update)

#         # =====================================================================
#         #  Write JAGS files
#         # =====================================================================
#         self.make_zinb_model(p_hyper=self.p_update)
#         self.make_zinb_run(nchains, round(niter/nchains), nburn)
        
#         # =====================================================================
#         #  Run JAGS
#         # =====================================================================
#         subprocess.run(f"mkdir -p {self.outdir}", shell=True)

#         print("Writing JAGS files..")

#         self.write_jags()
        
#         print("Running JAGS..")
        
#         self.run_jags(nchains, niter)        
        
#         # =====================================================================
#         #  Calc summary statistics
#         # =====================================================================
        
#         # self.chain_arr = np.array((self.chain['alpha'],
#         #                             self.chain['p'],
#         #                             self.chain['pi']))

#         prob = np.sum(self.chain[:,0] < max_pi)/len(self.chain[:,0])
        
#         return prob

#     def run_jags(self, nchains, niter):
#         '''
#         Run JAGS on command line
#         '''
#         # print(f'bash -c "source activate base; jags {self.outdir}/{self.runfile}"')
#         subprocess.run(f'bash -c "source activate base; jags {self.outdir}/{self.run_filename}"', 
#                        shell=True)
        
#         self.params = {}
        
#         with open(f'{self.outdir}/jags_index.txt') as file:
            
#             for line in file:
#                 ls = line.rstrip('\n').split(' ')
#                 # params[ls[0]] = ls[1:]
#                 self.params[ls[0]] = [int(s) for s in ls[1:]]
                
        
#         chain_ls = []
        
#         for chain_num in range(nchains):
            
#             jags_chain = np.loadtxt(f'{self.outdir}/jags_chain{chain_num+1}.txt')[:,1]
            
#             chain_tmp = np.zeros(round(niter))
            
#             for pidx, param in enumerate(self.params.keys()):
                
#                 chain_idxs = self.params[param]
                
#                 #-1 to convert 1-index to 0-index
#                 chain_tmp = jags_chain[chain_idxs[0]-1:chain_idxs[1]]
                
#                 chain_ls.append(chain_tmp)
        
#         self.chain = np.vstack(chain_ls).T

#     def write_jags(self):
#         '''
#         Write files to run JAGS on command line
#         '''
        
#         with open(f'{self.outdir}/{self.run_filename}','w') as f:

#             f.write(self.zinb_run)
        
#         with open(f'{self.outdir}/{self.model_filename}','w') as f:
#             f.write(self.zinb_model)
            
#         with open(f'{self.outdir}/{self.data_filename}','w') as f:
#             f.write(f'"y" <- c({",".join(self.counts.astype(str))})\n')
#             f.write(f'"N" <- {len(self.counts)}')
#             # f.write(f'"mean" <- {max(2,round(np.mean(counts)))}')

#     # def fit_mle(self, cts):
#     #     '''
#     #     Maximum likelihood fit to zero-inflated poisson
#     #     '''
#     #     # Inital guess
#     #     lambda_init = cts.mean()
#     #     excess_zeros = ((cts == 0).mean() - 
#     #                     stats.poisson.pmf(0, lambda_init) )
#     #     pi_init = excess_zeros if excess_zeros>0 else 0
        
#     #     # Fit parameters
#     #     result = minimize(self._zip_nloglike, 
#     #                       [lambda_init, pi_init], 
#     #                       args=cts)
#     #     return result.x
    
#     def zip_fit_mle(self, cts):
        
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             model = ZeroInflatedPoisson(cts)
#             fit = model.fit()
#             pi, lambda_ = fit.params
        
#         # return lambda_, pi
#             # CI_pi, CI_lambda_ = fit.conf_int()
#             # # range_CI_pi = CI_pi[1] - CI_pi[0]
#             # range_CI_lambda_ = CI_lambda_[1] - CI_lambda_[0]
        
#         return np.array([lambda_, pi])
        
#     def make_zinb_run(self,
#                       nchain,
#                       niter,
#                       nburn):
#         '''
#         JAGS format code to run MCMC on zero-inflated binomial model
#         '''
#         self.zinb_run = (
#             f'model in "{self.outdir}/{self.model_filename}"\n'
#             f'data in "{self.outdir}/{self.data_filename}"\n'
#             f'compile, nchains({nchain})\n'
#             f'initialize\n'
#             f'update {nburn}\n'
#             f'monitor pi\n'
#             f'monitor p\n'
#             f'monitor alpha\n'
#             f'monitor beta\n'
#             f'update {niter}\n'
#             f'coda *, stem("{self.outdir}/jags_")\n'
#             'exit')
    
#     def make_zinb_model(self, 
#                         alpha_hyper=(0.001,0.001), 
#                         pi_hyper=(1.001,1.001),
#                         p_hyper=(1.001,1.001)):
#         '''
#         JAGS format code to create a zero-inflated binomial model
#         Parameterized as follows:
#         Y ~ Gamma*Z
#         Z ~ Bernoulli(1-pi)
#         Gamma ~ Nbinom(alpha,p)
#         '''
        
#         self.zinb_model = ("""model {
#         ## Likelihood ##
#         for( i in 1:N ) {
#           y[i] ~ dpois( lambda.corrected[i] )
#           lambda.corrected[i] <- zero[i] * gamma[i] + 0.00001
#           zero[i] ~ dbern(1 - pi) # Bernoulli component for zero-inflation
#           gamma[i] ~ dgamma(alpha, beta)
#         }
#         beta <- p/(1 - p)
#         """
#         "## Priors ##\n"
#         f"    alpha ~ dgamma({alpha_hyper[0]},{alpha_hyper[1]})\n"
#         f"    pi ~ dbeta({pi_hyper[0]},{pi_hyper[1]}) # Beta prior on pi\n"
#         f"    p ~ dbeta({p_hyper[0]},{p_hyper[1]})\n"
#         "}"
#             )

#     def models(self):
#         '''
#         Old code to hold all the models
#         '''
#         self.zinb_model = """model {
#         ## Likelihood ##
#         for( i in 1:N ) {
#           y[i] ~ dpois( lambda.corrected[i] )
#           lambda.corrected[i] <- zero[i] * gamma[i] + 0.00001
#           zero[i] ~ dbern(1 - pi) # Bernoulli component for zero-inflation
#           gamma[i] ~ dgamma(alpha, beta)
#         }
#         beta <- p/(1 - p)

#         ## Priors ##
#         alpha ~ dgamma(0.001,0.001)
#         p ~ dbeta(1.001,1.001)
#         pi ~ dbeta(1.001,1.001) # Beta prior on pi        
#         }
#         """
        
#         self.nbinom_model = """model {
#         ## Likelihood ##
#         for( i in 1:N ) {
#           y[i] ~ dnegbin( p , n )
#         }
        
#         ## Priors ##
#         p ~ dbeta(1.001,1.001)
#         n ~ dgamma(0.001,0.001)
        
#         ## NB mean and variance
#         mu <- n*(1-p)/p
#         variance <- n*(1-p)/(p*p)
#         }
#         """
        
#         self.nbinom_run = (
#         f'model in "{self.outdir}/{self.modelfile}"\n'
#         f'data in "{self.outdir}/{self.datafile}"\n'
#         f'compile, nchains(1)\n'
#         f'initialize\n'
#         f'update 1000\n'
#         f'monitor p\n'
#         f'monitor n\n'
#         f'monitor mu\n'
#         f'monitor variance\n'
#         f'update 100000\n'
#         f'coda *, stem("{self.outdir}/jags_")\n'
#         'exit'
#         )
        
#         self.zinb_run = (
#         f'model in "{self.outdir}/{self.modelfile}"\n'
#         f'data in "{self.outdir}/{self.datafile}"\n'
#         f'compile, nchains(1)\n'
#         f'initialize\n'
#         f'update 5000\n'
#         f'monitor pi\n'
#         f'monitor p\n'
#         f'monitor alpha\n'
#         f'monitor beta\n'
#         f'update 100000\n'
#         f'coda *, stem("{self.outdir}/jags_")\n'
#         'exit'
#         )

#     @staticmethod
#     def _calc_p_posterior(total_counts, scaleby,
#                          b=1,
#                          c=1):
#         '''
#         Empirical posterior distribution over p for Nbinom-distributed counts.
#         Calculate with a beta(b,c) prior over p
#         '''
        
#         #If underdispersed relative to Poisson skip calc because it will break
#         if np.var(total_counts) <= np.mean(total_counts):
            
#             return (np.sum(total_counts),1)
            
#         k = np.mean(total_counts)**2/(np.var(total_counts)-np.mean(total_counts))
#         T = len(total_counts)
#         alpha_update = (k * T) + b
#         beta_update = np.sum(total_counts) + c
        
#         # Empirical max value of p
#         # p_max = np.mean(total_counts)/np.var(total_counts)
        
#         return (alpha_update/scaleby, beta_update/scaleby)

    
#     @staticmethod
#     def _zip_nloglike(params, counts):
#         '''
#         Negative log-likelihood function for zero-inflated Poisson
#         '''
        
#         lambda_ = params[0]
#         pi = params[1]
        
#         if lambda_ <= 0:
#             return np.inf
#         if pi < 0 or pi > 1:
#             return np.inf
    
#         Y=len(counts) - np.count_nonzero(counts)
#         n=len(counts)
        
#         return -( Y*np.log(pi + (1 - pi)*np.exp(-lambda_)) +
#                  (n - Y) * np.log(1 - pi) - 
#                  (n - Y) * lambda_ + 
#                  n * np.mean(counts) * np.log(lambda_) -
#                  np.sum(np.log(np.array([np.math.factorial(int(c)) for c in counts], dtype=float))) )
#                  #np.log(np.product(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) ) 

class countsCSS_JAGS:
    '''
    Hold and model counts data covering a single set of cluster-specific SNPs.
    '''
    
    def __init__(self, path_to_outdir,
                 counts, total_counts):
        
        self.outdir = path_to_outdir

        self.counts = counts
        self.total_counts = total_counts
        
        # Names of temp files for JAGS
        self.model_filename = "model.bug" # Model specification
        self.init_filename = "initial.txt" # Initial values and seed
        self.data_filename = "data.txt" # Counts data in JAGS format
        self.run_filename = "jags.bash" # 
        
        # self.measured_alpha = measured_alpha
        # Maximum Likelihood fit
        self.counts_mle = self.zip_fit_mle(self.counts)
        self.total_mle = self.zip_fit_mle(self.total_counts)
        
        # Load in models
        # self.models()

    def fit(self, 
            nchains = 1,
            niter = 100000,
            nburn = 5000,
            max_pi = 0.3):
        '''
        Fit counts data to model using JAGS.
        '''

        # =====================================================================
        #  Calc posterior over p for total_counts (This is bad? 8/9/23)
        # =====================================================================
        # scaleby=10

        # self.p_update = self._calc_p_posterior(self.total_counts, scaleby)
        # print(self.p_update)

        self.measured_alpha = max((1e-6,
                             (np.var(self.total_counts)-np.mean(self.total_counts))/np.mean(self.total_counts)**2))
        
        # =====================================================================
        #  Write JAGS files
        # =====================================================================
        print("Writing JAGS files..")

        # self.make_zinb_model()
        self.make_zinb_model(alpha_hyper = ((1/self.measured_alpha),1))
        self.make_inits()

        self.make_zinb_run(nchains, round(niter/nchains), nburn)
        
        self.write_jags()

        # =====================================================================
        #  Run JAGS
        # =====================================================================

        print("Running JAGS..")
        
        self.run_jags(nchains, niter)        
        
        # =====================================================================
        #  Calc summary statistics
        # =====================================================================

        # self.chain_arr = np.array((self.chain['alpha'],
        #                             self.chain['p'],
        #                             self.chain['pi']))

        prob = np.sum(self.chain['pi'] < max_pi)/len(self.chain['pi'])
        
        return prob

    def run_jags(self, nchains, niter):
        '''
        Run JAGS on command line
        '''
        
        # Run JAGS on command line
        # print(f'bash -c "source activate base; jags {self.outdir}/{self.runfile}"')
        subprocess.run(f'bash -c "source activate base; jags {self.outdir}/{self.run_filename}"', 
                       shell=True)
        
        time.sleep(5)
        
        # Get order of parameters
        self.params = {}
        with open(f'{self.outdir}/jags_index.txt') as file:  
            for line in file:
                ls = line.rstrip('\n').split(' ')
                # params[ls[0]] = ls[1:]
                self.params[ls[0]] = [int(s) for s in ls[1:]]
                
        # Save chain to data object
        chain_ls = []
        self.chain = {}
        for chain_num in range(nchains):
            
            # Read in JAGS chain
            jags_chain = np.loadtxt(f'{self.outdir}/jags_chain{chain_num+1}.txt')[:,1]
            
            chain_tmp = np.zeros(round(niter))
            
            for pidx, param in enumerate(self.params.keys()):
                
                chain_idxs = self.params[param]
                #-1 to convert 1-index to 0-index
                # chain_tmp = jags_chain[chain_idxs[0]-1:chain_idxs[1]]
                
                self.chain[param] = jags_chain[chain_idxs[0]-1:chain_idxs[1]]
                # chain_ls.append(chain_tmp)
        
        # self.chain = np.vstack(chain_ls).T

    def write_jags(self):
        '''
        Write files to run JAGS on command line
        '''
        subprocess.run(f"mkdir -p {self.outdir}", shell=True)

        # Write run file
        with open(f'{self.outdir}/{self.run_filename}','w') as f:
            f.write(self.zinb_run)
        
        # Model specification file
        with open(f'{self.outdir}/{self.model_filename}','w') as f:
            f.write(self.zinb_model)
            
        # Data file
        with open(f'{self.outdir}/{self.data_filename}','w') as f:
            f.write(f'"y" <- c({",".join(self.counts.astype(str))})\n')
            f.write(f'"N" <- {len(self.counts)}')
            # f.write(f'"mean" <- {max(2,round(np.mean(counts)))}')
        
        #Inits and seed file
        with open(f'{self.outdir}/{self.init_filename}','w') as f:
            f.write(self.inits)
            
            
    def zip_fit_mle(self, cts):
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = ZeroInflatedPoisson(cts)
            fit = model.fit()
            pi, lambda_ = fit.params
        
        # return lambda_, pi
            # CI_pi, CI_lambda_ = fit.conf_int()
            # # range_CI_pi = CI_pi[1] - CI_pi[0]
            # range_CI_lambda_ = CI_lambda_[1] - CI_lambda_[0]
        
        return np.array([lambda_, pi])
        
    def make_zinb_run(self,
                      nchain,
                      niter,
                      nburn):
        '''
        JAGS format code to run MCMC on zero-inflated binomial model
        '''
        self.zinb_run = (
            f'model in "{self.outdir}/{self.model_filename}"\n'
            f'data in "{self.outdir}/{self.data_filename}"\n'
            f'compile, nchains({nchain})\n'
            f'inits in "{self.outdir}/{self.init_filename}"\n'
            f'initialize\n'
            f'update {nburn}\n'
            f'monitor pi\n'
            f'monitor a\n'
            f'monitor b\n'
            f'monitor p\n'
            f'update {niter}\n'
            f'coda *, stem("{self.outdir}/jags_")\n'
            'exit')

    def make_inits(self):
        
        # counts_i ~ poisson(lambda_)
        lambda_init = self.counts.mean()
        
        # lambda_i ~ gamma(alpha, beta)
        alpha_init = self.measured_alpha
        pi_init = (self.counts == 0).mean() - stats.poisson.pmf(0, lambda_init)
        pi_init = max((1e-6,pi_init))
        
        p_init = (alpha_init/(alpha_init+lambda_init))

        self.inits = (
            f'"pi" <- {pi_init}\n'
            f'"p" <- {p_init}\n'
            f'"a" <- {alpha_init}\n'
            f'".RNG.name" <- "base::Super-Duper"\n'
            f'".RNG.seed" <- 1')
        
    def make_zinb_model(self, 
                        alpha_hyper=(0.001,0.001), 
                        pi_hyper=(1.001,1.001),
                        p_hyper=(1.001,1.001)):
        '''
        JAGS format code to create a zero-inflated binomial model
        Parameterized as follows:
        Y ~ Gamma*Z
        Z ~ Bernoulli(1-pi)
        Gamma ~ Nbinom(alpha,p)
        '''
        
        self.zinb_model = ("""model {
        ## Likelihood ##
        for( i in 1:N ) {
          y[i] ~ dpois( lambda.corrected[i] )
          lambda.corrected[i] <- zero[i] * gamma[i] + 0.00001
          zero[i] ~ dbern(1 - pi) # Bernoulli component for zero-inflation
          gamma[i] ~ dgamma(a, b)
        }
        b <- p/(1 - p)
        """
        "## Priors ##\n"
        f"    a ~ dgamma({alpha_hyper[0]},{alpha_hyper[1]})\n"
        f"    pi ~ dbeta({pi_hyper[0]},{pi_hyper[1]}) # Beta prior on pi\n"
        f"    p ~ dbeta({p_hyper[0]},{p_hyper[1]})\n"
        "}"
            )

    @staticmethod
    def _calc_p_posterior(total_counts, scaleby,
                         b=1,
                         c=1):
        '''
        Empirical posterior distribution over p for Nbinom-distributed counts.
        Calculate with a beta(b,c) prior over p
        '''
        
        #If underdispersed relative to Poisson skip calc because it will break
        if np.var(total_counts) <= np.mean(total_counts):
            
            return (np.sum(total_counts),1)
            
        k = np.mean(total_counts)**2/(np.var(total_counts)-np.mean(total_counts))
        T = len(total_counts)
        alpha_update = (k * T) + b
        beta_update = np.sum(total_counts) + c
        
        # Empirical max value of p
        # p_max = np.mean(total_counts)/np.var(total_counts)
        
        return (alpha_update/scaleby, beta_update/scaleby)

    
    @staticmethod
    def _zip_nloglike(params, counts):
        '''
        Negative log-likelihood function for zero-inflated Poisson
        '''
        
        lambda_ = params[0]
        pi = params[1]
        
        if lambda_ <= 0:
            return np.inf
        if pi < 0 or pi > 1:
            return np.inf
    
        Y=len(counts) - np.count_nonzero(counts)
        n=len(counts)
        
        return -( Y*np.log(pi + (1 - pi)*np.exp(-lambda_)) +
                 (n - Y) * np.log(1 - pi) - 
                 (n - Y) * lambda_ + 
                 n * np.mean(counts) * np.log(lambda_) -
                 np.sum(np.log(np.array([np.math.factorial(int(c)) for c in counts], dtype=float))) )
                 #np.log(np.product(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) ) 



class ZeroInflatedPoisson(GenericLikelihoodModel):
    
    def __init__(self, endog, exog=None, **kwds):
        if exog is None:
            exog = np.zeros_like(endog)
            
        super(ZeroInflatedPoisson, self).__init__(endog, exog, **kwds)
    
    def nloglikeobs(self, params):
        pi = params[0]
        lambda_ = params[1]

        return -np.log(self._zip_pmf(self.endog, pi=pi, lambda_=lambda_))
    
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params is None:
            lambda_start = self.endog.mean()
            excess_zeros = (self.endog == 0).mean() - stats.poisson.pmf(0, lambda_start)
            pi_start = excess_zeros if excess_zeros>0 else 0
            start_params = np.array([pi_start, lambda_start])
            
        return super(ZeroInflatedPoisson, self).fit(start_params=start_params,
                                                    maxiter=maxiter, maxfun=maxfun, **kwds)
    
    @staticmethod
    def _zip_pmf(x, pi, lambda_):
        '''zero-inflated poisson function, pi is prob. of 0, lambda_ is the fit parameter'''
        if pi < 0 or pi > 1 or lambda_ <= 0:
            return np.zeros_like(x)
        else:
            return (x == 0) * pi + (1 - pi) * stats.poisson.pmf(x, lambda_)


# class countsCSS:
#     '''
#     Hold and model counts data covering a single set of cluster-specific SNPs.
#     '''
    
#     def __init__(self, counts):
        
#         self.counts = counts
        
#         #Maximum Likelihood fit
#         self.mle = self.mle_fit()
                
#     def bayesian_fit(self, max_pi=0.3, niter=100000, nburn=2000):
#         '''
#         Bayesian fit to zero-inflated poisson
#         '''
        
#         # Initialize chain
#         bayes_init = self.mle
        
#         # Run chain
#         chain, accept_rate, lnprob = self._MCMC(self.counts, bayes_init,
#                                                self._zip_loglike, 
#                                                self._gaussian_prop_fn, 
#                                                prop_fn_kwargs={'sigma':np.array([0.2,0.4])},
#                                                niter=niter)
        
#         # Burn off first part of chain
#         chain_burn = chain[nburn:]
        
#         # Bin chain to nearest 0.01
        
#         lambda_bins = np.arange(0,
#                                 max(chain_burn[:,0]) + 0.01,
#                                 0.01)
#         pi_bins = np.arange(0,1.01,0.01)

#         # # 'Bin' values go leftmost (e.g. 0.01->0.02 count get assigned to 0.01)
#         lambda_hist = np.histogram(chain_burn[:,0],
#                                     bins=lambda_bins)
#         pi_hist = np.histogram(chain_burn[:,1],
#                                     bins=pi_bins)
        
#         #Multiply bins & grab MAP value
#         jt_bincount = np.outer(lambda_hist[0], pi_hist[0])
        
#         MAP_idx = np.unravel_index(jt_bincount.argmax(), 
#                                    jt_bincount.shape)
#         MAP_lambda = lambda_bins[MAP_idx[0]]
#         MAP_pi = pi_bins[MAP_idx[1]]

#         # Probability of a set of parameters where pi < max_pi
#         prob = np.sum(chain_burn[:,1] < max_pi)/(len(chain)-nburn)
       
#         # Save outputs
#         self.map = tuple([MAP_lambda,MAP_pi])
#         self.prob = prob
#         self.chain = chain_burn
        
#         return prob
        
#     def mle_fit(self):
#         '''
#         Maximum likelihood fit to zero-inflated poisson
#         '''
#         # Inital guess
#         lambda_init = self.counts.mean()
#         excess_zeros = ((self.counts == 0).mean() - 
#                         stats.poisson.pmf(0, lambda_init) )
#         pi_init = excess_zeros if excess_zeros>0 else 0
        
#         # Fit parameters
#         result = minimize(self._zip_nloglike, 
#                           [lambda_init, pi_init], 
#                           args=self.counts)
#         return result.x
    
#     @staticmethod
#     def _MCMC(counts, x0,
#           lnprob_fn, prop_fn, 
#           prop_fn_kwargs={}, 
#           niter=100000):
#         '''Metropolis-Hastings MCMC sampler
        
#         Args:
#             counts (arr): Vector of counts data to model
#             x0 (arr): Initial vector of parameters.
#             lnprob_fn (fxn): Function to compute log-posterior probability.
#             prop_fn (fxn): Function to draw proposals.
#             prop_fn_kwargs (TYPE, optional): Arguments for proposal function.
#             niter (int, optional): Number of iterations to run chain.
        
#         Returns:
#             chain (arr): Vector of chain.
#             accept_rate (float): Acceptance rate.
#             lnprob (arr): Log-posterior chain.
        
#         '''
        
#         ndim = len(x0)
        
#         # Data structures
#         chain = np.zeros((niter, ndim))
#         lnprob = np.zeros(niter)
#         accept_rate = np.zeros(niter)
        
#         # Generate accept decisions all @ once
#         u = np.random.uniform(0,1,niter)
        
#         # Initialize first samples
#         x0[x0<=0]=0.0001 # change 0-valued/neg initials to some small number
#         x = np.log(x0)    # transform initial parameters (lognormal) to new (normal)
#         chain[0] = tuple(x)
#         lnprob0 = lnprob_fn(x0, counts)
#         lnprob[0] = lnprob0
        
#         # start loop
#         naccept = 0
#         for ii in range(1, niter):
            
#             # Counter
#             if ii%10000==0:
#                 print('.')
        
#             # Propose the following parameters
#             x_star, factor = prop_fn(tuple(x), **prop_fn_kwargs)
            
#             # Because we log-transformed we cannot use the normal factor
#             # q(xi|xi+1)/q(xi+1|xi) = xi+1/xi for lognormal distribution
#             factor = np.product(np.exp(x_star)/np.exp(x))
            
#             if np.exp(x_star[1]) >=1: #auto-reject if pi > 1
#                 chain[ii] = x
#                 lnprob[ii] = lnprob0
#                 accept_rate[ii] = naccept / ii
                
#                 continue
        
#             # Compute Hastings ratio
#             # Log-posterior prob.
#             lnprob_plus1 = lnprob_fn(np.exp(x_star), counts)
#             # Proposal prob.
#             H = np.exp(lnprob_plus1 - lnprob0) * factor 
#             # equivalent to H = exp(ln(prob_plus1/prob0))
        
#             # Accept/Reject step (update acceptance counter)
#             if u[ii-1] < H:
#                 x = x_star
#                 lnprob0 = lnprob_plus1
#                 naccept += 1
        
#             # update chain
#             chain[ii] = x
#             lnprob[ii] = lnprob0
#             accept_rate[ii] = naccept / ii
        
#         return np.exp(chain), accept_rate, lnprob

#     @staticmethod
#     def _zip_loglike(params, counts):
#         '''
#         Log-likelihood function for zero-inflated Poisson
#         '''
        
#         lambda_ = params[0]
#         pi = params[1]
        
#         if lambda_ <= 0:
#             return -np.inf
#         if pi < 0 or pi > 1:
#             return -np.inf
    
#         Y=len(counts) - np.count_nonzero(counts)
#         n=len(counts)
        
#         return ( Y*np.log(pi + (1 - pi)*np.exp(-lambda_)) +
#                (n - Y) * np.log(1 - pi) - 
#                (n - Y) * lambda_ + 
#                 n * np.mean(counts) * np.log(lambda_) -
#                 np.sum(np.log(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) )
#                 #np.log(np.product(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) ) 

#     @staticmethod
#     def _zip_nloglike(params, counts):
#         '''
#         Negative log-likelihood function for zero-inflated Poisson
#         '''
        
#         lambda_ = params[0]
#         pi = params[1]
        
#         if lambda_ <= 0:
#             return np.inf
#         if pi < 0 or pi > 1:
#             return np.inf
    
#         Y=len(counts) - np.count_nonzero(counts)
#         n=len(counts)
        
#         return -( Y*np.log(pi + (1 - pi)*np.exp(-lambda_)) +
#                  (n - Y) * np.log(1 - pi) - 
#                  (n - Y) * lambda_ + 
#                  n * np.mean(counts) * np.log(lambda_) -
#                  np.sum(np.log(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) )
#                  #np.log(np.product(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) ) 
    
#     @staticmethod
#     def _gaussian_prop_fn(x0, sigma=0.2):
#         """
#         Gaussian proposal distribution.
    
#         Propose new parameters based on Gaussian distribution with
#         mean at current position and standard deviation sigma.
    
#         Since the mean is the current position and the standard
#         deviation is fixed. This proposal is symmetric so the ratio
#         of proposal densities is 1.
        
#         :param x0: Parameter array
#         :param sigma:
#             Standard deviation of Gaussian distribution. Can be scalar
#             or vector of len(x0)
    
#         :returns: (new parameters, ratio of proposal densities)
#         """
    
#         # Propose new parameters based on gaussian
#         # (Every normal is a version of stdnormal stretched
#         # by sigma and translated by x)
#         if hasattr(x0, '__len__'):
#             x_star = x0 + np.random.randn(len(x0)) * sigma
            
#         else:
#             x_star = x0 + np.random.randn() * sigma
        
#         # proposal ratio factor is 1 since jump is symmetric
#         qxx = 1
    
#         return (x_star, qxx)
        
        
class PhyloLevel:
    '''
    Holds a set phylogenetic level in reference to a classifier.
    '''
    
    def __init__(self, levelin, allclades, allclade_names):
        
        # Import information about the tree
        # Note: ideally this will be replaced by an import of a classifier class
        self.__allclades = allclades
        self.__allclade_names = allclade_names
        
        # Import levels from file or string
        if os.path.isfile(levelin):
            print('Reading in levels file...')
            self.clades, self.clade_names, self.names = self._parse_file(levelin)
        
        else: 
            print('Searching for the following levels in classifier:')
            print(levelin)
            self.clades, self.clade_names, self.names = self._parse_str(levelin)
            
    def _parse_file(self, levelin_file, uncl='-1'):
        '''
        Parse a 2 column delimited levels file 
        ( e.g genome_name1, clade_ID1
              genome_name2, clade_ID2 )
        '''
        
        with open(levelin_file,'r') as file:
            firstline = file.readline()
       
        if len(firstline.strip().split('\t'))==2:
            dlim='\t'
        elif len(firstline.strip().split(','))==2:
            dlim=','


        clade_ids = np.loadtxt(levelin_file, delimiter=dlim, dtype=str)
        
        clades, names = self._reshape(clade_ids[:,0],
                                      clade_ids[:,1],
                                      '-1')
        
        # Check that level is valid within classifier
        clade_names = self._check_file(clades, names)
        
        # Check that levels are not direct ancestors or descendants
        # of each other.
        self._ancdesc(clades, clade_names)
        
        return clades, clade_names, names

    def _parse_str(self, levelin_str):
        '''
        Parse a comma-delimited list of clade names (e.g. C.1,C.2,C.3)
        '''
        clade_names = levelin_str.strip().split(',')
        
        # Check that level is valid within classifier
        clades = self._check_str(clade_names)
        names = clade_names
        
        # Check that levels are not direct ancestors or descendants
        # of each other.
        self._ancdesc(clades, clade_names)
        
        return clades, clade_names, names
    
    def _check_file(self, clades, names):
        '''
        Check that groupings specified by a file are valid clades in classifer.
        '''
        
        clade_names=[]
        allclades_ls = list(self.__allclades.values())
        allclade_names = np.array(list(self.__allclades.keys()))
        
        # Plus return the missing data (clade_names)
        for i, clade in enumerate(clades):
        
            # Check if grouping of genomes is actually a clade in classifier
            match_bool = [ set(clade) == set(genomes) for genomes in allclades_ls ]
            
            if np.sum(match_bool)==0:
                raise Exception(f"Error: Could not find a valid clade corresponding to {names[i]} in classifier.")
            
            # Ugly
            clade_names.append(str(allclade_names[match_bool][0]))
            
        return clade_names
    
    def _check_str(self, clade_names):
        '''
        Check if clade names specified are valid clades in classifier.
        '''
        clades = []
        
        for name in clade_names:

            if name not in self.__allclade_names:
                raise Exception(f"Error: {name} is not a valid clade in the classifier!")

            clades.append(self.__allclades[name])
            
        return clades
    
    @staticmethod
    def _ancdesc(clades, clade_names):
        '''
        Check that no two clades in a list are direct ancestors or descendants 
        of each other.
        '''

        for i, clade in enumerate(clades):
            
            # Get other clades to cp against
            cp = (clades[:i] + clades[i+1 :])
            
            # Are any genomes duplicated in 2 clades
            dup_bool = [len(set(clade).intersection(set(genomes)))>0 for genomes in cp]

            if np.sum(dup_bool) > 0:
                raise Exception(f"Error: {clade_names[i]} is either an ancestor or descendant of {clade_names[dup_bool.index(True)]}.",
                                "Note that within a level the same genome cannot be included in two clades.")

    @staticmethod
    def _reshape(names_long, clade_ids_long, uncl_marker):
        '''
        Reshape a two column 'long' levels array into list of arrays
        '''
        
        clades = []; names = []
        
        for c in np.unique(clade_ids_long):
            names.append(str(c))
            clades.append(names_long[np.in1d(clade_ids_long,c)])
        
        # Remove things designated as unclassifed
        if uncl_marker in names:
            # Get index of unclassified names
            idx = names.index(uncl_marker)
            # And remove
            del names[idx]; del clades[idx]
        
        else:
            print('Nothing was found as unclassified in clade IDs. Ignore if intentional!')
                
        return clades, names

#%% Test countsCSS class

# counts = clade_counts[7][0] + clade_counts[7][1]
# init = np.array([0,0.1])
# niter=20000
# nburn=2000

# def bayesian_fit(self, max_pi=0.3, niter=100000, nburn=2000):
#     '''
#     Bayesian fit to zero-inflated poisson.
#     '''
    
#     # Initialize chain
#     bayes_init = init
    
#     # Run chain
#     chain, accept_rate, lnprob = MCMC(counts, bayes_init,
#                                         zip_loglike, 
#                                         gaussian_prop_fn, 
#                                         prop_fn_kwargs={'sigma':np.array([0.2,0.4])},
#                                         niter=niter)
    
#     # Burn off first part of chain
#     chain_burn = chain[nburn:]
    
#     # Bin chain to nearest 0.01
#     lambda_bins = np.arange(0,
#                             max(chain_burn[:,0]) + 0.02,
#                             0.01)
#     pi_bins = np.arange(0,1.01,0.01)

#     # 'Bin' values go leftmost (e.g. 0.01->0.02 count get assigned to 0.01)
#     lambda_hist = np.histogram(chain_burn[:,0],
#                                 bins=lambda_bins)
#     pi_hist = np.histogram(chain_burn[:,1],
#                                 bins=pi_bins)

#     # lambda_bincount = np.bincount( np.digitize(chain_burn[:,0],
#     #                                             bins=lambda_bins) )
#     # pi_bincount = np.bincount( np.digitize(chain_burn[:,1],
#     #                                         bins=pi_bins) )
    
#     #Multiply bins & grab MAP value
#     jt_bincount = np.outer(lambda_hist[0], pi_hist[0])
#     # jt_bincount = np.outer(lambda_bincount, pi_bincount)

#     MAP_idx = np.unravel_index(jt_bincount.argmax(), 
#                                jt_bincount.shape)
    
#     MAP_lambda = lambda_bins[MAP_idx[0]]
#     MAP_pi = pi_bins[MAP_idx[1]]
#         # nothing is binned to pi=0 (want to fix)
    
#     # Probability of a set of parameters where pi < max_pi
#     prob = np.sum(chain_burn[:,1] < max_pi)/(len(chain)-nburn)
   
#     # Save outputs
#     self.map = tuple([MAP_lambda,MAP_pi])
#     self.prob = prob
    
#     return prob

# def zip_loglike(params, counts):
#     '''
#     Log-likelihood function for zero-inflated Poisson
#     '''
    
#     lambda_ = params[0]
#     pi = params[1]
    
#     if lambda_ <= 0:
#         return -np.inf
#     if pi < 0 or pi > 1:
#         return -np.inf

#     Y=len(counts) - np.count_nonzero(counts)
#     n=len(counts)
    
#     return ( Y*np.log(pi + (1 - pi)*np.exp(-lambda_)) +
#            (n - Y) * np.log(1 - pi) - 
#            (n - Y) * lambda_ + 
#             n * np.mean(counts) * np.log(lambda_) -
#             np.sum(np.log(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) )
#             #np.log(np.product(np.array([np.math.factorial(int(c)) for c in counts], dtype='float128'))) ) 

# def gaussian_prop_fn(x0, sigma=0.2):
#     """
#     Gaussian proposal distribution.

#     Propose new parameters based on Gaussian distribution with
#     mean at current position and standard deviation sigma.

#     Since the mean is the current position and the standard
#     deviation is fixed. This proposal is symmetric so the ratio
#     of proposal densities is 1.
    
#     :param x0: Parameter array
#     :param sigma:
#         Standard deviation of Gaussian distribution. Can be scalar
#         or vector of len(x0)

#     :returns: (new parameters, ratio of proposal densities)
#     """

#     # Propose new parameters based on gaussian
#     # (Every normal is a version of stdnormal stretched
#     # by sigma and translated by x)
#     if hasattr(x0, '__len__'):
#         x_star = x0 + np.random.randn(len(x0)) * sigma
        
#     else:
#         x_star = x0 + np.random.randn() * sigma
    
#     # proposal ratio factor is 1 since jump is symmetric
#     qxx = 1

#     return (x_star, qxx)

# def MCMC(counts, x0,
#       lnprob_fn, prop_fn, 
#       prop_fn_kwargs={}, 
#       niter=100000):
#     '''Metropolis-Hastings MCMC sampler
    
#     Args:
#         counts (arr): Vector of counts data to model
#         x0 (arr): Initial vector of parameters.
#         lnprob_fn (fxn): Function to compute log-posterior probability.
#         prop_fn (fxn): Function to draw proposals.
#         prop_fn_kwargs (TYPE, optional): Arguments for proposal function.
#         niter (int, optional): Number of iterations to run chain.
    
#     Returns:
#         chain (arr): Vector of chain.
#         accept_rate (float): Acceptance rate.
#         lnprob (arr): Log-posterior chain.
    
#     '''
    
#     ndim = len(x0)
    
#     # Data structures
#     chain = np.zeros((niter, ndim))
#     lnprob = np.zeros(niter)
#     accept_rate = np.zeros(niter)
    
#     # Generate accept decisions all @ once
#     u = np.random.uniform(0,1,niter)
    
#     # Initialize first samples
#     x0[x0<=0]=0.0001 # change 0-valued/neg initials to some small number
#     x = np.log(x0)    # transform initial parameters (lognormal) to new (normal)
#     chain[0] = tuple(x)
#     lnprob0 = lnprob_fn(x0, counts)
#     lnprob[0] = lnprob0
    
#     # start loop
#     naccept = 0
#     for ii in range(1, niter):
        
#         # Counter
#         if ii%10000==0:
#             print('.')
    
#         # Propose the following parameters
#         x_star, factor = prop_fn(tuple(x), **prop_fn_kwargs)
        
#         # Because we log-transformed we cannot use the normal factor
#         # q(xi|xi+1)/q(xi+1|xi) = xi+1/xi for lognormal distribution
#         factor = np.product(np.exp(x_star)/np.exp(x))
        
#         if np.exp(x_star[1]) >=1: #auto-reject if pi > 1
#             chain[ii] = x
#             lnprob[ii] = lnprob0
#             accept_rate[ii] = naccept / ii
            
#             continue
    
#         # Compute Hastings ratio
#         # Log-posterior prob.
#         lnprob_plus1 = lnprob_fn(np.exp(x_star), counts)
#         # Proposal prob.
#         H = np.exp(lnprob_plus1 - lnprob0) * factor 
#         # equivalent to H = exp(ln(prob_plus1/prob0))
    
#         # Accept/Reject step (update acceptance counter)
#         if u[ii-1] < H:
#             x = x_star
#             lnprob0 = lnprob_plus1
#             naccept += 1
    
#         # update chain
#         chain[ii] = x
#         lnprob[ii] = lnprob0
#         accept_rate[ii] = naccept / ii
    
#     return np.exp(chain), accept_rate, lnprob

# #%% OLD

# def classify(path_to_cts_file, path_to_classifier, level_input, 
#               path_to_output_frequencies, path_to_output_data=False,
#               max_perc_diff=0.3, max_snp_diff=False, min_snps=10):
#     '''Type strains in a metagenome sample given a PhLAMe classifier file.

#     Args:
#         path_to_cts_file (str): Path to PhLAMe counts file.
#         path_to_classifier (str): Path to PhLAMe classifier file.
#         level_input (str): Path to file defining the phylogenetic level
#         to type strains at. This file can either list clade names comma-delimited,
#         or group genomes into clades with labels in a 2 column .tsv.
#         path_to_output_frequencies (str): Path to output frequencies file (.csv).
#         path_to_output_data (str, optional): If True, outputs a data file with counts
#         and modeling information at the defined level. Defaults to False.
#         max_perc_diff (TYPE, optional): Maximum . Defaults to 0.3.
#         max_snp_diff (TYPE, optional): If True, . Defaults to False.
#         min_snps (int, optional): Minimum number of csSNPs with >0 counts to call
#                                   a lineage as present. Defaults to 10.

#     Returns:
#         cts_cov (TYPE): DESCRIPTION.
#         fit_info (TYPE): DESCRIPTION.
#         sample_frequencies (TYPE): DESCRIPTION.

#     '''
    
#     print("Reading in counts file...")
#     with gzip.open(path_to_cts_file,'rb') as f:
#         counts, pos = pickle.load(f)
                
#     print("Reading in classifier file...")
#     with gzip.open(path_to_classifier, 'rb') as f:
#         cssnp_dct = pickle.load(f)
        
#         all_csSNPs = cssnp_dct['cssnps']
#         csSNP_pos = cssnp_dct['cssnp_pos']
#         all_clades = cssnp_dct['clades']
#         all_clade_names = cssnp_dct['clade_names']
        
        
#     print("Reading in levels...")
#     mylevel = PhyloLevel(level_input, all_clades, all_clade_names)
    
#     #grab csSNPs for level clades
#     csSNPs = all_csSNPs[:,np.in1d(all_clade_names, mylevel.clade_names)]

#     print("Classifying...")
#     cts_cov, fit_info, frequencies = counts_CSS(counts, pos, 
#                                                 csSNPs, csSNP_pos,
#                                                 mylevel.names,
#                                                 min_snps=min_snps)
#     #Save output
#     frequencies.to_csv(path_to_output_frequencies, sep=',')
    
#     if path_to_output_data:
#         with gzip.open(path_to_output_data,'wb') as f:
#             pickle.dump([cts_cov, fit_info],f)
    
#     return cts_cov, fit_info, frequencies


# def counts_CSS(counts, pos, csSNPs, csSNP_pos, clade_names, 
#                 min_snps=10, max_pi=0.3, min_prob=0.85):
#     '''
#     Args:
#         counts (arr): Array of allele counts across informative positions (8xp).
#         pos (arr): Corresponding positions on the reference genome for counts (p x 1).
#         csSNPs (arr): p x c matrix giving csSNPs for each clade in typing scheme (01234=NATCG).
#         csSNP_pos (arr): Corresponding positions on the reference for csSNPs.
#         clade_names (list): List of clade names to label output frequencies file.
#         min_snps (int): Minimum number SNPs with non-zero counts in sample to count a clade as present.

#     Returns:
#         counts_coverage (TYPE): DESCRIPTION.
#         filter_counts_coverage (TYPE): DESCRIPTION.
#         fit_info (TYPE): DESCRIPTION.
#         sample_frequencies (TYPE): DESCRIPTION.

#     '''
    
#     # =========================================================================
#     #     Make some data structures to index through arrays & store results
#     # =========================================================================
    
#     # positions on reference for this level
#     # counts_CSS_pos = csSNP_pos[np.sum(csSNPs,1)>0]
#     counts_CSS_pos = csSNP_pos[np.nonzero(csSNPs)[0]]

#     # corresponding index on counts mat
#     # counts_CSS_idx = np.arange(0,len(pos))[np.sum(csSNPs,1)>0]
#     counts_CSS_idx = np.arange(0,len(pos))[np.nonzero(csSNPs)[0]]

#     # cssnp allele 01234=NATCG 
#     # (note will be larger than counts_CSS_pos because pos can have >1 allele)
#     alleles = csSNPs[np.nonzero(csSNPs)]
#     # corresponding clade index
#     clades = np.nonzero(csSNPs)[1]

#     # Check
#     counts_CSS_bool = np.in1d(pos,counts_CSS_pos)
#     if np.sum(counts_CSS_bool) != len(np.unique(counts_CSS_pos)):
#         raise Exception('Counts and classifier positions do not match. ',
#                         'Check that reference genomes are same.')
                
#     # =========================================================================
#     #     Grab just the allele info from counts
#     # =========================================================================
    
#     np.seterr(divide='ignore', invalid='ignore') #ignore divide by 0 warnings
    
#     cts_f = np.zeros(len(counts_CSS_pos))
#     cts_r = np.zeros(len(counts_CSS_pos))
#     tot_f = np.zeros(len(counts_CSS_pos))
#     tot_r = np.zeros(len(counts_CSS_pos))
    
#     for i, p in enumerate(counts_CSS_pos):
        
#         # -1 +3 is to convert 01234 NATCG to index 
#         cts_f[i] = counts[counts_CSS_idx[i],int(alleles[i]-1)]
#         cts_r[i] = counts[counts_CSS_idx[i],int(alleles[i]+3)]
        
#         tot_f[i] = np.sum(counts[counts_CSS_idx[i],:4])
#         tot_r[i] = np.sum(counts[counts_CSS_idx[i],4:])

#     # def get_allele_counts(self):
        
#     #     # =====================================================================
#     #     #     Grab just the allele info from counts
#     #     # =====================================================================
        
#     #     np.seterr(divide='ignore', invalid='ignore') #ignore divide by 0 warnings
        
#     #     counts = self.countsmat.counts
#     #     counts_idx = self.informative_pos_idx
#     #     alleles = self.level_cfr.alleles
        
#     #     #Initialize data structures
#     #     cts_f = np.zeros(len(counts_idx))
#     #     cts_r = np.zeros(len(counts_idx))
#     #     tot_f = np.zeros(len(counts_idx))
#     #     tot_r = np.zeros(len(counts_idx))
        
#     #     for i, p in enumerate(self.informative_pos):
            
#     #         # -1 +3 is to convert 01234 NATCG to index 
#     #         cts_f[i] = counts[counts_idx[i],
#     #                           int(alleles[i]-1)]
            
#     #         cts_r[i] = counts[counts_idx[i],
#     #                           int(alleles[i]+3)]
            
#     #         tot_f[i] = np.sum(counts[counts_idx[i],:4])
#     #         tot_r[i] = np.sum(counts[counts_idx[i],4:])
            
#     #     self.allele_counts = tuple([cts_f,cts_r,tot_f,tot_r])


#     # =========================================================================
#     #     For each clade - fit data to expected counts model
#     # =========================================================================
    
#     frequencies = np.zeros(len(clade_names))
#     save_cts_f = []; save_cts_r = []
#     save_tot_f = []; save_tot_r = []
#     save_cts_pos = []


#     # Reshape data to group by clade
#     for c in range(len(clade_names)):
        
#         # Grab fwd and rev counts for this clade
#         clade_cts_f = cts_f[np.where(clades==c)]
#         clade_cts_r = cts_r[np.where(clades==c)]
#         clade_tot_f = tot_f[np.where(clades==c)]
#         clade_tot_r = tot_r[np.where(clades==c)]
        
#         clade_cts_pos = counts_CSS_pos[np.where(clades==c)]

#         # To have a frequency, a call must have > min_SNP pos with 
#         # nonzero counts in either fwd OR rev reads
#         cts2model = clade_cts_f + clade_cts_r
#         total2model = clade_tot_f + clade_tot_r

#         if np.count_nonzero(cts2model) > min_snps:
            
#             print(f"Fit results for clade: {clade_names[c]}")
                        
#             # Zero-inflated Poisson model
                
#             cts_lambda, cts_pi, _ = zip_fit(cts2model)
#             tot_lambda, tot_pi, _ = zip_fit(total2model)
            
#             print(f"MLE FIT: lambda={cts_lambda:.2f}, pi={cts_pi:.2f}")
#             # Only count clades with good probability
#             if cts_pi < max_pi:
#                 frequencies[c] = cts_lambda/tot_lambda
        
#         # Otherwise is zero
#         else:
#             print(f"clade: {clade_names[c]} does not have not enough csSNPs to model")
        
#         # Save data and fit information
#         save_cts_f.append( clade_cts_f )
#         save_cts_r.append( clade_cts_r )
#         save_tot_f.append( clade_tot_f )
#         save_tot_r.append( clade_tot_r )
#         save_cts_pos.append( clade_cts_pos )

                
#     # Output everything as nice data structures
#     frequencies[np.isnan(frequencies)] = 0
#     sample_frequencies = pd.DataFrame(frequencies)
#     sample_frequencies.columns = clade_names

#     data = [save_cts_f, save_cts_r, 
#             save_tot_f, save_tot_r,
#             save_cts_pos]
    
#     fit_info = 'placeholder'
#     # fit_info = [cts_lambda_arr, cts_pi_arr,
#     #             total_lambda_arr, total_pi_arr]

#     return data, fit_info, sample_frequencies

# def zip_pmf(x, pi, lambda_):
#     '''zero-inflated poisson function, pi is prob. of 0, lambda_ is the fit parameter'''
#     if pi < 0 or pi > 1 or lambda_ <= 0:
#         return np.zeros_like(x)
#     else:
#         return (x == 0) * pi + (1 - pi) * stats.poisson.pmf(x, lambda_)

# class ZeroInflatedPoisson(GenericLikelihoodModel):
#     def __init__(self, endog, exog=None, **kwds):
#         if exog is None:
#             exog = np.zeros_like(endog)
            
#         super(ZeroInflatedPoisson, self).__init__(endog, exog, **kwds)
    
#     def nloglikeobs(self, params):
#         pi = params[0]
#         lambda_ = params[1]

#         return -np.log(zip_pmf(self.endog, pi=pi, lambda_=lambda_))
    
#     def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
#         if start_params is None:
#             lambda_start = self.endog.mean()
#             excess_zeros = (self.endog == 0).mean() - stats.poisson.pmf(0, lambda_start)
#             pi_start = excess_zeros if excess_zeros>0 else 0
#             start_params = np.array([pi_start, lambda_start])
            
#         return super(ZeroInflatedPoisson, self).fit(start_params=start_params,
#                                                     maxiter=maxiter, maxfun=maxfun, **kwds)

# def zip_fit(data):
    
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         model = ZeroInflatedPoisson(data)
#         fit = model.fit()
#         pi, lambda_ = fit.params
    
#     # return lambda_, pi
#         CI_pi, CI_lambda_ = fit.conf_int()
#         # range_CI_pi = CI_pi[1] - CI_pi[0]
#         range_CI_lambda_ = CI_lambda_[1] - CI_lambda_[0]
    
#     return lambda_, pi, range_CI_lambda_

# def plot_chain(counts, chain, 
#                accept_rate, lnprob,
#                nburn=2000):
    
#     # Burn initial samples    
#     chain_burn = chain[nburn:]
#     accept_rate_burn = accept_rate[nburn:]
#     lnprob_burn = lnprob[nburn:]

#     # Bin chain parameters
#     lambda_bins=np.arange(0,max(chain_burn[:,0])+0.01,0.01)
#     lambda_bincount = np.bincount(np.digitize(chain_burn[:,0], bins=lambda_bins))
    
#     pi_bins=np.arange(0,1,0.01)
#     pi_bincount = np.bincount(np.digitize(chain_burn[:,1], bins=pi_bins))

#     # PLOT
#     fig, axs = plt.subplots(4,2)
    
#     # 2D histogram of draws
#     x_bins=np.arange(min(chain_burn[:,0])-0.5,max(chain_burn[:,0])+0.5,0.01)
#     y_bins=np.arange(0,1,0.01)
#     axs[0,0].hist2d(chain_burn[:,0], chain_burn[:,1], bins = [x_bins, y_bins], cmap='inferno')
#     axs[0,0].set_xlabel('Lambda'); axs[0,0].set_ylabel('Pi')
#     axs[0,0].set_title('2D Histogram of draws')
    
#     # Histogram of counts
#     axs[0,1].hist(counts,bins=np.arange(0,max(counts)+2),color='k', alpha=0.5)
#     axs[0,1].set_title('Histogram of counts')
#     axs[0,1].legend()
    
#     # Lambda posterior estimate marginalized over pi
#     axs[1,0].plot(lambda_bins[:len(lambda_bincount)], lambda_bincount/np.sum(lambda_bincount), color='b', label='MCMC')
#     axs[1,0].set_xlabel('Lambda'); axs[1,1].set_ylabel('Probability')
#     axs[1,0].set_title('Lambda posterior pdf')
#     axs[1,0].axvline(chain[0,0], color='r', label='MLE Lambda estimate')
#     axs[1,0].axvline(np.mean(chain_burn[:,0]), color='k', label='MAP Lambda estimate')
#     axs[1,0].legend()

#     # Pi posterior estimate marginalized over lambda
#     axs[1,1].plot(np.arange(0,1,0.01)[:len(pi_bincount)], pi_bincount/np.sum(pi_bincount), color='g', label='MCMC')
#     axs[1,1].set_xlabel('Pi'); axs[1,1].set_ylabel('Probability')
#     axs[1,1].set_title('Pi posterior pdf')
#     axs[1,1].axvline(chain[0,1], color='r', label='MLE Pi estimate')
#     axs[1,1].axvline(np.mean(chain_burn[:,1]), color='k', label='MAP Pi estimate')
#     axs[1,1].legend()    

#     # Lambda trace plot
#     axs[2,0].plot(np.arange(0,nburn), chain[:nburn,0], linewidth=0.3, color='r')
#     axs[2,0].plot(np.arange(nburn,len(chain)), chain_burn[:,0], linewidth=0.3, color='b')
#     axs[2,0].set_xlabel('Iteration number'); axs[2,0].set_ylabel('Lambda')
#     # axs[2,0].set_xlim(2700,2900)
#     axs[2,0].set_title('Lambda trace plot')
#     # Pi trace plot
#     axs[2,1].plot(np.arange(0,nburn), chain[:nburn,1], linewidth=0.3, color='r')
#     axs[2,1].plot(np.arange(nburn,len(chain)), chain_burn[:,1],linewidth=0.3, color='g')
#     axs[2,1].set_xlabel('Iteration number'); axs[2,1].set_ylabel('Pi')
#     axs[2,1].set_title('Pi trace plot')
#     # axs[2,1].set_xlim(2700,2900)
    
#     # Acceptance Rate
#     axs[3,0].plot(np.arange(0,len(accept_rate)), accept_rate, color='k', linewidth=0.3)
#     axs[3,0].plot(np.arange(0,nburn), accept_rate[:nburn], color='r', linewidth=0.3)
#     axs[3,0].set_xlabel('Iteration number'); axs[3,0].set_ylabel('Acceptance rate')
#     axs[3,0].set_title('Acceptance Rate')
#     # Log-posterior
#     axs[3,1].plot(np.arange(0,len(lnprob)), lnprob, color='k', linewidth=0.3)
#     axs[3,1].plot(np.arange(0,nburn), lnprob[:nburn], color='r', linewidth=0.3)
#     axs[3,1].set_xlabel('Iteration number'); axs[3,1].set_ylabel('Log-posterior')
#     axs[3,1].set_title('Log-posterior plot')

#     return fig