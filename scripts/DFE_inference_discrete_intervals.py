
"""
Uses dadi to infer the demographic history using syn-SFS and fitdadi to infer the DFE of a nonsynonymous sfs.
For now trying to use the example data to compute the SFS. Later we will test different h coefficients

JCM 201907011

usage: python DFE_inference_discrete_intervals.py 0.5


in which 0.5 is the dominance_coefficient 
"""

import sys
import os
import logging
import time
import argparse
import warnings
import numpy
import dadi
import Selection
import scipy.stats.distributions
import scipy.integrate
import scipy.optimize
import pandas as pd


# dadi characterizes genotype fitnesses as:
# 1, 1 + 2sh, and 1+2s, where 1+2sh
h_coefficient = sys.argv[1]
times = 5 # iteractions across the different DFE models 
print("Using the dominance coefficient equal to")
print(h_coefficient)

#####################################
# Defining the functions used
#####################################

def gamma_dist(mgamma, alpha, beta):
	"""Define a gamma distribution.

	mgamma: float which describes the mean value of this gamma
		 distribution.
	alpha: shape parameter of the gamma distribution.
	beta: scale parameter of the gamma distribution.
	"""
	return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta)	

def two_epoch_sel(params, ns, pts, h=0.5):
	nu, T, gamma = params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)
	phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
	fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
	return fs

def discrete_not_complementary(mgamma, p1, p2, p3, p4): #, alpha, beta
	"""Define a discrete not complementary distribtuion.

	mgamma: float which describes the mean value of this gamma distribution.
	p1, p2, p3, p4: are the breakpoints of a discrete distribution

	five bins, (the discrete DFE) with breaks at s = [0, 10-5, 10-4, 10-3, 10-2, 1]
	# bins are discretized by scaling s coefficients to 2N_anc*s (ancestral 2N is 100000)
	"""
	two_Nanc=10000
	mgamma=-mgamma
	#assume anything with gamma<1e-4 is neutral
	# Assuming anything with s = 1e-5 is bin 1 2Ns = (10000*1e-5)
	if (two_Nanc*0 <= mgamma) and (mgamma < two_Nanc*1e-5):
		bin1 = p1/(0.1 - 0)
		return(bin1)

	# Assume anything with s => 1e-5 and s < 1e-4 is bin 2 (2Ns = 0.1 - 1)
	if (two_Nanc*1e-5 <= mgamma) and (mgamma < two_Nanc*1e-4):
		bin2 = p2/(1 - 0.1)
		return(bin2)

	# Assume anything with s [1e-4, < 1e-3) is bin 3 (2Ns = 1 - 10)
	if (two_Nanc*1e-4 <= mgamma) and (mgamma < two_Nanc*1e-3):
		bin3 = p3/(10 - 1)
		return(bin3)

	# Assume anything with s [1e-3, < 1e-2) is bin 4 (2Ns = 10 - 100)
	if (two_Nanc*1e-3 <= mgamma) and (mgamma < two_Nanc*1e-2):
		bin4 = p4/(100 - 10)
		return(bin4)
	else:
		return (1-(p1+p2+p3+p4))/(100-1000)


def discrete_not_complementary(mgamma, a1, a2, a3, a4):
    #made this range overlap with int_breaks for convenience
    Nanc=10000.
    x0 = 0. * 2 * Nanc
    x1 = 1e-5 * 2 * Nanc
    x2 = 1e-4 * 2 * Nanc
    x3 = 1e-3 * 2 * Nanc
    x4 = 1e-2 * 2 * Nanc
    x5 = 1 * 2 * Nanc
    a5 = 1. - a1 - a2 - a3 - a4
    mgamma = -mgamma
    if (mgamma >= x0) and (mgamma <= x1):
        return (a1/(x1-x0))
    elif (mgamma > x1) and (mgamma <= x2):
        return (a2/(x2-x1))
    elif (mgamma > x2) and (mgamma <= x3):
        return (a3/(x3-x2))
    elif (mgamma > x3) and (mgamma <= x4):
        return (a4/(x4-x3))
    elif (mgamma > x4) and (mgamma <= x5):
        return (a5/(x5-x4))
    else:
        return 0

def discretedfe(mgamma, p1, p2, p3, p4, p5): #, alpha, beta
	"""Define a discrete uniform distribution with 5 discrete bins, in which p_i (i:1-5) are the variables for each selection coefficient bin.

	five bins, (the discrete DFE) with breaks at s = [0, 10-5, 10-4, 10-3, 10-2, 1]
	# bins are discretized by scaling s coefficients to 2N_anc*s (ancestral 2N is 100000)
	"""
	mgamma = -mgamma
	two_Nanc=10000
	# Assuming anything with s = 1e-5 is bin 1 2Ns = (10000*1e-5)
	if (two_Nanc*0 <= mgamma) and (mgamma < two_Nanc*1e-5):
		bin1 = p1/(two_Nanc*1e-5 - 0)
		return(bin1)

	# Assume anything with s => 1e-5 and s < 1e-4 is bin 2 (2Ns = 0.1 - 1)
	if (two_Nanc*1e-5 <= mgamma) and (mgamma < two_Nanc*1e-4):
		bin2 = p2/(two_Nanc*1e-4 - two_Nanc*1e-5)
		return(bin2)

	# Assume anything with s [1e-4, < 1e-3) is bin 3 (2Ns = 1 - 10)
	if (two_Nanc*1e-4 <= mgamma) and (mgamma < two_Nanc*1e-3):
		bin3 = p3/(two_Nanc*1e-3 - two_Nanc*1e-4)
		return(bin3)

	# Assume anything with s [1e-3, < 1e-2) is bin 4 (2Ns = 10 - 100)
	if (two_Nanc*1e-3 <= mgamma) and (mgamma < two_Nanc*1e-2):
		bin4 = p4/(two_Nanc*1e-2 - two_Nanc*1e-3)
		return(bin4)

	if (two_Nanc*1e-2 <= mgamma) and (mgamma < two_Nanc*1e-1):
		bin5 = p5/(two_Nanc*1e-1 - two_Nanc*1e-2)
		return(bin5)

	if (mgamma > 100):
		return(0)

# define a constraint function. At MLE, consfunc = 0, for which all the 4 discrete parameters should sum up to 1
def consfunc(x, *args):
	return 1-sum(x) #[0:-2]

def integrating_over_multiple_gammas(demog_params = [2, 0.05], 
									theta_ns = 4000., 
									ns = numpy.array([250]), 
									pts_l = [600, 800, 1000], 
									tiny_e = 0.000000000001,
									int_breaks = [1e-5, 0.1, 1, 10, 100, 500]):
	# Npts log-spaced points
	# -1e-9 is just adding tiny values so it is not exactly in the border of the break.
	int_breaks = [i - tiny_e for i in int_breaks]
	print(int_breaks)

	spectra = Selection.spectra(demog_params, ns, two_epoch_sel, pts_l=pts_l,
								int_breaks=int_breaks, Npts=200,
								echo=True, mp=True)
	return(spectra)


########################################################
# Integrating over a range of gammas and set breaks
########################################################

#set demographic parameters and theta. this is usually inferred from the synonymous sites, we will do it when not using simulated data
# demog_params = [2, 0.05]

theta_ns = 4000. # NS = 4000 = 4NALNS.
spectra = integrating_over_multiple_gammas()


########################################################
# Loading sample data
########################################################

data = dadi.Spectrum.from_file('example.sfs')
# This sample data has the 2anc size = 10000
# A two fold expansion was simulated (10000 -> 20000)
# demographic parameters are: 
# Sample data was drawn from a gamma distribution with alpha = 0.18 and beta = 686.7

########################################################
# Fitting a gamma DFE to the data with a grid search not complementary function
########################################################
gamma_vec = numpy.frompyfunc(gamma_dist, 3, 1)
sel_params = [0.1, 1000] # alpha, beta
lower_bound = [1e-3, 1e-3] # lower bound of alpha, lower bound of beta
upper_bound = [1, 50000.] # upper bound of alpha, upper bound of beta

gamma_max_likelihoods = []
gamma_guesses = dict()
for i in range(times):
	p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound, upper_bound=upper_bound)
	print("Beginning optimization with guess, {0}.".format(p0))

	popt = Selection.optimize_log(p0, data, spectra.integrate, gamma_vec, 
								   theta_ns, lower_bound=lower_bound, 
								   upper_bound=upper_bound, verbose=len(sel_params), 
								   maxiter=50)
	
	gamma_max_likelihoods.append(popt[0])
	gamma_guesses[popt[0]] = popt

gamma_max_likelihoods.sort()

########################################################
# Outputting guesses 
########################################################

results_discrete = []
print(gamma_guesses.keys())
for i in range(len(gamma_guesses.keys())):
	best_popt_discrete = gamma_guesses[gamma_max_likelihoods[i]]
	results_discrete.append([best_popt_discrete[0], best_popt_discrete[1][0], best_popt_discrete[1][1]]) # likelihood, alpha, beta

# Writing results as csvs:
df_discrete = pd.DataFrame(results_discrete, columns =['Likelihood', 'alpha', 'beta'], dtype = float) #, 'alpha', 'beta'
df_discrete.to_csv('DFE_inference_gamma_{h}_.csv'.format(h=h_coefficient), index=False)

########################################################
# Fitting a discrete DFE to the data with a grid search not complementary function
########################################################

# Vectorizing the DFE
discrete_vec = numpy.frompyfunc(discrete_not_complementary, 5, 1)
print(discrete_vec)

# Selection parameters in terms of 2Nanc*s
sel_params = [0.1, 1, 10, 100] # 0.2, 1000. #bin1, # bin2, #bin3, #bin4 #alpha, #beta
lower_bound = [1e-2, 1e-2, 1e-2, 1e-2] #, 1e-2, 1e-2 #bin1, # bin2, #bin3, #bin4, #bin5, #alpha, #beta
upper_bound = [1, 1, 1, 1] #1, 50000.

discrete_max_likelihoods = []
discrete_guesses = dict()
for i in range(times):
	p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound, upper_bound=upper_bound)
	print("Beginning optimization with guess, {0}.".format(p0))
	popt = Selection.optimize_log(p0, data, spectra.integrate, discrete_vec, 
								   theta_ns, lower_bound=lower_bound, 
								   upper_bound=upper_bound, verbose=len(sel_params), 
								   maxiter=50)
	
	discrete_max_likelihoods.append(popt[0])
	discrete_guesses[popt[0]] = popt

discrete_max_likelihoods.sort()

########################################################
# Outputting guesses 
########################################################

results_discrete = []
print(discrete_guesses.keys())
for i in range(len(discrete_guesses.keys())):
	best_popt_discrete = discrete_guesses[discrete_max_likelihoods[i]]
	results_discrete.append([best_popt_discrete[0], best_popt_discrete[1][0], best_popt_discrete[1][1], best_popt_discrete[1][2], best_popt_discrete[1][3]]) #,  best_popt_discrete[1][4],  best_popt_discrete[1][5]]

# Writing results as csvs:
df_discrete = pd.DataFrame(results_discrete, columns =['Likelihood', 'bin1', 'bin2', 'bin3', 'bin4'], dtype = float) #, 'alpha', 'beta'
df_discrete.to_csv('DFE_inference_discrete_{h}_.csv'.format(h=h_coefficient), index=False)


########################################################
# Fitting a discrete DFE to the data with a grid search COMPLEMENTARY
########################################################

# Vectorizing the DFE
discrete_vec = numpy.frompyfunc(discretedfe, 6, 1)

# Selection parameters in terms of 2Nanc*s
sel_params = [0.1, 1, 10, 100, 1000] # 0.2, 1000. #bin1, # bin2, #bin3, #bin4 #bin5
lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2] #, 1e-2, 1e-2 #bin1, # bin2, #bin3, #bin4, #bin5
upper_bound = [1, 1, 1, 1, 1] #1, 50000.

discrete_max_likelihoods = []
discrete_guesses = dict()
for i in range(times):
	p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound, upper_bound=upper_bound)
	print("Beginning optimization with guess, {0}.".format(p0))
	popt = Selection.optimize_cons(p0, data, spectra.integrate, discrete_vec, 
								   theta_ns, lower_bound=lower_bound, 
								   upper_bound=upper_bound, verbose=len(sel_params), 
								   maxiter=100, constraint=consfunc)
	
	discrete_max_likelihoods.append(popt[0])
	discrete_guesses[popt[0]] = popt

discrete_max_likelihoods.sort()

########################################################
# Outputting guesses 
########################################################

results_discrete = []
print(discrete_guesses.keys())
for i in range(len(discrete_guesses.keys())):
	best_popt_discrete = discrete_guesses[discrete_max_likelihoods[i]]
	print(best_popt_discrete)
	results_discrete.append([best_popt_discrete[0], best_popt_discrete[1][0], best_popt_discrete[1][1], best_popt_discrete[1][2], best_popt_discrete[1][3], best_popt_discrete[1][4]]) #,  best_popt_discrete[1][4],  best_popt_discrete[1][5]]

# Writing results as csvs:
df_discrete = pd.DataFrame(results_discrete, columns =['Likelihood', 'bin1', 'bin2', 'bin3', 'bin4', 'bin5'], dtype = float) #, 'alpha', 'beta'
df_discrete.to_csv('DFE_inference_discrete_constrained_{h}_.csv'.format(h=h_coefficient), index=False)
