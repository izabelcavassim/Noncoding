
"""
Uses dadi to infer the demographic histoy using syn-SFS and fitdadi to infer the DFE of a nonsynonymous sfs.

JCM 201907011
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

#####################################
# Defining the functions used
#####################################

def gamma_dist(mgamma, alpha, beta):
	"""Define a gamma distribution.

	self: reference to this instance of a gamma distribution.
	mgamma: float which describes the mean value of this gamma
		 distribution.
	alpha: shape parameter of the gamma distribution.
	beta: scale parameter of the gamma distribution.
	"""
	return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta)	

#the demography + selection function. single size change and selection.
def two_epoch_sel(params, ns, pts, h=0.5):
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def neugamma(mgamma, pneu, pgamma, alpha, beta):
    mgamma=-mgamma
    #assume anything with gamma<1e-4 is neutral
    if (0 <= mgamma) and (mgamma < 1e-4):
        return pneu/(1e-4) + pgamma*Selection.gamma_dist(-mgamma, alpha, beta)
    else:
        return Selection.gamma_dist(-mgamma, alpha, beta) * pgamma

def discretegamma(mgamma, p1, p2, p3, p4, p5, alpha, beta):
	"""Define a gamma distribution with 5 discrete bins, in which p_i (i:1-5) are the variables for each selection coefficient bin.

	five bins, (the discrete DFE) with breaks at s = [0, 10-5, 10-4, 10-3, 10-2, 1]

	self: reference to this instance of a discrete-gamma distribution.
	mgamma: float which describes the mean value for each bin of the discrete-gamma distribution
	alpha: shape parameter of the non-neutral elements of the neutral-gamma distribution.
	beta: scale parameter of the non-neutral elements of the neutral-gamma distribution.
	
	"""
	mgamma = -mgamma
	# Assume anything with gamma < 1e-4 is bin 1
	if (0 <= mgamma) and (mgamma < 1e-5):
		bin1 = (1e-5 - 0)/p1
		#bin1 =  p1 * Selection.gamma_dist(-mgamma, alpha, beta)
		return(bin1)

	# Assume anything with gamma [1e-5, < 1e-4) is bin 2
	if (1e-5 <= mgamma) and (mgamma < 1e-4):
		bin2 = (1e-4 - 1e-5)/p2
		#bin2 =  p2 * Selection.gamma_dist(-mgamma, alpha, beta)
		return(bin2)

	# Assume anything with gamma [1e-4, < 1e-3) is bin 3
	if (1e-4 <= mgamma) and (mgamma < 1e-3):
		bin3 = (1e-3 - 1e-4)/p3
		#bin3 =  p3 * Selection.gamma_dist(-mgamma, alpha, beta)
		return(bin3)

		# Assume anything with gamma [1e-3, < 1e-2) is bin 4
	if (1e-3 <= mgamma) and (mgamma < 1e-2):
		bin4 = (1e-2 - 1e-3)/p4
		#bin4 =  p4 * Selection.gamma_dist(-mgamma, alpha, beta)
		return(bin4)

	if (1e-2 <= mgamma) and (mgamma < 1):
		bin5 = (1-1e-2)/p5
		#bin5 =  p5 * Selection.gamma_dist(-mgamma, alpha, beta)
		return(bin5)

	if (mgamma > 1):
		return(0)
	#else:
	#	bin5 = Selection.gamma_dist(-mgamma, alpha, beta) * (1 -(p1+p2+p3+p4))

# define a constraint function. At MLE, consfunc = 0, for which all the 5 discrete gamma parameters should sum up to 1
def consfunc(x, *args):
    return 1-sum(x[0:-2])

#set demographic parameters and theta. this is usually inferred from
#synonymous sites
demog_params = [2, 0.05]
theta_ns = 4000.
ns = numpy.array([250])


#integrate over a range of gammas and set breaks
pts_l = [600, 800, 1000]
# -1e-9 is just adding tiny values so it is not exactly in the border of the break.
#int_breaks = [0, 0.00001-1e-9, 0.0001-1e-9, 0.001-1e-9, 0.01-1e-9, 1-1e-9
#int_breaks = [1e-5, 0.00001, 0.0001, 0.001, 0.01, 1, 500]
int_breaks = [1e-5-1e-9, 1e-4-1e-9, 1e-3-1e-9, 1e-2-1e-9, 0.1-1e-9]
#int_breaks = [1e-5, 1e-4, 1e-3, 1e-2, 0.1]
spectra = Selection.spectra(demog_params, ns, two_epoch_sel, pts_l=pts_l,
                            int_breaks=int_breaks, Npts=300,
                            echo=True, mp=True)

#load sample data
data = dadi.Spectrum.from_file('example.sfs')

########################################################
# Fitting a discrete gamma to the data with a grid search
########################################################
neugamma_vec = numpy.frompyfunc(discretegamma, 8, 1)

print(neugamma_vec)


sel_params = [1e-5-1e-9, 1e-4-1e-9, 1e-3-1e-9, 1e-2-1e-9, 0.1, 0.2, 1000.] #bin1, # bin2, #bin3, #bin4, #bin5, #alpha, #beta
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-2] #bin1, # bin2, #bin3, #bin4, #bin5, #alpha, #beta
upper_bound = [1, 1, 1, 1, 1, 1, 50000.]

p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound, upper_bound=upper_bound)
print("This is the initial values")
print(p0)
popt = Selection.optimize_cons(p0, data, spectra.integrate, neugamma_vec, 
                               theta_ns, lower_bound=lower_bound, 
                               upper_bound=upper_bound, verbose=len(sel_params), 
                               maxiter=100, constraint=consfunc)
print(popt)