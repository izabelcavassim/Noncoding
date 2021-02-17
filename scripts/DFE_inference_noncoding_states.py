
"""
Uses dadi and fitdadi to infer the DFE of a given synonynmous sfs.

Script written by Jon 201907011, and adapted by Maria Izabel Cavassim Alves
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
import matplotlib.pyplot as pyplot
import pylab
import pandas as pd


# transcribed_and_enhancer_states
# transcription_states
# promoters_states


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
	"""Like *argparse.ArgumentParser*, but prints help when no arguments."""

	def error(self, message):
		"""Print error message, then help."""
		sys.stderr.write('error: %s\n\n' % message)
		self.print_help()
		sys.exit(2)


class DemographicAndDFEInference():
	"""Wrapper class to allow functions to reference each other."""

	def ExistingFile(self, fname):
		"""Return *fname* if existing file, otherwise raise ValueError."""
		if os.path.isfile('/Users/PM/Dropbox/Noncoding/results/VCF_expression/'+fname):
			return fname
		else:
			raise ValueError("%s must specify a valid file name" % fname)

	def Existingtreatment(self, treatment):
		"""Return *treatment* if existing, else raise ValueError.

		# 1) Dexamethasone in PBMC (that's some sort of blood cell) -> "PBMC_Dex"
		# 2) Vitamin D in melanocytes -> "Mel_VitD"
		# 3) Aldosterone in muscle -> "SMC_Ald"
		# 4) nicotine in PBMC -> "PBMC_Nicotine"
		# 5) Vitamin A in melanocytes or UBC -> "Mel_RetAc" or
		""" 

		if treatment == "transcribed_and_enhancer_states":
			return treatment
		elif treatment == "transcription_states":
			return treatment
		elif treatment == "promoters_states":
			return treatment
		else:
			raise ValueError('%s must specify a valid treatment.' % treatment)

	def inferDFEParser(self):
		"""Return *argparse.ArgumentParser* for ``fitdadi_infer_DFE.py``."""
		parser = ArgumentParserNoArgHelp(
			description=(
				'Given the number of individuals in population one and two, '
				'this script outputs a `*pops_file.txt` in the format '
				'specified for use by the python package, `easySFS.py`.'),
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)
		parser.add_argument(
			'syn_input_sfs', type=self.ExistingFile,
			help=('full path to the vcf file to compute the "Synonynomous" site-frequency spectrum from which the '
				  'demographic parameters should be inferred.'))
		parser.add_argument(
			'nonsyn_input_sfs', type=self.ExistingFile,
			help=('full path to the vcf file to compute the Nonsynonynomous site-frequency spectrum from which the '
				  'distribution of fitness effects should be inferred.'))
		parser.add_argument(
			'outprefix', type=str,
			help='The file prefix for the output `*inferred_demography.txt`.')
		parser.add_argument(
			'--combined', type=str,
			help=('Boolean parameter to compute the DFE by combining deg and nodeg treatments'), 
			default='False')
		parser.add_argument(
			'--treatment', type=self.Existingtreatment,
			help=('The treatment of organism from which the input `.vcf` is drawn from. '
				'Must be "transcribed_and_enhancer_states" \r\n'
				'Must be "transcription_states" \r\n'
				'Must be "promoters_states" \r\n'),
			default="transcription_states")
		return parser

	def snm(self, notused, ns, pts):
		"""Return a standard neutral model.
		ns = (n1, )
			n1: Number of samples in resulting Spectrum object.
		pts: Number of grid points to use in integration.
		"""
		xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
		phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

		# Construct Spectrum object.
		fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
		return fs

	def two_epoch(self, params, ns, pts):
		"""Define a two-epoch demography, i.e., an instantaneous size change.

		params = (nu, T)
			nu: Ratio of contemporary to ancient population size.
			T: Time in the past at which size change occured,
				in units of 2*N_a.
		ns = (n1, )
			n1: Number of samples in resulting Spectrum object.
		pts: Number of grid points to use in integration.
		"""
		nu, T = params  # Define given parameters.
		xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
		phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

		phi = dadi.Integration.one_pop(phi, xx, T, nu)  # Integrate.

		# Construct Spectrum object.
		fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
		return fs

	def two_epoch_sel(self, params, ns, pts):
		"""Define a two-epoch demography, i.e., an instantaneous size change.

		This method incorporates a gamma parameter.

		params = (nu, T, gamma)
			nu: Ratio of contemporary to ancient population size.
			T: Time in the past at which size change occured,
				in units of 2*N_a.
			gamma: Parameter tuple describing a gamma distribution.
		ns = (n1, )
			n1: Number of samples in resulting Spectrum object.
		pts: Number of grid points to use in integration.
		"""
		nu, T, gamma = params  # Define given parameters.
		xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
		phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)  # Define initial phi

		phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)

		# Construct Spectrum object.
		fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
		return fs

	def gamma_dist(self, mgamma, alpha, beta):
		"""Define a gamma distribution.

		self: reference to this instance of a gamma distribution.
		mgamma: float which describes the mean value of this gamma
			 distribution.
		alpha: shape parameter of the gamma distribution.
		beta: scale parameter of the gamma distribution.
		"""
		return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta)

	def neugamma(self, mgamma, pneu, alpha, beta):
		"""Define a neutral-gamma distribution.

		self: reference to this instance of a neutral-gamma distribution.
		mgamma: float which describes the mean value of this neutral-gamma
			distribution.
		pneu: proportion of elements which are assumed to be neutral, i.e.,
			equal to 0.
		alpha: shape parameter of the non-neutral elements of the
			neutral-gamma distribution.
		beta: scale parameter of the non-neutral elements of the
			neutral-gamma distribution.
		"""
		mgamma = -mgamma
		# Assume anything with gamma < 1e-4 is neutral
		if (0 <= mgamma) and (mgamma < 1e-4):
			return pneu / (1e-4) + (1 - pneu) * self.gamma_dist(
				-mgamma, alpha, beta)
		else:
			return self.gamma_dist(-mgamma, alpha, beta) * (1 - pneu)

	def compute_sfs(self, syn_deg, non_deg):
		"""Return both syn SFS and non-syn SFS by vcf files provided in dadi format.

		"""
		pop_ids, ns = ['YRI'],[108] #, 'CEU'], [20,24]

		# Parse the VCF file to generate a data dictionary
		dd_syn = dadi.Misc.make_data_dict_vcf('/Users/PM/Dropbox/Noncoding/data/iza_expression/'+syn_deg, '/Users/PM/Dropbox/Noncoding/software/dadi/examples/fs_from_data/1KG.YRI.CEU.popfile.txt')
		dd_nonsyn = dadi.Misc.make_data_dict_vcf('/Users/PM/Dropbox/Noncoding/data/iza_expression/'+non_deg, '/Users/PM/Dropbox/Noncoding/software/dadi/examples/fs_from_data/1KG.YRI.CEU.popfile.txt')

		# If we had outgroup information, we could unfold the fs.
		fs_folded_syn = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns, polarized=False)
		fs_folded_nonsyn = dadi.Spectrum.from_data_dict(dd_nonsyn, pop_ids, ns, polarized=False)

		return([fs_folded_syn, fs_folded_nonsyn])

	def main(self):
		"""Execute main function."""
		# Parse command line arguments
		parser = self.inferDFEParser()
		args = vars(parser.parse_args())
		prog = parser.prog

		# Assign arguments
		syn_input_sfs = args['syn_input_sfs']
		nonsyn_input_sfs = args['nonsyn_input_sfs']
		outprefix = args['outprefix']
		treatment = args['treatment']
		combined = args['combined']

		# Numpy options
		numpy.set_printoptions(linewidth=numpy.inf)

		# create output directory if needed
		outdir = os.path.dirname(args['outprefix'])
		if outdir:
			if not os.path.isdir(outdir):
				if os.path.isfile(outdir):
					os.remove(outdir)
				os.mkdir(outdir)

		# Output files: logfile
		# Remove output files if they already exist
		underscore = '' if args['outprefix'][-1] == '/' else '_'
		inferred_demography = \
			'{0}{1}inferred_demography.txt'.format(
				args['outprefix'], underscore)
		inferred_DFE = \
			'{0}{1}inferred_DFE.txt'.format(args['outprefix'], underscore)
		logfile = '{0}{1}log.log'.format(args['outprefix'], underscore)
		to_remove = [logfile, inferred_demography, inferred_DFE]
		for f in to_remove:
			if os.path.isfile(f):
				os.remove(f)

		# Set up to log everything to logfile.
		logging.shutdown()
		logging.captureWarnings(True)
		logging.basicConfig(
			format='%(asctime)s - %(levelname)s - %(message)s',
			level=logging.INFO)
		logger = logging.getLogger(prog)
		warning_logger = logging.getLogger("py.warnings")
		logfile_handler = logging.FileHandler(logfile)
		logger.addHandler(logfile_handler)
		warning_logger.addHandler(logfile_handler)
		formatter = logging.Formatter(
			'%(asctime)s - %(levelname)s - %(message)s')
		logfile_handler.setFormatter(formatter)
		logger.setLevel(logging.INFO)

		# print some basic information
		logger.info('Beginning execution of {0} in directory {1}\n'.format(
			prog, os.getcwd()))
		logger.info('Progress is being logged to {0}\n'.format(logfile))
		logger.info('Parsed the following arguments:\n{0}\n'.format(
			'\n'.join(['\t{0} = {1}'.format(*tup) for tup in args.items()])))


		# Construct initial Spectrum object from input synonymous sfs.
		print(syn_input_sfs)
		print(nonsyn_input_sfs)

		if combined == 'False':
			syn_data, nonsyn_data  = self.compute_sfs(syn_input_sfs, nonsyn_input_sfs)
			#syn_data = dadi.Spectrum.from_file(syn_input_sfs)

		# Adding the SFSs of deg and nodeg treatments for hypothesis testing 
		if combined == 'True':
			base_name = syn_input_sfs.split("deg")[0]
			print(base_name)
			syn_data, nonsyn_data = self.compute_sfs_both_treatments(base_name)	
			treatment = treatment+"combined"

		print(syn_data)
		syn_ns = syn_data.sample_sizes  # Number of samples.
		pts_l = [800, 1000, 1200] # Number of grid points to use in the integration

		# Optomize parameters for this model.
		# First set parameter bounds for optimization
		upper_bound = [8, 3]
		lower_bound = [1e-3, 0]
		initial_guess = [0.5, 0.1]

		############################################################################
		#         Infering the demographic model based on the syn sites            #
		############################################################################
 
		with open(inferred_demography, 'w') as f:
			f.write('Beginning with demographic inference.\n')
			max_likelihood = -1e25
			for i in range(25):

				# Start at initial guess
				p0 = initial_guess

				# Randomly perturb parameters before optimization.
				p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
				
				# Make the extrapolating version of demographic model function.
				func_ex = dadi.Numerics.make_extrap_log_func(self.two_epoch)
				
				logger.info('Beginning optimization with guess, {0}.'.format(p0))
				
				popt = dadi.Inference.optimize_log_lbfgsb(
					p0=p0, data=syn_data, model_func=func_ex, pts=pts_l,
					lower_bound=lower_bound, upper_bound=upper_bound,
					verbose=len(p0), maxiter=30)
				
				logger.info('Finished optimization with guess, ' + str(p0) + '.')

				logger.info('Best fit parameters: {0}.'.format(popt))
				
				# Calculate the best-fit model allele-frequency spectrum (AFS).
				# Note, this spectrum needs to be multiplied by "theta".
				non_scaled_spectrum = func_ex(popt, syn_ns, pts_l)
				
				# Likelihood of the data given the model AFS.
				multinomial_ll_non_scaled_spectrum = dadi.Inference.ll_multinom(model=non_scaled_spectrum, data=syn_data)

				logger.info('Maximum log composite likelihood: {0}.'.format(multinomial_ll_non_scaled_spectrum))

				# Calculating theta: 4N\muL (synonymous).
				theta = dadi.Inference.optimal_sfs_scaling(non_scaled_spectrum, syn_data)

				logger.info('Optimal value of theta: {0}.'.format(theta))

				if multinomial_ll_non_scaled_spectrum > max_likelihood:
					best_params = popt
					best_non_scaled_spectrum = non_scaled_spectrum
					max_likelihood = multinomial_ll_non_scaled_spectrum
					theta_syn = theta


			# Plotting data and the fit SFS
			fig = pylab.figure(1)
			fig.clear()
			dadi.Plotting.plot_1d_comp_multinom(best_non_scaled_spectrum, syn_data)
			fig.savefig('1d_comp_{treatment}.pdf'.format(treatment=treatment))

			# Scaling the spectrum by theta
			best_scaled_spectrum = theta_syn * best_non_scaled_spectrum
			
			# Here we multiply it by 2.31 because of transition transversion bias and CpG mutation rate (Eyre-Walker 2006)
			theta_nonsyn = theta_syn * 2.31 
			
			poisson_ll = dadi.Inference.ll(model=best_scaled_spectrum, data=syn_data)

			f.write('Best fit parameters: {0}.\n'.format(best_params))
			f.write('Maximum multinomial log composite likelihood: {0}.\n'.format(max_likelihood))
			f.write('Maximum poisson log composite likelihood: {0}.\n'.format(poisson_ll))
			f.write('Non-scaled best-fit model spectrum: {0}.\n'.format(best_non_scaled_spectrum))
			f.write('Optimal value of theta_syn: {0}.\n'.format(theta_syn))
			f.write('Optimal value of theta_nonsyn: {0}.\n'.format(theta_nonsyn))
			f.write('Scaled best-fit model spectrum: {0}.\n'.format(best_scaled_spectrum))
			logger.info('Finished demographic inference.')

			############################################################################
			#                   Fitting the DFE to the SFS 							   #
			############################################################################
			logger.info('Finished demographic inference.')
			logger.info('Beginning DFE inference.')
			nonsyn_ns = nonsyn_data.sample_sizes

			demog_params = best_params

			# This is the exon length for all genes within each treatment
			Mel_RetAcdeg=9708623
			Mel_RetAcnodeg=47292474
			#Together:
			
			Mel_RetAc=Mel_RetAcdeg+Mel_RetAcnodeg

			Mel_VitDdeg=17481424
			Mel_VitDnodeg=45531498

			#Together: 
			Mel_VitD=Mel_VitDdeg+Mel_VitDnodeg

			PBMC_Dexdeg=34208438
			PBMC_Dexnodeg=44831273
			PBMC_Dex=PBMC_Dexdeg+PBMC_Dexnodeg

			PBMC_Nicotinedeg=13067737
			PBMC_Nicotinenodeg=51844665

			PBMC_Nicotine=PBMC_Nicotinedeg+PBMC_Nicotinenodeg

			SMC_Alddeg=11928678
			SMC_Aldnodeg=57159923

			SMC_Ald=SMC_Alddeg+SMC_Aldnodeg

			All_genes=96196469

			Total_length = {"PBMC_Dexdeg":PBMC_Dexdeg,"PBMC_Dexnodeg":PBMC_Dexnodeg, "Mel_VitDdeg":Mel_VitDdeg,"Mel_VitDnodeg":Mel_VitDnodeg, "SMC_Alddeg":SMC_Alddeg, 
			"SMC_Aldnodeg":SMC_Aldnodeg,"PBMC_Nicotinedeg": PBMC_Nicotinedeg, "PBMC_Nicotinenodeg":PBMC_Nicotinenodeg,"Mel_RetAcdeg":Mel_RetAcdeg, "Mel_RetAcnodeg":Mel_RetAcnodeg, "All_genes":All_genes}

			Total_length_combined = {"PBMC_Dex":PBMC_Dex,"Mel_VitD":Mel_VitD,"SMC_Ald":SMC_Ald,"PBMC_Nicotine":PBMC_Nicotine,"Mel_RetAc":Mel_RetAc}

			# Estimating the Length of synonymous sites based on 2.31 ratio		
			Lsyn_dict = {k: v /(2.31+1) for k, v in Total_length.items()}
			Lsyn_dict_combined = {k: v /(2.31+1) for k, v in Total_length_combined.items()}

			if combined == "False":
				# Theta parameters = Thete = 4*Na*u*Lsyn
				Lsyn = Lsyn_dict[treatment]  # Length of synonymous sites.
			else:
				Lsyn = Lsyn_dict_combined[base_name]

				print("This is basename")
				print(base_name)
				print(Lsyn)

			u = 1.5e-08
			###u = 2.5x10-8, 1.5e-08 # Huber et al 2016
			#u= 1.5x10-8 or 1.8x10-8 (because of compatibiity with Bokyo et al 2020)

			### Effective population size
			Na = theta_syn / (4 * u * Lsyn)

			# Reasoning behind this???
			#max_s = 0.5
			#max_gam = max_s * 2 * Na
			int_bounds=(1e-5 ,500)

			ptsl = [600, 800, 1000]

			# Generating the spectra for a range of gammas:
			spectra = Selection.spectra(demog_params, nonsyn_ns,
									self.two_epoch_sel,
									pts_l=pts_l, int_bounds=int_bounds,
									Npts=300, echo=True, mp=True)

			#BETAinit = max_gam / 3
			#initial_guess = [0.09, BETAinit]
			#upper_beta = 10 * max_gam
			#lower_bound = [1e-3, 0]
			#upper_bound = [1, upper_beta]

			sel_params = [0.2, 1000.]
			lower_bound = [1e-3, 1e-2]
			upper_bound = [1, 50000.]


			gamma_max_likelihoods = []
			gamma_guesses = dict()

			############################################################################ 
			#                  Fitting a gamma model to the DFEs			           #
			############################################################################
			for i in range(25):
				p0 = sel_params
				p0 = dadi.Misc.perturb_params(p0, lower_bound=lower_bound,
										  upper_bound=upper_bound)
				logger.info('Beginning optimization with guess, {0}.'.format(p0))
				popt = numpy.copy(Selection.optimize_log(p0, nonsyn_data,
													 spectra.integrate,
													 self.gamma_dist,
													 theta_nonsyn,
													 lower_bound=lower_bound,
													 upper_bound=upper_bound,
													 verbose=len(p0),
													 maxiter=30))
				logger.info('Finished optomization, results are {0}.'.format(popt))

				gamma_max_likelihoods.append(popt[0])
				gamma_guesses[popt[0]] = popt

			# Sorted data
			gamma_max_likelihoods.sort()

			############################################################################
			#        Fitting more complex models to the DFEs : neutral+gamma           #
			############################################################################
			neugamma_vec = numpy.frompyfunc(self.neugamma, 4, 1)
			sel_params = [0.2 , 0.2 , 1000.]
			#initial_guess = [0.999999999, 0.09, BETAinit]
			lower_bound = [1e-3, 1e-3, 1e-2]
			upper_bound = [1, 1, 50000.] # proportion of neutral alleles 
			neugamma_max_likelihoods = []
			neugamma_guesses = dict()
			for i in range(25):
				p0_neugamma = sel_params # instead of initial guesses
				p0_neugamma = dadi.Misc.perturb_params(p0_neugamma,
												   lower_bound=lower_bound,
												   upper_bound=upper_bound)
				logger.info('Beginning optimization with guess, {0}.'.format(
				p0_neugamma))
				popt = numpy.copy(Selection.optimize_log(p0_neugamma, nonsyn_data,
													 spectra.integrate,
													 neugamma_vec,
													 theta_nonsyn,
													 lower_bound=lower_bound,
													 upper_bound=upper_bound,
													 verbose=len(p0_neugamma),
													 maxiter=30))
				logger.info('Finished optimization, results are {0}.'.format(popt))
				
				neugamma_max_likelihoods.append(popt[0])
				neugamma_guesses[popt[0]] = popt

			# Sorted data
			neugamma_max_likelihoods.sort()

			logger.info('Integrating expected site-frequency spectrum.')

			logger.info('Outputing results.')

			results_gamma = list()
			results_neugamma = list()
			with open(inferred_DFE, 'w') as f:
				f.write('Assuming a gamma-distributed DFE...\n')
				f.write('Outputting best 5 MLE estimates.\n')

				# Outputting gamma guesses 
				print(gamma_guesses)
				for i in range(len(gamma_guesses.keys())-1):
					best_popt = gamma_guesses[gamma_max_likelihoods[i]]
					results_gamma.append([best_popt[0], best_popt[1][0], best_popt[1][1]])
					expected_sfs = spectra.integrate(best_popt[1], self.gamma_dist, theta_nonsyn)
					f.write('The population-scaled best-fit parameters: {0}.\n'.format(
							best_popt))
					# Divide output scale parameter by 2 * N_a
					f.write('The non-scaled best-fit parameters: ''[{0}, array({1})].\n'.format(
						best_popt[0],
						numpy.divide(best_popt[1], numpy.array([1, 2 * Na]))))
				f.write('The expected SFS is: {0}.\n\n'.format(expected_sfs))
				f.write('Assuming a neutral-gamma-distributed DFE...\n')
				f.write('Outputting best 25 MLE estimates.\n')
				
				# Outputting neutral + gamma guesses 
				for i in range(len(neugamma_guesses.keys())-1):
					best_popt_neugamma = neugamma_guesses[neugamma_max_likelihoods[i]]
					expected_sfs_neugamma = spectra.integrate(best_popt_neugamma[1], neugamma_vec, theta_nonsyn)
					results_neugamma.append([best_popt_neugamma[0], best_popt_neugamma[1][0], best_popt_neugamma[1][1], best_popt_neugamma[1][2]])
					f.write('The population-scaled best-fit parameters: {0}.\n'.format(best_popt_neugamma))
					
					# Divide output scale parameter by 2 * N_a
					f.write(
						'The non-scaled best-fit parameters: ''[{0}, array({1})].\n'.format(
						best_popt_neugamma[0],numpy.divide(best_popt_neugamma[1],numpy.array([1, 1, 2 * Na]))))
				
					f.write('The expected SFS is: {0}.\n\n'.format(expected_sfs_neugamma))
			

			# Writing results as csvs:
			df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood', 'alpha', 'beta'], dtype = float) 
			df_gamma.to_csv('DFE_inference_gamma_{treatment}.csv'.format(treatment=treatment), index=False)

			df_neu_gamma = pd.DataFrame(results_neugamma, columns =['Likelihood', 'pneu', 'alpha', 'beta'], dtype = float)
			df_neu_gamma.to_csv('DFE_inference_neugamma_{treatment}.csv'.format(treatment=treatment), index=False)

			logger.info('Pipeline executed succesfully.')

if __name__ == '__main__':
	DemographicAndDFEInference().main()

## These are the optimal parameters when the spectrum is folded. They can be
## found simply by passing fold=True to the above call to optimize_log.
#pfold =  array([1.907,  0.073,  1.830,  0.899,  0.425,  0.113])
#
## The interface to hessian computation is designed for general functions, so we
## need to define the specific functions of interest here. These functions
## calculate -ll given the logs of the parameters. (Because we work in log
## parameters, the uncertainties we estimate will be *relative* parameter
## uncertainties.)
#from dadi.Inference import ll_multinom
#func = lambda lp: -ll_multinom(func_ex(numpy.exp(lp), ns, pts_l), data)
#foldfunc = lambda lp: -ll_multinom(func_ex(numpy.exp(lp), ns, pts_l).fold(), 
#                                   data.fold()) 
#
## Calculate the two hessians
#h = dadi.Hessian.hessian(func, numpy.log(params), 0.05)
#hfold = dadi.Hessian.hessian(foldfunc, numpy.log(pfold), 0.05)
#
## Now we calculate the *relative* parameter uncertainties.
#uncerts = numpy.sqrt(numpy.diag(numpy.linalg.inv(h)))
#uncerts_folded = numpy.sqrt(numpy.diag(numpy.linalg.inv(hf)))
#
## The increase in uncertainty is not too bad. Tp increasing by 50% is the only
## substantial one.
#print uncerts_folded/uncerts - 1