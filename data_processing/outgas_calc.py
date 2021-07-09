import numpy as np

def get_initial_concentration(solubility, abundance, pressure, initial_temp):
	"""Returns the initial concentration of gas in the sample, in atoms per cubic centimeter.

	Keyword arguments:

	solubility -- fractional volume of polymer that the gas occupies
	pressure -- internal gas pressure within sample before outgassing in pascals
	initial_temp -- temperature of sample (in kelvin) before heating
	"""
	k_b = 1.381e-23 # J/k
	# recall that 1 Pa / 1 J = 1 m^-3
	return solubility * abundance * pressure / (k_b * initial_temp) / 1000000 # factor of 1000000 to convert to cm^-3

def arrhenius_relation(D_0, E_a, T):
	"""Returns the diffusion for each temeperature in T. Makes use of the Boltzmann constant in eV per
	kelvin [approximated to 9 decimal places, NOT EXACT].

	Keyword arguments:
	
	D_0 -- The diffusion constant at infinite temperature in m^2/s
	E_a -- the activation energy of the gas/polymer combination in eV
	T -- array-like: temperatures of the sample in kelvin at each time value
	"""
	k_b = 8.617333262e-5 # eV/K
	return D_0 * np.exp((-1.0)*E_a/(k_b*T))

def get_concentrations(c_0, diff_const, thickness, times):
	"""Returns the concentration of gas in the sample at each corresponding time value.

	Keyword arguments:
	c_0 -- initial concentration
	diff_const -- numpy array of calculated diffusion constants of gas/polymer pair
	thickness -- thickness of sample
	times -- numpy array of time values to use as inputs
	temperatures -- numpy array of temperature values corresponding to times
	"""
	concentrations = np.zeros(len(times))

	for idx, t in enumerate(times):

		coef = 8.0 * c_0 * thickness / (np.pi**2)
		
		# calculate sum
		tot = 0
		for n in range(1000):
			tot += np.exp(-1.*(((np.pi*(2*n+1))/thickness)**2)*diff_const[idx]*t)/((2*n+1)**2)

		tot = coef * tot
		concentrations[idx] = tot

	return concentrations

def get_flow_rate(c_0, diff_const, thickness, times):
	"""Returns the outgassing flow rate in the sample at each corresponding time value.

	Keyword arguments:
	c_0 -- initial concentration
	diff_const -- diffusion constant of gas
	thickness -- thickness of sample
	times -- numpy array of time values to use as inputs
	"""
	rates = np.zeros(len(times))

	for idx, t in enumerate(times):

		coef = 4 * c_0 * diff_const[idx] / thickness

		# calculate sum
		tot = 0
		for n in range(1000):
			tot += np.exp(-((np.pi*(2*n+1))/thickness)**2 * diff_const[idx] * t)

		tot = coef * tot
		rates[idx] = tot

	return rates
