from scipy.interpolate import UnivariateSpline

# In example below, using 'Orion' as cluster. Can change the name for any cluster.

percentage_to_probability_fit = UnivariateSpline((np.cumsum(p_binary_fraction['Orion']) / np.sum(p_binary_fraction['Orion']))[p_binary_fraction['Orion']>1E-10], binary_fraction['Orion'][p_binary_fraction['Orion'] > 1E-10], s=0)

# Notes : We are making the Inverse cumulative probability distribution. So, if you feed it 50%, it will pop out the right probability.
# the > 1E-10 is just a trick to ignore all values that have probability less than 1E-10
# s = 0 is the smoothing factor. Setting it to zero ensures that the interpolater goes through all the points.

# To get the 16th, 50th, and 84th percentiles:

#50th percentile
percentage_to_probability_fit(0.5)
# and so on for the 84th and 16th; repeat for all clusters.


