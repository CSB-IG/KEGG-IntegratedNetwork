import powerlaw
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--network", help="Network file", dest='network',required=True, metavar="AllKEGG.ncbi.filter.txt")

args = parser.parse_args()

network_file = args.network

####### NETWORK DEFINED BY ALL CONNECTED NODES
##########################################
##########################################

data = pd.read_csv(network_file)
indegree = data['Indegree'].values
outdegree = data['Outdegree'].values
degree = indegree + outdegree


fit = powerlaw.Fit(degree)

print("ALPHA")
print(fit.power_law.alpha)

print("SIGMA")
print(fit.power_law.sigma)

print("EXPONENTIAL")
R, p = fit.distribution_compare('power_law', 'exponential')
print (R, p)

print("TRUNCATED POWER LAW")
R, p = fit.distribution_compare('power_law', 'truncated_power_law')
print (R, p)

print("LOGNORMAL")
R, p = fit.distribution_compare('power_law', 'lognormal')
print (R, p)

print("LOGNORMAL POSITION")
R, p = fit.distribution_compare('power_law', 'lognormal_positive')
print (R, p)

