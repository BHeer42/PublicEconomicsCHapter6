# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:43:17 2020

@author: heerburk

Ch6_social_security1.py

Chapter 6.3., Public Economics (2019) by Burkhard Heer

computes the transition and welfare effects of a
PAYG pension in a simple 2-period OLG model

"""
# part 1: import libraries
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['figure.figsize'] = (10,6)
log = math.log
e = math.e

# part 2: parameterization
labor = 0.30    # constant labor supply
alpha = 0.36        # production elasticity of capital
b = 0.40        # discount factor
n = 0.10        # population growth rate
nu0 = 257.15    # disutility of labor
nu1 = 3.33      # inverse of Frisch elasticity
tau = (0, 0.30) # parameter values for pension replacement rate (tuple)
nt=10           # number of transition periods   
             
# part 3: defining of the functions
#
# wage function
def wage(k,l):
    return (1-alpha) * k**alpha * l**(-alpha)

# interest rate function
def interest_rate(k,l):
    return alpha * k**(alpha - 1) * l**(1-alpha)

# steady state computation
# input: guess for k
# output: steady-state condition (=0 in steady state)
def ksteady(k):
    w = wage(k,labor)           # wage
    r = interest_rate(k,labor)  # interest rate
    d = tau0 * w * labor        # pension
    # non-linear equation (6.19) in Heer (2019)
    y = (1+n)*k-w*labor*b/(1+b) + 1/(1+b)*(1+b+b*r+n)/(1+r)*d
    # reformulatoin of (6.19) so that k=0 is not a solution
    # we are looking for non-trivial steady state k>0
    y= y/k
    return y

# function that computes k_t+1 given k_t
def kdyn(k1):
    w0 = wage(k0,labor)          # wage in t
    w1 = wage(k1,labor)          # wage in t+1
    r1 = interest_rate(k1,labor) # interest rate in t+1
    d1 = tau1 * w1 * labor        # pension in t+1
    y = (1+n)*k1 - (1-tau0)*w0*labor*b/(1+b) + 1/(1+b)*(1+n)/(1+r1) * d1
    return y
    
# part 4: main program
tau0 = tau[0]       # initial steady state: no pension  


# 4.1 
# come up with an initial guess for steady state value
# by searching over a grid of values for k
kgrid = np.linspace(0.001, 1, 100)
ygrid = np.zeros(100)   # values of non-linear eq for k
for i in range(100):
    ygrid[i] = ksteady(kgrid[i])

plt.plot(kgrid, ygrid, 'b-', linewidth=2)
plt.show()

ygrid = abs(ygrid)
print("Minimum of ygrid: " + str(ygrid.min()))

i0 = np.where(ygrid == ygrid.min()) # find index of y with minimum
kguess = kgrid[i0]
print("Index of Minimum: " + str(kguess))
print("ksteady for initial guess: " + str(ksteady(kguess)))

# 4.2
# compute initial steady state with tau=0 and steady-state lifetime utility
kss = scipy.optimize.fsolve(ksteady, kguess)
print(str(kss))
wss = wage(kss,labor)
rss = interest_rate(kss, labor)
dss = tau0 * wss * labor
c1ss = 1/(1+b) * (wss*labor + (n-rss)/(1+rss)*dss)
c2ss = b*c1ss*(1+rss)
utilss = log(c1ss) + b*log(c2ss)-nu0*labor**(1+nu1)/(1+nu1)
print(utilss)


# 4.3
# compute final steady state with tau=0.3 and steady-state lifetime utility
tau0 = tau[1]
# find new initial value
for i in range(100):
    ygrid[i] = ksteady(kgrid[i])

ygrid = abs(ygrid)

i0 = np.where(ygrid == ygrid.min()) # find index of y with minimum
kguess = kgrid[i0]
kssd = scipy.optimize.fsolve(ksteady, kguess)
print(ksteady(kssd))
print(str(kssd))
wssd = wage(kssd,labor)
rssd = interest_rate(kssd, labor)
dssd = tau0 * wssd * labor
c1ssd = 1/(1+b) * (wssd*labor + (n-rssd)/(1+rssd)*dssd)
c2ssd = b*c1ssd*(1+rssd)
utilssd = log(c1ssd) + b*log(c2ssd)-nu0*labor**(1+nu1)/(1+nu1)
print(utilssd)

# 4.4
# compute steady state welfare effect
cec = e**((utilssd-utilss)/(1+b) )-1

print("Welfare effect in steady state from the introduction of a pension:")
print("consumption equivalent change: " + str(cec))
#print("cec: " + str(cec))

# 4.5
# computation of the dynamics
# from initial to final steady state

kt = np.zeros(nt+1)     # time path for capital stok in periot t
utilt = np.zeros(nt)  # lifetime utility of generation born in t
periods = range(nt)   # range of period 0,..., nt
cect = np.zeros(nt)   # welfare effects, % of consumption

kt[0] = kssd            # initial value of capital stock in period 0
utilt[0] = utilssd
cect[0] = e**((utilt[0]-utilss)/(1+b) )-1

tau0 = 0.3    # tax rate at young age
tau1 = 0.3    # tax rate at old age


# i=0: Period 0 with steady state, tau=30% for young and old
# i=1: Period 1, young generation still has to finance old agents 
# (tau0=30% so that d_1=tau_0 w_1 l), but will not receive
# a pension in old age in period 2 (i=2) with tau1=0% so that d_2=0

for i in periods:
    k0=kt[i]
    print(i,k0)
    if i==2: 
        tau1 = 0  # first generation still has to pay contributions 
    if i==3: 
        tau0 = 0  # first generation which does not have to pay contributions
	
    k1 = scipy.optimize.fsolve(kdyn, k0)
    kt[i+1] = k1
    
    w0 = wage(k0,labor)
    w1 = wage(k1,labor)
    r1 = interest_rate(k1, labor)
    d0 = tau0 * w0 * labor    
    d1 = tau1 * w1 * labor

    c1 = 1/(1+b)*(w0*labor-d0+d1*(1+n)/(1+r1))
    c2 = b * c1 * (1+r1)
    utilt[i] = log(c1) + b*log(c2) - nu0*labor**(1+nu1) / (1+nu1)
    cect[i] = e**((utilt[i]-utilss) / (1+b) ) -1
	
fig, axes = plt.subplots(2, 1, figsize=(8, 16))
axes[0].set_xlabel('time')
axes[0].set_ylabel('capital')
axes[0].plot(kt[1:nt])
axes[1].set_xlabel('time')
axes[1].set_ylabel('welfare')
axes[1].plot(cect[1:nt])
plt.show()                