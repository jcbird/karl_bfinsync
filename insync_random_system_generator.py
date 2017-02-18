import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib as mlt
import astropy
from scipy.stats import lognorm
import h5py
import random

grav = 6.67e-8
mass_sol = 1.99e33
rad_sol = 6.955e10

min_abs_period = 0.3 * 86400.
max_abs_period = 1.e6 * 86400.

min_abs_distance = 0.014 * 1.496e13
max_abs_distance = 391.42 * 1.496e13

min_abs_mass = 0.01
max_abs_mass = 1.0

def distance_calc(period,m1,m2):

    g = grav
    mass1 = m1 * mass_sol
    mass2 = m2 * mass_sol

    p1 = period
    G = g**0.5

    M2 = mass2**(3./2.)
    PI = np.pi**(-1.)
    M1 = mass1**(-1.)
    SUM = (1. + (mass2 / mass1)) ** (-1.)

    r1 = (p1 * G * M2 * PI * M1 * SUM) ** (2./3.)

    distance = r1 * (1. + (mass1/mass2))
    
    
    return distance



def velocity_calc(distance,m1,m2):
    dist = distance
    mass1 = m1*mass_sol
    mass2 = m2*mass_sol

    prim_dist = dist / ( 1. + (mass1/mass2) )
    sec_dist = dist - prim_dist
    
    a = grav
    b = mass2**2.
    c = ( sec_dist * mass1 )
    d = ( 1. + ( mass2/mass1 )) ** 2.

    velocity2 = np.sqrt(a*b*(c**-1))

    return velocity2


def roche_lobe_check(period,m1,m2,**options):
    sys_validity = False
    if options.get("action") == 'check':
        separation = distance_calc(period,m1,m2)

    mass1 = m1*mass_sol
    mass2 = m2*mass_sol
    q = mass2 / mass1 
    
    ratio = ( 0.49 * q**(-2./3.) ) / ( .6 * q**(-2./3.) + np.log(1 + q**(-1./3.)))

    roche_lobe = 2. * rad_sol
    LHS = roche_lobe / separation
    RHS = ratio
    sys_validity = True
    if LHS >= RHS:
##        print "Your system is not a physically distinct binary. \nMass accretion effects present."
        sys_validity = False
    return sys_validity



def binary_system_generator(N):
    m1 = 1.0
    period_arr = []
    distance_arr = []
    companion_arr = []
    velocity_arr = []
    phase_arr = []
    count=0
    for i in range(np.int(5*N)):
        periodx = (10.**(np.random.normal(loc=4.8,scale=2.3,size=1))) * 84600.
        companionx =  np.random.power(a=2.8, size=1)

        condition1 = np.logical_and(periodx >= min_abs_period,
                                    periodx <= max_abs_period)
        
        condition2 = companionx >= 0.02

        if condition1 == False or condition2 == False:
            continue

        if condition1 == True and condition2 == True:
            sys_separation = distance_calc(periodx,m1,companionx)
            
        sys_validity = roche_lobe_check(periodx,m1,companionx,action='check')

        if sys_validity == True:
            sys_velocity = velocity_calc(sys_separation,m1,companionx)
            
        if sys_validity == False :
            continue
            
        if condition1 == True and condition2 == True:
            count += 1
            period_arr.append(periodx)
            distance_arr.append(sys_separation)
            velocity_arr.append(sys_velocity)
            companion_arr.append(companionx)
            phase_arr.append(random.uniform(0.,2.*np.pi))
            
        print np.float(i)*100. / (5.*N), '% simulations possible'            
        print np.float(count)*100./np.float(N),'% Useful systems produced \033[2A'

        if len(period_arr) ==  N:
            break
    return period_arr,distance_arr,velocity_arr,companion_arr,phase_arr
      

def binary_system_generator_v2(N):
    period_arr = []
    distance_arr = []
    primary_arr = []
    companion_arr = []
    velocity_arr = []
    phase_arr = []
    count=0
    for i in range(np.int(5*N)):
        periodx = lognorm.rvs(2.3,scale = np.exp(4.8), size = 1) * 86400.
        companionx = 0.158*0.578*lognorm.rvs(0.69,scale=np.exp(0.08),size = 1) #primary mass
        companiony = 0.158*0.578*lognorm.rvs(0.69,scale=np.exp(0.08),size = 1) #companion mass

        condition1 = np.logical_and(periodx >= min_abs_period,
                                    periodx <= max_abs_period)
        
        condition2 = np.logical_and(companionx >= min_abs_mass,
                                    companionx <= max_abs_mass)

        condition3 = companiony/companionx <= 1.0

        if condition1 == False or condition2 == False or condition3 == False:
            continue

        if condition1 == True and condition2 == True and condition3 == True:
            sys_separation = distance_calc(periodx,companionx,companiony)
            
        sys_validity = roche_lobe_check(periodx,companionx,companiony,action='check')

        if sys_validity == True:
            sys_velocity = velocity_calc(sys_separation,companionx,companiony)
            
        if sys_validity == False :
            continue
            
        if condition1 == True and condition2 == True:
            count += 1
            period_arr.append(periodx)
            distance_arr.append(sys_separation)
            velocity_arr.append(sys_velocity)
            primary_arr.append(companionx)
            companion_arr.append(companiony)
            phase_arr.append(random.uniform(0.,2.*np.pi))
            
        print np.float(i)*100. / (5.*N), '% simulations possible'            
        print np.float(count)*100./np.float(N),'% Useful systems produced \033[2A'

        if len(period_arr) ==  N:
            break
    return period_arr,distance_arr,velocity_arr,primary_arr,companion_arr,phase_arr          
            


        
