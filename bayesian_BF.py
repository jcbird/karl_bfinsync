""" bayesian_BF

Module to compute BF using Bayesian Inference

Functions:
    read_pdetect

"""
import numpy as np


def prior(BF):
    """
    Computes Prior of BF calculation according to Equation 12 of Clark+ (2012).
    From Clark, based on Sivia & Skilling 2006:
    Prior is the 'Jeffrey's Prior'; appropriate because BF is scale parameter.
    Prior is uniform in log BF.


    Parameters
    ----------------

    BF : float in [0,1]
        Binary Fraction

    Returns
    ----------------
    P(BF | B) where B is 'background information'

    """
    return 1./BF


def likelihood_D(BF, pdetect, binarity):
    """
    Likelihood function P(D| BF, B) where BF is binary fraction, D is the
    set of stars detected as binaries and B is background information.

    See equations 7 and 12 of Clark+ (2012)

    Note: This function only computes the P(D| BF, B) given ONE set of detected
    binaries (D). If we break up the sample into sets according to time
    baseline, the P(D | BF, B) becomes a sum over all likelihood_D.

    Parameters
    ----------------

    BF : float in [0,1]
        Binary fraction
    pdetect : array of length N
        Probability that each star could be detected as binary if it was one.
        This is the result of the sumulations. pdectect_i must be in [0,1].
    binarity : array of length N
        Representation of whether or not each star was detected as binary.
        Assumed 1 for detection; 0 for non-detection.

    Returns
    ----------------
    P(D | BF, B) as float

    """
    prob_binary = BF * pdetect + (1. - BF) * 1E-3
    isbinary = binarity == 1.0
    binary_term = np.prod(prob_binary[isbinary])
    nonbinary_term = np.prod(1. - prob_binary[~isbinary])
    return binary_term * nonbinary_term


def prob_BF(likelihood, prior):
    """
    Compute P( BF | D,B) according to equation 12 of Clark+ (2012).
    Note, if we break up sample, likelood needs to be a sum over all
    likelihood_Ds.

    Parameters
    ----------------

    likelihood : float
        P(D | BF, B) from some combinations of likelihood_D
    prior : float
        P(BF | B) from prior function

    Returns
    ----------------
    P(BF | D, B) as float

    """
    return likelihood * prior


def normalized_pBF(pdetect, binarity, BFsteps=1000):
    """
    Returns array of normalized probabilities given a set of BFs and
    their associated prob_BF.
    Normalization is fromn equation 13.
    int_0^1 P(BF | D, B) dBF = 1

    Calculates probability of BF in [0,1] in BFsteps using
    prior, likelihood_D, prob_BF functions.

    Parameters
    ----------------

    pdetect : array of length N
        Probability that each star could be detected as binary if it was one.
        This is the result of the sumulations. pdectect_i must be in [0,1].
    binarity : array of length N
        Representation of whether or not each star was detected as binary.
        Assumed 1 for detection; 0 for non-detection.
    BFsteps : int
        Number of steps to take in BF range from 0 to 1

    Returns
    ----------------
    BFs, prob_BFs arrays
    BFs is assumed BF.
    prob_BFs is associated normalized P(BF | D, B).
    """

    init_BF = 1.0/BFsteps
    BFs = np.linspace(init_BF, 1.0, BFsteps)
    dBF = np.diff(BFs)[0]
    lh_bfs = np.array([likelihood_D(BF, pdetect, binarity) for BF in BFs])
    p_BFs = np.array([prob_BF(lh, prior(BF)) for lh, BF in zip(lh_bfs, BFs)])
    normalization = np.nansum(p_BFs * dBF)  # equation 13 Clark
    return BFs, p_BFs / normalization
