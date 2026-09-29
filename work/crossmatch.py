import numpy as np

def haversine_distance(mu_alphaO, mu_deltaO, mu_alphaX, mu_deltaX):
    '''
    This function computes the angular distance between two vectors on the unit sphere, via a 
    trigonometric identity that avoids numerical problems.
    '''
    mu_alphaO_rad = np.deg2rad(mu_alphaO)
    mu_deltaO_rad = np.deg2rad(mu_deltaO)
    mu_alphaX_rad = np.deg2rad(mu_alphaX)
    mu_deltaX_rad = np.deg2rad(mu_deltaX)
    Delta_mu_alpha_rad = mu_alphaO_rad - mu_alphaX_rad
    Delta_mu_delta_rad = mu_deltaO_rad - mu_deltaX_rad
    a = np.sin(Delta_mu_delta_rad / 2)**2 + np.cos(mu_deltaO_rad) * np.cos(mu_deltaX_rad) * np.sin(Delta_mu_alpha_rad / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return np.rad2deg(c) * 3600 # units: arcsec

def matrix_rotation(sigma_xprime, sigma_yprime, critical_value, beta):
    '''
    This function performs an ellipse rotation in its matrix representation of an angle beta.
    sigma_xprime and sigma_yprime are the first and second diagonal components of the
    covariance matrix before rotation (i.e. along x' and y'), could be equal when errors
    are circles, and must be given in arcsec.
    Before the rotation, the rescaling must be performed, if necessary. If not, set
    critical_value to 1. beta must be in decimal degrees.
    '''
    beta = np.deg2rad(beta)
    sigma_xprime = sigma_xprime / critical_value**0.5
    sigma_yprime = sigma_yprime / critical_value**0.5
    a11 = sigma_xprime**2 * np.cos(beta)**2 + sigma_yprime**2 * np.sin(beta)**2
    a12 = (sigma_yprime**2 - sigma_xprime**2) * np.cos(beta) * np.sin(beta)
    a21 = a12
    a22 = sigma_xprime**2 * np.sin(beta)**2 + sigma_yprime**2 * np.cos(beta)**2
    return a11, a12, a21, a22

def Mahalanobis_distance_squared(O11, O12, O21, O22, X11, X12, X21, X22, mu_alphaO, mu_deltaO, mu_alphaX, mu_deltaX):
    '''
    This function computes the Mahalanobis distance squared, dM2, given the elements of the 
    covariance matrices of the sources in catalogues O and X and the mean values of the positions
    in equatorial celestial coordinates given in decimal degrees.
    '''
    dH = haversine_distance(mu_alphaO, mu_deltaO, mu_alphaX, mu_deltaX)
    if dH == 0:
        return 0.0
    Omega_X, Omega_O = Omega_function(mu_alphaO, mu_deltaO, mu_alphaX, mu_deltaX, formula="AppendixB")
    Omega_X = np.deg2rad(Omega_X)
    Omega_O = np.deg2rad(Omega_O)
    cos_Omega_X = np.cos(Omega_X)
    sin_Omega_X = np.sin(Omega_X)
    cos_Omega_O = np.cos(Omega_O)
    sin_Omega_O = np.sin(Omega_O)
    X11_rot = X11 * cos_Omega_X**2 + (X12 + X21) * cos_Omega_X * sin_Omega_X + X22 * sin_Omega_X**2
    X12_rot = (X22 - X11) * cos_Omega_X * sin_Omega_X + X12 * cos_Omega_X**2 - X21 * sin_Omega_X**2
    X21_rot = (X22 - X11) * cos_Omega_X * sin_Omega_X + X21 * cos_Omega_X**2 - X12 * sin_Omega_X**2
    X22_rot = X11 * sin_Omega_X**2 - (X12 + X21) * cos_Omega_X * sin_Omega_X + X22 * cos_Omega_X**2
    O11_rot = O11 * cos_Omega_O**2 + (O12 + O21) * cos_Omega_O * sin_Omega_O + O22 * sin_Omega_O**2
    O12_rot = (O22 - O11) * cos_Omega_O * sin_Omega_O + O12 * cos_Omega_O**2 - O21 * sin_Omega_O**2
    O21_rot = (O22 - O11) * cos_Omega_O * sin_Omega_O + O21 * cos_Omega_O**2 - O12 * sin_Omega_O**2
    O22_rot = O11 * sin_Omega_O**2 - (O12 + O21) * cos_Omega_O * sin_Omega_O + O22 * cos_Omega_O**2
    cov11 = O11_rot + X11_rot
    cov12 = O12_rot + X12_rot
    cov21 = O21_rot + X21_rot
    cov22 = O22_rot + X22_rot
    det_cov = cov11 * cov22 - cov12 * cov21
    inv_cov11 = cov22 / det_cov
    dM2 = dH**2 * inv_cov11
    return dM2
    
def Omega_function(mu_alphaO, mu_deltaO, mu_alphaX, mu_deltaX, formula="AppendixB"):
    '''
    This function computes the angle Omega needed to perform the rotation of the covariance matrix
    of the X-ray source due to high values of declination, where the spheric geometry becomes
    non negligible.
    All arguments must be in decimal degrees.
    Check Appendix A of Pineau et al. 2011.
    '''
    n_P11 = np.deg2rad(haversine_distance(mu_alphaO, mu_deltaO, mu_alphaX, mu_deltaX) / 3600) # decimal degrees to radians
    N_P11 = np.deg2rad(mu_alphaO - mu_alphaX) # decimal degrees to radians
    x_P11 = np.deg2rad(90 - mu_deltaO) # decimal degrees to radians
    o_P11 = np.deg2rad(90 - mu_deltaX) # decimal degrees to radians
    s_P11 = (n_P11 + x_P11 + o_P11)/2
    arg_cos = np.sin(s_P11)*np.sin(s_P11 - x_P11)/(np.sin(n_P11)*np.sin(o_P11))
    arg_cos_O = np.sin(s_P11)*np.sin(s_P11 - o_P11)/(np.sin(n_P11)*np.sin(x_P11))
    """
    When n_P11 (angular separation) is very small, the spherical triangle
    is nearly degenerate and s_P11 - x_P11 -> 0 mathematically, but floating
    point rounding can yield a small negative residual (~-1e-16) instead of
    exact zero. This makes arg_cos slightly negative, which has no physical meaning. 
    clip() corrects this numerical noise without affecting the physical result:
    """
    cosX_2 = np.clip(arg_cos, 0.0, None)**0.5
    cosO_2 = np.clip(arg_cos_O, 0.0, None)**0.5
    
    if formula == "AppendixB":
        arg = np.sin(x_P11)*np.sin(N_P11)/np.sin(n_P11) # modified formula
        arg_O = np.sin(o_P11)*np.sin(N_P11)/np.sin(n_P11)
        """
        Same issue as above: for very small n_P11, arg can slightly exceed
        arcsin's [-1, 1] domain due to floating point rounding
        (observed residuals ~1e-13). 
        clip() corrects this without affecting the physical result:
        """
        arcsin = np.arcsin(np.clip(arg, -1.0, 1.0))
        arcsin_O = np.arcsin(np.clip(arg_O, -1.0, 1.0))
    elif formula == "A6":
        arg = x_P11*np.sin(N_P11)/n_P11 # original formula from Pineau et al. 2011
        arg_O = o_P11*np.sin(N_P11)/n_P11
        if -1 <= arg <= 1 and -1 <= arg_O <= 1:
            arcsin = np.arcsin(arg)
            arcsin_O = np.arcsin(arg_O)
            print("A6 succeded at arcsin({}), ".format(arg), "in degrees: {}°".format(arg*180/np.pi))
        else: 
            print("A6 failed at arcsin({})".format(arg), "in degrees: {}°".format(arg*180/np.pi))
            return
    else:
        raise ValueError(formula)
  
    condition1 = cosX_2 > 2**0.5/2
    condition2 = arcsin >= 0
    condition1_O = cosO_2 > 2**0.5/2
    condition2_O = arcsin_O >= 0
    if condition1 and condition2:
        X_P11 = np.pi - arcsin
    elif condition1 and not condition2:
        X_P11 = -np.pi - arcsin
    else:
        X_P11 = arcsin
    if condition1_O and condition2_O:
        O_P11 = np.pi - arcsin_O
    elif condition1_O and not condition2_O:
        O_P11 = -np.pi - arcsin_O
    else:
        O_P11 = arcsin_O
    if X_P11 >= 0:
        Omega_X = np.pi/2 - X_P11
    else:
        Omega_X = -np.pi/2 - X_P11
    if O_P11 >= 0:
        Omega_O = np.pi/2 - O_P11
    else:
        Omega_O = -np.pi/2 - O_P11
    return np.rad2deg(Omega_X), np.rad2deg(Omega_O)
