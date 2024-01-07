import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import pandas as pa
from numpy import cos, sin, tan, exp, log, sinh, cosh, tanh, arcsin, arccos, arctan, arccosh, arcsinh, arctanh, rint, \
                    floor, ceil, pi, sqrt

def main(V_string, V_name, BC, N):

    L = 10
    dx = 2*L/N
    x = np.arange(-L, L, dx)
    X = np.copy(x)

    V = eval(V_string)                                                      #potential array
    if type(V) == type(0) or type(V) == type(0.1):
        V = np.copy(np.ones(len(x)) * V)

    vec = np.ones(N)
    H = (1/dx**2+V) * np.diag(vec) - 1/(2*dx**2) * np.diag(vec[:-1], -1) \
    - 1/(2*dx**2) * np.diag(vec[:-1], 1)                                    #hamiltonian in position representation       

    if int(BC) == 1 :                                                       #periodic hamiltonian
        H[0,-1] = - 1/(2*dx**2)
        H[-1,0] = - 1/(2*dx**2)

    (vp, Tvec_p) = linalg.eigh(H)                      #200 first eigenvalues and eigenstates
    vec_p = np.transpose(Tvec_p)

    df = pa.DataFrame(vec_p, columns = vp)
    df.to_csv('VP/1D/' + V_name + '.csv')       #save into csv files
    df2 = pa.DataFrame(vp)
    df2.to_csv('VP/1D/' + V_name + '_vp.csv')

def dynamique(potential, iw, x0, k0, T, pas, sigma):

    df_ = "VP/1D/" + potential + ".csv"
    df = pa.read_csv(df_, index_col=0)
    vec_p = df.to_numpy()                       #eigenstates

    df_vp_ = "VP/1D/" + potential + "_vp.csv"
    df_vp = pa.read_csv(df_vp_, index_col=0)
    vp = df_vp.to_numpy().flatten()             #eigenvalues

    L = 10
    N = len(vp)
    dx = 2*L/N
    x = np.arange(-L, L, dx, dtype="complex128")
    p = np.linspace(-np.pi/dx, np.pi/dx, len(x), dtype="complex128")
                                                            
    onde_init_n = x*0                     

    if iw == "Gaussian" :
        onde_init_n = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-((x-x0)/(np.sqrt(2)*sigma))**2) * np.exp(-1j * k0 * x)

    elif iw == "Lorentzian" :
        onde_init_n = 1/np.pi * sigma / ((x-x0)**2 + sigma**2) * np.exp(-1j * k0 * x)

    elif iw == "Uniform density" :
        onde_init_n = np.ones(len(x))

    onde_init_n = onde_init_n/np.linalg.norm(onde_init_n)

    an = x*0
    an = np.dot(vec_p, np.conjugate(onde_init_n))                           #wave's coefficients

    psi = [onde_init_n]                                                     #initial wave
    psi_norm = [np.sum(np.conjugate(onde_init_n)*onde_init_n).real]         #norm of the wave
    x_mean = [np.sum(x * np.conjugate(onde_init_n) * onde_init_n).real]
    x_inc = [np.sqrt(np.sum(x*x * np.conjugate(onde_init_n) * onde_init_n) - x_mean[0]**2).real]

    psi_p_0 = p * 0                                                         #initial momentum wave
    psi_p_0 = 1/np.sqrt(2*np.pi) * np.sum(onde_init_n * np.exp(-1j*np.outer(p, x)), axis = 1)
    psi_p = [psi_p_0/np.linalg.norm(psi_p_0)]
    p_mean = [np.sum(p * np.conjugate(psi_p_0) * psi_p_0).real]
    p_inc = [np.sqrt(np.sum(p*p * np.conjugate(psi_p_0) * psi_p_0) - p_mean[0]**2).real]

    t = np.arange(0, T, pas)
    for T in t:
        psi_T = x*0
        psi_p_T = p*0

        psi_T = np.sum(an[:, np.newaxis] * np.exp(-1j * vp[:, np.newaxis] * T) * vec_p, axis = 0)   #wave at T
        psi_p_T = 1/np.sqrt(2*np.pi) * np.sum(psi_T * np.exp(-1j*np.outer(p, x)), axis = 1)         #momentum wave at T
        psi_p_T = psi_p_T/np.linalg.norm(psi_p_T)

        psi.append(psi_T)
        x_mean.append(np.sum(x * np.conjugate(psi_T) * psi_T).real)
        x_inc.append(np.sqrt(np.sum(x*x * np.conjugate(psi_T) * psi_T) - x_mean[-1]**2).real)

        psi_p.append(psi_p_T)
        p_mean.append(np.sum(p * np.conjugate(psi_p_T) * psi_p_T).real)
        p_inc.append(np.sqrt(np.sum(p*p * np.conjugate(psi_p_T) * psi_p_T) - p_mean[-1]**2).real)

        psi_norm.append(np.sum(np.conjugate(psi_T)*psi_T).real)

    return psi, psi_norm, psi_p, t, x_mean, p_mean, x_inc, p_inc