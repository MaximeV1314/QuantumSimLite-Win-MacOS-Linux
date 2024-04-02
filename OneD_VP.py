import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import pandas as pa
from numpy import cos, sin, tan, exp, log, sinh, cosh, tanh, arcsin, arccos, arctan, arccosh, arcsinh, arctanh, rint, \
                    floor, ceil, pi, sqrt
import matplotlib.animation as animation
import tkinter  as tk
from tkinter import ttk
import sys

def main(V_string, V_name, BC, N, d):

    L = 10
    dx = dy = 2*L/N
    x = np.linspace(-L, L, N)
    X = np.copy(x)

    if d == 1 :

        try :
            V = eval(V_string)
        except :
            return 0
        
        if type(V) == type(0) or type(V) == type(0.1):
            V = np.copy(np.ones(len(x)) * V)

        V[V == np.Inf] = np.exp(650)
        V[V == -np.Inf] = -np.exp(650)
        
        vec = np.ones(N)    
        H = (1/dx**2 + V) * np.diag(vec) + (-1/(2*dx**2)) * np.diag(vec[:-1], -1) \
        + (-1/(2*dx**2)) * np.diag(vec[:-1], 1)                                 # hamiltonian in position representation   

        if int(BC) == 1 :                                                       #periodic hamiltonian
            H[0,-1] = - 1/(2*dx**2)
            H[-1,0] = - 1/(2*dx**2)

        (vp, Tvec_p) = linalg.eigh(H)                      #200 first eigenvalues and eigenstates
        vec_p = np.transpose(Tvec_p)

        limite = None

    elif d == 2 :
        #y = np.arange(-L, L, dx)
        y = np.linspace(-L,L, N)
        Y = np.copy(x)

        x, y = np.meshgrid(x, y)
        X, Y = x, y

        try :
            V = eval(V_string)
        except :
            return 0

        if type(V) == type(0) or type(V) == type(0.1):
            V = np.ones((N,N)) * V

        V[V == np.Inf] = np.exp(650)
        V[V == -np.Inf] = -np.exp(650)

        vec = np.ones(N**2)
        V_f = V.flatten()

        H = (1/dx**2 + 1/dy**2 + V_f) * np.diag(vec) - 1/(2*dx**2) * np.diag(vec[:-1], -1) - 1/(2*dx**2) * np.diag(vec[:-1], 1) -\
        1/(2*dy**2) * np.diag(vec[:-N], N) - 1/(2*dy**2) * np.diag(vec[:-N], -N)

        for i in range(N, N**2, N):
            H[i, i-1] = 0
            H[i-1, i] = 0
        
        (vp, Tvec_p) = linalg.eigh(H)                      
        vec_p = np.transpose(Tvec_p)

        limite = None

    df = pa.DataFrame(vec_p[0:limite])
    df.to_csv("VP/" + str(d) + "D/" + V_name + ".bz2")       #save into csv files
    df_50 = pa.DataFrame(vec_p[0:50])
    df_50.to_csv("VP/" + str(d) + "D/show/" + V_name + ".bz2")       #save into csv files
    df2 = pa.DataFrame(vp[0:limite])
    df2.to_csv("VP/" + str(d) + "D/" + V_name + "_vp.csv")

    return 1


def dynamique(potential, iw, x0, kx0, T, pas, sigma, y0, ky0, name):

    popup = tk.Toplevel()

    progress_bar = ttk.Progressbar(popup, orient='horizontal', mode='determinate', length=280, maximum=T)
    progress_bar.grid(column=0, row=0, columnspan=2, padx=10, pady=10)

    value_label = ttk.Label(popup, text="Reading eigenstates...")
    value_label.grid(column=0, row=1, columnspan=2, padx = 10, pady = 10)

    df_pot = pa.read_csv('VP/potential_data_name.csv', index_col=0)
    d = int(df_pot[potential][1])

    df_ = "VP/" + str(d) + "D/" + potential + ".bz2"
    vec_p = pa.read_csv(df_, index_col=0).to_numpy(dtype=np.complex128) 

    df_vp_ = "VP/" + str(d) + "D/" + potential + "_vp.csv"
    df_vp = pa.read_csv(df_vp_, index_col=0)
    vp = df_vp.to_numpy().flatten()             #eigenvalues

    try : 
        value_label['text'] = "Simulation calculation\nCurrent progress : 0.0%"
    except :
        0

    if d == 1 :

        L = 10
        N = len(vp)
        dx = 2*L/N
        x = X = np.linspace(-L, L, N, dtype="complex128")
        Px = np.linspace(-np.pi/dx, np.pi/dx, len(x), dtype="complex128")
                                                                                   
        if iw == "Gaussian" :
            onde_init_n = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-((x-x0)/(np.sqrt(2)*sigma))**2) * np.exp(-1j * kx0 * x)

        elif iw == "Lorentzian" :
            onde_init_n = 1/np.pi * sigma / ((x-x0)**2 + sigma**2) * np.exp(-1j * kx0 * x)

        elif iw == "Uniform density" :
            onde_init_n = np.ones(len(x))

        onde_init_n = onde_init_n/np.linalg.norm(onde_init_n)
        psi_p_0 = 1/np.sqrt(2*np.pi) * np.sum(onde_init_n * np.exp(-1j*np.outer(Px, x)), axis = 1)   #initial momentum wave

        RS = -1

        x_mean = [np.sum(X * np.conjugate(onde_init_n) * onde_init_n).real]
        x_inc = [np.sqrt(np.sum(X*X * np.conjugate(onde_init_n) * onde_init_n) - x_mean[0]**2).real]
                                                    
        psi_p = [(psi_p_0/np.linalg.norm(psi_p_0))]
        Px_mean = [np.sum(Px * np.conjugate(psi_p_0) * psi_p_0).real]
        Px_inc = [np.sqrt(np.sum(Px*Px * np.conjugate(psi_p_0) * psi_p_0) - Px_mean[0]**2).real]
        
        y_mean = Py_mean = y_inc = Py_inc = []

    elif d == 2 :

        L = 10
        N = int(np.sqrt(len(vec_p[0])))
        dx = 2*L/N
        x = y = np.linspace(-L, L, N)
        X, Y = np.meshgrid(x, y)

        px = py = np.linspace(-np.pi/dx, np.pi/dx, len(x), dtype="complex128")
        Px, Py = np.meshgrid(px, py)

        if iw == "Gaussian" :
            onde_init_n = ( 2 * np.pi * sigma**2)**(-1/4) * np.exp(-((X-x0)/(np.sqrt(2)*sigma))**2) * np.exp(-1j * kx0 * X) * \
                                                np.exp(-((Y-y0)/(np.sqrt(2)*sigma))**2) * np.exp(-1j * ky0 * Y)

        elif iw == "Lorentzian" :
            onde_init_n = 1/np.pi**2 * sigma / ((X-x0)**2 + sigma**2) * np.exp(-1j * kx0 * X) * \
                                    sigma / ((Y-y0)**2 + sigma**2) * np.exp(-1j * ky0 * Y)

        elif iw == "Uniform density" :
            onde_init_n = np.ones((N,N))

        onde_init_n = (onde_init_n/np.linalg.norm(onde_init_n)).flatten()
        psi_p_0 = 1/np.sqrt(2*np.pi) * np.sum(onde_init_n * np.exp(-1j*np.outer(Px, X)-1j*np.outer(Py, Y)), axis = 1)

        X = X.flatten()
        Y = Y.flatten()
        Px = Px.flatten()
        Py = Py.flatten()

        x_mean = [np.sum(X * np.conjugate(onde_init_n) * onde_init_n).real]
        x_inc = [np.sqrt(np.sum(X*X * np.conjugate(onde_init_n) * onde_init_n) - x_mean[0]**2).real]
        y_mean = [np.sum(Y * np.conjugate(onde_init_n) * onde_init_n).real]
        y_inc = [np.sqrt(np.sum(Y*Y * np.conjugate(onde_init_n) * onde_init_n) - y_mean[0]**2).real]
                                                    
        psi_p = [(psi_p_0/np.linalg.norm(psi_p_0))]
        Px_mean = [np.sum(Px * np.conjugate(psi_p_0) * psi_p_0).real]
        Px_inc = [np.sqrt(np.sum(Px*Px * np.conjugate(psi_p_0) * psi_p_0) - Px_mean[0]**2).real]
        Py_mean = [np.sum(Py * np.conjugate(psi_p_0) * psi_p_0).real]
        Py_inc = [np.sqrt(np.sum(Py*Py * np.conjugate(psi_p_0) * psi_p_0) - Py_mean[0]**2).real]

    an = np.dot(vec_p, np.conjugate(onde_init_n))                           #wave's coefficients

    psi = [onde_init_n]                                                     #initial wave
    psi_norm = [np.sum(np.conjugate(onde_init_n)*onde_init_n).real]         #norm of the wave

    t = np.arange(0, T, pas)

    if d==1 :
        for t_i in t:

            try: 
                progress_bar["value"] += pas
            except :
                return 0
            value_label["text"] = "Simulation calculation\nCurrent progess : %.1f" % (t_i/T * 100) + "%"

            psi_T = np.sum(an[:, np.newaxis] * np.exp(-1j * vp[:, np.newaxis] * t_i) * vec_p, axis = 0)   #wave at T
            psi_p_T = 1/np.sqrt(2*np.pi) * np.sum(psi_T * np.exp(-1j*np.outer(Px, x)), axis = 1)         #momentum wave at T
            psi_p_T = psi_p_T/np.linalg.norm(psi_p_T)

            psi.append(psi_T)
            x_mean.append(np.sum(X * np.conjugate(psi_T) * psi_T).real)
            x_inc.append(np.sqrt(np.sum(X*X * np.conjugate(psi_T) * psi_T) - x_mean[-1]**2).real)

            psi_p.append(psi_p_T)
            Px_mean.append(np.sum(Px * np.conjugate(psi_p_T) * psi_p_T).real)
            Px_inc.append(np.sqrt(np.sum(Px*Px * np.conjugate(psi_p_T) * psi_p_T) - Px_mean[-1]**2).real)

            psi_norm.append(np.sum(np.conjugate(psi_T)*psi_T).real)

        pa.DataFrame(psi).to_csv("VP/" + str(d) + "D/dynamic/" + potential + "/" + name + ".csv")
        pa.DataFrame(psi_p).to_csv("VP/" + str(d) + "D/dynamic/" + potential + "/" + name + "_momentum.csv")
        pa.DataFrame([x_mean, x_inc, Px_mean, Px_inc, psi_norm]).to_csv("VP/" + str(d) + "D/dynamic/" + potential + "/" + name + "_info.csv")

    
    else :
        for t_i in t :

            try: 
                progress_bar["value"] += pas
            except :
                return 0
            value_label["text"] = "Simulation calculation\nCurrent progess : %.1f" % (t_i/T * 100) + "%"

            psi_T = np.sum(an[:, np.newaxis] * np.exp(-1j * vp[:, np.newaxis] * t_i) * vec_p, axis = 0)   #wave at T
            psi_p_T = 1/(2*np.pi) * np.sum(psi_T * np.exp(-1j*np.outer(Px, X)-1j*np.outer(Py, Y)), axis = 1)         #momentum wave at T
            psi_p_T = psi_p_T/np.linalg.norm(psi_p_T)   

            psi.append(psi_T)
            x_mean.append(np.sum(X * np.conjugate(psi_T) * psi_T).real)
            x_inc.append(np.sqrt(np.sum(X*X * np.conjugate(psi_T) * psi_T) - x_mean[-1]**2).real)
            y_mean.append(np.sum(Y * np.conjugate(psi_T) * psi_T).real)
            y_inc.append(np.sqrt(np.sum(Y*Y * np.conjugate(psi_T) * psi_T) - y_mean[-1]**2).real)

            psi_p.append(psi_p_T)
            Px_mean.append(np.sum(Px * np.conjugate(psi_p_T) * psi_p_T).real)
            Px_inc.append(np.sqrt(np.sum(Px*Px * np.conjugate(psi_p_T) * psi_p_T) - Px_mean[-1]**2).real)
            Py_mean.append(np.sum(Py * np.conjugate(psi_p_T) * psi_p_T).real)
            Py_inc.append(np.sqrt(np.sum(Py*Py * np.conjugate(psi_p_T) * psi_p_T) - Py_mean[-1]**2).real)

            psi_norm.append(np.sum(np.conjugate(psi_T)*psi_T).real)

        pa.DataFrame(psi).to_csv("VP/" + str(d) + "D/dynamic/" + potential + "/" + name + ".csv")
        pa.DataFrame(psi_p).to_csv("VP/" + str(d) + "D/dynamic/" + potential + "/" + name + "_momentum.csv")
        pa.DataFrame([x_mean, y_mean, x_inc, y_inc, Px_mean, Py_mean, Px_inc, Py_inc, psi_norm]).to_csv("VP/" + str(d) + "D/dynamic/" + potential + "/" + name + "_info.csv")

    popup.destroy()

#######################################################
#######   QuantumSimLite made by @maximev131   ########
####    Reach me out at maxime.vinteler@yahoo.fr   ####
#######################################################