import OneD_VP
import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import askquestion, showinfo
from PIL import Image, ImageTk
import os
import pandas as pa
import potential_creation_gui as pcg
import numpy as np
from idlelib.tooltip import Hovertip
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
import matplotlib.animation as animation
from numpy import cos, sin, tan, exp, log, sinh, cosh, tanh, arcsin, arccos, arctan, arccosh, arcsinh, arctanh, rint, \
                    floor, ceil, pi, sqrt
from time import strftime

matplotlib.use('TkAgg')

def def_new_potential(potential_cb, create):
    if create:
        pcg.create_new_potential()
    potential_list = [K for K in pa.read_csv('VP/1D/potential_data_name.csv', index_col=0).columns.values]

    potential_cb.config(values=potential_list)
    potential_cb.current(len(potential_list)-1)

def create_param_VP(rooot):

    """
    Notebook : Eigenstates part
    """

    def update_vp(event):

        """
        Update graphic when modify parameters
        """

        df_ = "VP/1D/" + potential_cb.get() + ".csv"
        df = pa.read_csv(df_, index_col=0)
        df_vp_ = "VP/1D/" + potential_cb.get() + "_vp.csv"
        df_vp = pa.read_csv(df_vp_, index_col=0)
        vp = df_vp.to_numpy().flatten()

        index = int(spin_box.get())
        if type(event) == type(0):
            index = event
    
        L = 10
        N = len(vp)
        dx = 2*L/N
        x = np.arange(-L, L, dx)
        X = np.copy(x)

        ax.clear()
        ax.grid()
        ax.set_xlim(-10,10)
        ax.plot([-L,L], [0,0], linestyle = "dotted", color = "black")
        ax.set_title("Energy level = %.2f" % np.array(df_vp.iloc[index])[0])

        if potential_checkb_var.get()==True:
            df_pot = pa.read_csv('VP/1D/potential_data_name.csv', index_col=0)
            try :
                V = eval(df_pot[potential_cb.get()][0])
            except :
                V = eval(str(df_pot[potential_cb.get()][0])) * np.ones(len(x))

            if np.max(np.abs(V))==0 :
                ax.plot(x, 0*x, "--", color = "green",label = "Potential")
            else:
                ax.plot(x, V * np.max(np.abs(df.iloc[index])**2)/np.max(np.abs(V)),\
                        '--', color = "green", label = "Potential")
                
        if Wave_part_cb.get() == "Real & Imag":
            ax.plot(x, df.apply(np.real).iloc[index], color = "blue", label = "Real part")
            ax.plot(x, df.apply(np.imag).iloc[index], color = "red", label = "Imaginary part")
        elif Wave_part_cb.get() == "Density of Prob.":
            ax.plot(x, df.iloc[index]**2, color = "black", linewidth = 4, label = "Density of Prob.")
        else :
            ax.plot(x, df.apply(np.real).iloc[index], color = "blue", label = "Real part")
            ax.plot(x, df.apply(np.imag).iloc[index], color = "red", label = "Imaginary part")
            ax.plot(x, np.abs(df.iloc[index])**2, color = "black", linewidth = 4, label = "Density of Prob.")

        ax2.clear()
        ax2.grid()
        ax2.set_xlim(-0.5,20.5)
        ax2.set_ylim(vp[0]-1, vp[20]+1)
        ax2.plot(np.arange(len(vp)), vp, "+", color = "blue")
        ax2.plot(index, vp[index], "+", color = "red", markersize = 10, mew = 5)

        if legend_checkb_var.get() == True :
            ax.legend(loc='upper right')
        
        if type(event) == type(0):
            figure.savefig("media/" + potential_cb.get() + "/Eigenstates/" + str(event) + ".png", bbox_inches='tight', dpi=300)
        else:
            figure_canvas.draw()

    frame = tk.Frame(rooot, bg="white")

    frame.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=3)
    frame.columnconfigure(2, weight=8)

    image = Image.open(os.path.dirname(__file__) + "\info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    ttk.Label(frame, text = "Potential : ", background="white").grid(column=0, row=0, sticky=tk.W, ipadx=5, ipady=5)

    potential_list = [K for K in pa.read_csv('VP/1D/potential_data_name.csv', index_col=0).columns.values]

    potential_cb = ttk.Combobox(frame, values = potential_list)
    potential_cb.grid(column=1, row=0, sticky= tk.W, padx = 5, pady = 5)
    potential_cb.bind("<<ComboboxSelected>>", update_vp)
    potential_cb.current(0)

    new_potential = ttk.Button(frame, text = "New Potential", command = lambda : def_new_potential0())
    new_potential.grid(column=1, row=1, sticky=tk.W, padx = 5, pady = 5)

    def def_new_potential0():
        new_potential.state(['disabled'])
        def_new_potential(potential_cb, True)
        new_potential.state(['!disabled'])

    delete_potential = ttk.Button(frame, text = "Delete", command = lambda : delete())
    delete_potential.grid(column=1, row = 1, sticky=tk.E)

    def delete():
        """
        Function for delete a specific potential
        """
        result = askquestion("Delete", "Are you sure to delete " + potential_cb.get() + \
                             "?\n\nNote that media will be still available in media/" + potential_cb.get() + \
                             " (deleted 'date'). You can still delete this file manually.\n\n" + \
                "Ensure that you are not in the file that will be deleted during the process to avoid errors.", icon='warning')
        if result == 'yes' :

            if potential_cb.get() in {"Box", "Ring", "Harmonic Oscillator", "Barrier", "1D lattice"} :
                showinfo(title='Information', message="You cannot delete the following default potentials : \n" + \
                         "Box, Ring, Harmonic Oscillator, Barrier, 1D lattice.")
            else :
                df_d = pa.read_csv('VP/1D/potential_data_name.csv', index_col=0)
                del df_d[potential_cb.get()]
                df_d.to_csv('VP/1D/potential_data_name.csv')
                os.remove("VP/1D/" + potential_cb.get() + ".csv")
                os.remove("VP/1D/" + potential_cb.get() + "_vp.csv")
                os.rename("media/" + potential_cb.get(), "media/" + potential_cb.get() + " (deleted " + strftime("%Y-%m-%d %H-%M-%S") + ")")
                def_new_potential(potential_cb, False)
                update_vp("")
        
    Eig_label = ttk.Label(frame, text = "Eigenstates :", image = img, compound="right", background="white")
    Eig_label.image = img
    Eig_label.grid(column=0, row=3, sticky=tk.W, ipadx=5, ipady=5)
    spin_box = ttk.Spinbox(frame, from_=0, to=20, command = lambda : update_vp(""))
    spin_box.grid(column=1, row = 3, pady=5)
    spin_box.insert(0, "0")

    Tip_Eig = Hovertip(Eig_label,'Displays the eigenstate (upper graphic) in position representation based on its eigenenergy.' + \
                       '\nThe eigeenergy ​​are arranged in ascending order, from the lowest energy to the twentieth lowest (lower graphic).')

    ttk.Label(frame, text = "Wave type :", background="white").grid(column=0, row=4, sticky=tk.W, ipadx=5, ipady=5)
    Wave_part_cb = ttk.Combobox(frame, values = ["Real & Imag", "Density of Prob.", "All"])
    Wave_part_cb.grid(column=1, row=4, sticky= tk.W, padx = 5, pady = 5)
    Wave_part_cb.current(2)
    Wave_part_cb.bind("<<ComboboxSelected>>", update_vp)

    potential_checkb_var = tk.BooleanVar(value = "True")
    potential_checkb = ttk.Checkbutton(frame, text = "Show Potential", \
                                        variable = potential_checkb_var, command=lambda:update_vp(""))
    potential_checkb.grid(column=1, row=5, sticky= tk.W, padx = 5, pady = 5)

    legend_checkb_var = tk.BooleanVar(value = "True")
    legend_checkb = ttk.Checkbutton(frame, text = "Show legend", \
                                        variable = legend_checkb_var, command=lambda:update_vp(""))
    legend_checkb.grid(column=0, row=5, sticky= tk.E, padx = 5, pady = 5)

    save_eig_button = ttk.Button(frame, text = "Save All", image=img, compound="right", command = lambda : save())
    save_eig_button.image = img
    save_eig_button.grid(column=1, row=7, pady = 10)
    Tip_save_eig = Hovertip(save_eig_button, "Save all eigenstates corresponding to the selected potential as .png files in the media/'potential name'/Eigenstates directory.\n\
Saving all eigenstates may take several seconds, during which the application may experience a temporary freeze. Please, wait until the process is complete.")

    def save():
        """Save eigenstates in media"""

        save_eig_button.state(['disabled'])
        for i in range(0,21):
            update_vp(i)
        save_eig_button.state(['!disabled'])

    frame2 = tk.Frame(frame, bg="white")
    frame2.grid(column=2, row = 0, rowspan=100)
    figure = Figure(figsize=(7, 7), dpi=100)
    figure_canvas = FigureCanvasTkAgg(figure, frame2)
    figure_canvas.get_tk_widget().pack()
    ax = figure.add_axes([0.1,0.45,0.8,0.5])
    ax2 = figure.add_axes([0.1,0.1,0.8,0.3])
    ax2.grid()
    ax.grid()
        
    return frame

def create_param_DYN(rooot):

    """
    Notebook : dynamic part
    """

    frame = tk.Frame(rooot, bg="white")

    frame.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=3)
    frame.columnconfigure(2, weight=7)

    image = Image.open(os.path.dirname(__file__) + "\info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    ttk.Label(frame, text = "Potential : ", background="white").grid(column=0, row=0, sticky=tk.E, ipadx=5, ipady=5)

    potential_list = [K for K in pa.read_csv('VP/1D/potential_data_name.csv', index_col=0).columns.values]

    potential_choice_dyn = tk.StringVar()
    potential_cb_dyn = ttk.Combobox(frame, textvariable = potential_choice_dyn, values = potential_list,\
                                    postcommand = lambda : def_new_potential(potential_cb_dyn, False))
    potential_cb_dyn.grid(column=1, row=0, sticky= tk.W, padx = 5, pady = 5)
    potential_cb_dyn.current(0)

    ttk.Label(frame, text = "Initial wave : ", background="white").grid(column=0, row=3, sticky=tk.E, ipadx=5, ipady=5)

    initial_wave_list = ["Gaussian", "Lorentzian", "Uniform density"]

    initial_wave_choice = tk.StringVar()
    initial_wave_cb = ttk.Combobox(frame, textvariable = initial_wave_choice, values = initial_wave_list)
    initial_wave_cb.grid(column=1, row=3, sticky= tk.W, padx = 5, pady = 5)
    initial_wave_cb.current(0)

    slider_eps_value = tk.DoubleVar(value=2)
    slider_eps_label = ttk.Label(frame, text = "σ = ", image=img, compound="right", background="white")
    slider_eps_label.image = img
    slider_eps_label.grid(column=0, row=4, sticky=tk.E, ipadx=5, ipady=5)
    slider_eps_value_label = ttk.Label(frame, text = "2", background="white")
    slider_eps_value_label.grid(column=1, row=4, sticky=tk.W, ipadx=5, ipady=5)

    Tip_slider_eps = Hovertip(slider_eps_label, "Width at half maximum of the Gaussian/Lorentzian.")
    
    slider_eps = ttk.Scale(frame, from_=0.3, to=4, variable=slider_eps_value, 
                          command = lambda event : slider_eps_value_label.configure(text = "%.1f" % slider_eps_value.get()))
    slider_eps.grid(column=1, row=4, padx=5, pady=5)

    
    slider_x0_value = tk.DoubleVar()
    slider_x0_label = ttk.Label(frame, text = "Initial position : " + str(), background="white")
    slider_x0_label.grid(column=0, row=5, sticky=tk.E, ipadx=5, ipady=5)
    slider_x0_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_x0_value_label.grid(column=1, row=5, sticky=tk.W, ipadx=5, ipady=5)
    
    slider_x0 = ttk.Scale(frame, from_=-10, to=10, variable=slider_x0_value, 
                          command = lambda event : slider_x0_value_label.configure(text = "%.1f" % slider_x0_value.get()))
    slider_x0.grid(column=1, row=5, padx=5, pady=5)


    slider_k0_value = tk.DoubleVar()
    slider_k0_label = ttk.Label(frame, text = "Initial momentum : " + str(), background="white")
    slider_k0_label.grid(column=0, row=6, sticky=tk.E, ipadx=5, ipady=5)
    slider_k0_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_k0_value_label.grid(column=1, row=6, sticky=tk.W, ipadx=5, ipady=5)
    
    slider_k0 = ttk.Scale(frame, from_=-15, to=15, variable=slider_k0_value, 
                          command = lambda event : slider_k0_value_label.configure(text = "%.1f" % slider_k0_value.get()))
    slider_k0.grid(column=1, row=6, padx=5, pady=5)

    slider_time_value = tk.DoubleVar()
    slider_time_label = ttk.Label(frame, text = "Time : ", image = img, compound= "right", background="white")
    slider_time_label.image = img
    slider_time_label.grid(column=0, row=7, sticky=tk.E, ipadx=5, ipady=5)
    slider_time_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_time_value_label.grid(column=1, row=7, sticky=tk.W, ipadx=5, ipady=5)
    
    slider_time = ttk.Scale(frame, from_=5, to=60, variable=slider_time_value,
                          command = lambda event : slider_time_value_label.configure(text = str(slider_time_value.get())[:2]))
    slider_time.grid(column=1, row=7, padx=5, pady=5)
    slider_time.set(15)

    Tip_time = Hovertip(slider_time_label, "One unit of time = 100 frames and is approximately 4 seconds.")

    legend_checkb_var_dyn = tk.BooleanVar(value = "True")
    legend_checkb_dyn = ttk.Checkbutton(frame, text = "Show legend", variable = legend_checkb_var_dyn)
    legend_checkb_dyn.grid(column=1, row=8, sticky= tk.W, padx = 5, pady = 5)

    start_button = ttk.Button(frame, text = "Start", image = img, compound="right", command = lambda:start_dynamique(""))
    start_button.image = img
    start_button.grid(column=0, row=9, columnspan=2, sticky=tk.E, pady = 10)

    Tip_start = Hovertip(start_button, "Initiate dynamic evolution.\n\nThe upper graphic displays the state evolution in position representation, \n\
providing frame count, norm, mean, and standard deviation indicators. \n\n\
The lower graphic illustrates the state evolution in momentum representation, \n\
featuring mean, standard deviation, and Heisenberg inequality indicators.\n\n\
You can pause the animation by clicking on it. \n\n\
The animation may take several seconds, depending on the values of t and the discretization level.")

    save_button = ttk.Button(frame, text = "Save", image=img, compound="right", command = lambda:start_dynamique("save"))
    save_button.image = img
    save_button.grid(column=0, row=9, columnspan=2, padx = 100, pady = 10)

    Tip_save_dyn = Hovertip(save_button, "Save the animation as a mp4 in media/'potential name'/Dynamic directory.\n\
Saving all eigenstates may take several seconds/minutes (depending on the values of t and the discretization level),\n\
during which the application may experience a temporary freeze. Please, wait until the process is complete.")

    frame2 = tk.Frame(frame, bg="white")
    frame2.grid(column=2, row = 0, rowspan=100)
    figure = Figure(figsize=(7, 7), dpi=100)
    figure_canvas = FigureCanvasTkAgg(figure, frame2)
    ax = figure.add_axes([0.1,0.45,0.8,0.5])
    ax2 = figure.add_axes([0.1,0.1,0.8,0.3])
    ax.grid()
    ax2.grid()
    figure_canvas.get_tk_widget().pack()

    def start_dynamique(event):

        """
        Function that initialize and start / save the animation
        """

        ax.cla()
        ax2.cla()

        try :
            psi, psi_norm, psi_p, t, x_mean, p_mean, x_inc, p_inc = \
            OneD_VP.dynamique(potential_choice_dyn.get(), initial_wave_choice.get(), \
            slider_x0_value.get(), slider_k0_value.get(), slider_time_value.get(), 0.01, slider_eps_value.get())
        except :
            showinfo("Information", "This potential does not exist anymore.")
            return 0
        
        L = 10
        N = len(psi[0])
        dx = 2*L/N
        x = np.arange(-L, L, dx)
        X = np.copy(x)
        p = np.linspace(-np.pi/dx, np.pi/dx, len(x))
        
        line_real, = ax.plot([], [], color = "blue", label = "Real part")
        line_imag, = ax.plot([], [], color = "red", label = "Imaginary part")
        line_proba, = ax.plot([], [], linewidth = 4, color = "black", label = "Density of prob.")

        df_pot = pa.read_csv('VP/1D/potential_data_name.csv', index_col=0)
        V = eval(df_pot[potential_cb_dyn.get()][0])
        if np.max(np.abs(V))==0 :
            ax.plot(x, 0*x, "--", color = "green",label = "Potential")
        else:
            ax.plot(x, V * np.max(np.abs(np.array(psi))**2)/np.max(np.abs(V)),\
                    '--', color = "green", label = "Potential")

        if legend_checkb_var_dyn.get() == True:
            ax.legend(loc='lower right')

        ax.set_xlim(-L, L)
        ax.set_ylim(np.min(np.array(psi).real), np.max(np.max(np.array(psi).real)))
        ax.grid()
        ax.plot([-L,L], [0,0], linestyle = "dotted", color = "black")
        box = ax.text(0.12,0.05, "", bbox={'facecolor':'w', 'alpha':0.8, 'pad':5},
                transform=ax.transAxes, ha="center")
        

        line_real_p, = ax2.plot([], [], color = "blue")
        line_imag_p, = ax2.plot([], [], color = "red")
        line_proba_p, = ax2.plot([], [], linewidth = 4, color = "black")

        ax2.set_xlim(-20, 20)   #provient du momentum max qui est de +/- 10
        ax2.set_ylim(np.min(np.array(psi_p).real), np.max(np.array(psi_p).real))
        ax2.grid()
        ax2.plot([-np.pi/dx, np.pi/dx], [0,0], linestyle = "dotted", color = "black")
        box2 = ax2.text(0.13,-0.65, "", bbox={'facecolor':'w', 'alpha':0.8, 'pad':5},
                transform=ax.transAxes, ha="center")
        
        anim_running = True

        def onClick(event):
            """Function that pause the animation when clicking on it"""
            nonlocal anim_running
            if anim_running:
                ani.event_source.stop()
                anim_running = False
            else:
                ani.event_source.start()
                anim_running = True
        
        def animate(i, X, psi, psi_norm, psi_p, p):
                """Animation function"""

                line_real.set_data(X, np.real(psi[i]))
                line_imag.set_data(X, np.imag(psi[i]))
                line_proba.set_data(X, np.real(psi[i] * np.conjugate(psi[i])))
                box.set_text("Frame : " + str(i) + "\nNorm : " + "%.2f" % psi_norm[i] + \
                               "\n<X> = " + "%.2f" % x_mean[i] + "\nΔX = " + "%.2f" % x_inc[i])

                line_real_p.set_data(p, np.real(psi_p[i]))
                line_imag_p.set_data(p, np.imag(psi_p[i]))
                line_proba_p.set_data(p, np.real(psi_p[i] * np.conjugate(psi_p[i])))
                box2.set_text("<P> = " + "%.2f" % p_mean[i] + "\nΔP = " + "%.2f" % p_inc[i] + \
                              "\nΔXΔP = " + "%.2f" % (x_inc[i] * p_inc[i]))

                return line_real, line_imag, line_proba, box, box2, line_real_p, line_imag_p, line_proba_p

        ani = animation.FuncAnimation(figure, animate, frames=len(t), fargs = (x, psi, psi_norm, psi_p, p),
                                        interval= 10, blit=True)
        
        if event == 'save':
            start_button.state(['disabled'])
            save_button.state(['disabled'])
            writervideo = animation.FFMpegWriter(fps=40)
            ani.save('media/' + potential_choice_dyn.get() + '/dynamic/' + strftime("%Y-%m-%d %H-%M-%S") + '.mp4',writer=writervideo, dpi = 150)
            start_button.state(['!disabled'])
            save_button.state(['!disabled'])

        figure.canvas.mpl_connect('button_press_event', onClick)
        figure_canvas.draw()

    return frame

def main():

    """
    Main
    """
    root = tk.Tk()

    notebook_param = ttk.Notebook(root)
    notebook_param.pack()

    root.title('QM Software by @Maximev1314')
    root.geometry("1000x720+50+50")
    root.resizable(False, False)
    root.iconbitmap('maximelogo.ico')

    frame1 = create_param_VP(notebook_param)
    frame1.pack()
    frame2 = create_param_DYN(notebook_param)
    frame2.pack()

    notebook_param.add(frame1, text = "Eigenstates")
    notebook_param.add(frame2, text = "Dynamic")

    root.mainloop()

main()