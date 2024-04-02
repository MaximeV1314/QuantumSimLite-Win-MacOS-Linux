import OneD_VP
import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import askquestion, showinfo
from PIL import Image, ImageTk
import os
import pandas as pa
import potential_creation_gui as pcg
import dynamic_creation_gui as dcg
import numpy as np
from idlelib.tooltip import Hovertip
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
import matplotlib.animation as animation
from numpy import cos, sin, tan, exp, log, sinh, cosh, tanh, arcsin, arccos, arctan, arccosh, arcsinh, arctanh, rint, \
                    floor, ceil, pi, sqrt
from time import strftime
import glob
import threading
import warnings

matplotlib.pyplot.switch_backend('agg')
matplotlib.use('TkAgg')
warnings.filterwarnings("ignore",category=matplotlib.MatplotlibDeprecationWarning)

def potential_cb_update(event, potential_cb, create, dynamic_choice_cb):

    """
    Create a new potential and/or update potential list from comboboxes.
    """

    if create:
        pcg.create_new_potential()

    elif create == False and type(dynamic_choice_cb) != type(0):
        df_pot = pa.read_csv('VP/potential_data_name.csv', index_col=0)
        d = int(df_pot[potential_cb.get()][1])
    
        dynamic_list_ = os.listdir("VP/" + str(d) + "D/dynamic/" + potential_cb.get())

        if dynamic_list_ == []:
            dynamic_list = [""]
            dynamic_choice_cb.config(values = dynamic_list)
            dynamic_choice_cb.current(0)

        else :
            dynamic_list = []
            for a in dynamic_list_ :
                if a[-9:] == "_info.csv" :
                    dynamic_list.append(a[0:-9])
            dynamic_choice_cb.config(values = dynamic_list)
            dynamic_choice_cb.current(len(dynamic_list)-1)

    potential_list = [K for K in pa.read_csv('VP/potential_data_name.csv', index_col=0).columns.values]
    potential_cb.config(values=potential_list)


def create_param_VP(rooot):

    """
    Notebook : Eigenstates part
    """

    def update_vp(event):

        """
        Update graphic when modify parameters
        """

        df_pot = pa.read_csv('VP/potential_data_name.csv', index_col=0)
        d = int(df_pot[potential_cb.get()][1])

        df_ = "VP/" + str(d) + "D/show/" + potential_cb.get() + ".bz2"
        df = pa.read_csv(df_, index_col=0)
        df = df.to_numpy(dtype = np.complex128)
        df_vp_ = "VP/" + str(d) + "D/" + potential_cb.get() + "_vp.csv"
        df_vp = pa.read_csv(df_vp_, index_col=0)
        vp = df_vp.to_numpy().flatten()

        try : 
            index = int(spin_box.get())
        except :
            showinfo("Information", "The value entered in the eigenstate spin box is incorrect." \
                     + "Please enter a number between 0 and 49.")
            index = 0

        if index >= 50 :
            showinfo("Information", "The app displays only the first fifty eigenstates")
            index = 49
            spin_box.delete(0,"end")
            spin_box.insert(49, "49")
            
        if type(event) == type(0):
            index = event

        ax.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()

        L = 10
        
        if d == 1 :

            cmap_cb.state(['disabled'])
            interpo_cb.state(['disabled'])
            
            if type(event) != type(0) or event == 39:
                
                Wave_part_cb.state(['!disabled'])
                potential_checkb.state(['!disabled'])
                legend_checkb.state(['!disabled'])

                Wave_part_cb.state(['readonly'])
                potential_checkb.state(['readonly'])
                legend_checkb.state(['readonly'])

            N = len(vp)
            dx = 2*L/N
            x = np.arange(-L, L, dx)
            X = np.copy(x)

            ax3.set_axis_off()
            ax4.set_axis_off()
            ax5.set_axis_off()

            ax2.grid()
            ax.grid()

            ax.set_xlim(-10,10)
            ax.plot([-L,L], [0,0], linestyle = "dotted", color = "black")
            ax.set_title("Eigenstate")

            if potential_checkb_var.get()==True:
                try :
                    V = eval(df_pot[potential_cb.get()][0])

                except :
                    V = eval(str(df_pot[potential_cb.get()][0])) * np.ones(len(x))

                if np.max(np.abs(V))==0 :
                    ax.plot(x, 0*x, "--", color = "green",label = "Potential")
                else:
                    ax.plot(x, V * np.max(np.abs(df[index])**2)/np.max(np.abs(V)),\
                            '--', color = "green", label = "Potential")
                    
            if Wave_part_cb.get() == "Real & Imag":
                ax.plot(x, np.real(df[index]), color = "blue", label = "Real part")
                ax.plot(x, np.imag(df[index]), color = "red", label = "Imaginary part")
            elif Wave_part_cb.get() == "Density of Prob.":
                ax.plot(x, np.real(df[index] * np.conjugate(df[index])), color = "black", linewidth = 4, label = "Density of Prob.")
            else :
                ax.plot(x, np.real(df[index]), color = "blue", label = "Real part")
                ax.plot(x, np.imag(df[index]), color = "red", label = "Imaginary part")
                ax.plot(x, np.real(df[index] * np.conjugate(df[index])), color = "black", linewidth = 4, label = "Density of Prob.")

            ax2.set_title("Energy level = %.2f" % np.array(df_vp.iloc[index])[0])
            ax2.set_xlim(-0.5,50.5)
            ax2.set_ylim(vp[0]-1, vp[49])
            ax2.plot(np.arange(len(vp)), vp, "+", color = "blue")
            ax2.plot(index, vp[index], "+", color = "red", markersize = 10, mew = 5)

            if legend_checkb_var.get() == True :
                ax.legend(loc='upper right')

        elif d == 2 :

            Wave_part_cb.state(['disabled'])
            potential_checkb.state(['disabled'])
            legend_checkb.state(['disabled'])

            if type(event) != type(0) or event == 39 :
                interpo_cb.state(['!disabled'])
                cmap_cb.state(['!disabled'])
                cmap_cb.state(['readonly'])
                interpo_cb.state(['readonly'])

            L = 10
            N = int(np.sqrt(len(df[0])))
            dx = 2*L/N
            x = np.linspace(-L, L, N)
            y = np.copy(x)
            x, y = np.meshgrid(x, y)
            X, Y = x, y

            ax.set_axis_off()
            ax2.set_axis_off()

            cax = ax3.inset_axes([1.01, 0, 0.05, 1], transform=ax3.transAxes)
            colorbar = figure.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap_cb.get()), ax=ax3,cax=cax,)
            colorbar.set_ticks([])

            cax4 = ax4.inset_axes([-0.08, 0, 0.05, 1], transform=ax4.transAxes)
            colorbar4 = figure.colorbar(matplotlib.cm.ScalarMappable(cmap="viridis"), ax=ax4,cax=cax4)
            colorbar4.set_ticks([])

            ax4.yaxis.tick_right()
            ax5.yaxis.tick_right()

            ax3.set_xticks([0, int(N/4), int(2*N/4), int(3*N/4), N-1], [-10, -5, 0, 5, 10])
            ax3.set_yticks([0, int(N/4), int(2*N/4), int(3*N/4), N-1], [-10, -5, 0, 5, 10])
            ax4.set_xticks([0, int(N/4), int(2*N/4), int(3*N/4), N-1], [-10, -5, 0, 5, 10])
            ax4.set_yticks([0, int(N/4), int(2*N/4), int(3*N/4), N-1], [-10, -5, 0, 5, 10])

            ax3.set_title("Probability density of eigenstate")
            ax4.set_title("Potential")
            ax5.set_title("Energy level = %.2f" % np.array(df_vp.iloc[index])[0])

            try :
                V = eval(df_pot[potential_cb.get()][0]).reshape(N,N)

            except :
                V = eval(str(df_pot[potential_cb.get()][0])) * np.ones((N,N))

            vec_2 = np.real(np.real(df[index].reshape(N, N) * np.conjugate(df[index].reshape(N, N))))
            ax3.imshow(vec_2, cmap = cmap_cb.get(), interpolation = interpo_cb.get(), origin='lower')

            ax4.imshow(V, interpolation="gaussian", origin='lower')

            ax5.grid()
            ax5.set_xlim(-0.5,50.5)
            ax5.set_ylim(vp[0], vp[49])
            ax5.plot(np.arange(len(vp)), vp, "+", color = "blue")
            ax5.plot(index, vp[index], "+", color = "red", markersize = 10, mew = 5)

        
        if type(event) == type(0):
            figure.savefig("media/" + potential_cb.get() + "/Eigenstates/" + str(event) + ".png", bbox_inches='tight', dpi=300)
        else:
            figure_canvas.draw()


    frame = tk.Frame(rooot, bg="white")

    frame.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=3)
    frame.columnconfigure(2, weight=8)

    image = Image.open(os.path.dirname(__file__) + "/info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    ttk.Label(frame, text = "Potential : ", background="white").grid(column=0, row=0, sticky=tk.W, ipadx=5, ipady=5)

    potential_list = [K for K in pa.read_csv('VP/potential_data_name.csv', index_col=0).columns.values]

    potential_cb = ttk.Combobox(frame, values = potential_list)
    potential_cb.grid(column=1, row=0, sticky= tk.W, padx = 5, pady = 5)
    potential_cb.bind("<<ComboboxSelected>>", update_vp)
    potential_cb.current(0)
    potential_cb['state'] = 'readonly'

    new_potential = ttk.Button(frame, text = "New Potential",  image = img, compound="right", command = lambda : def_new_potential_start())
    new_potential.grid(column=1, row=1, sticky=tk.W, ipadx = 6)

    Tip_new_potential = Hovertip(new_potential, "Create a new potential and calculate their eigenstates and eigenenergies.")

    def def_new_potential_start():
        new_potential.state(['disabled'])
        potential_cb_update("", potential_cb, True, 0)
        new_potential.state(['!disabled'])

    delete_potential = ttk.Button(frame, text = "Delete", image=img, compound="right", command = lambda : delete())
    delete_potential.image = img
    delete_potential.grid(column=1, row = 2, sticky=tk.W)

    Tip_Eig = Hovertip(delete_potential,"Delete a potential by choosing it using the potential combobox.")

    def delete():
        """
        Function for delete a specific potential
        """
        result = askquestion("Delete", "Are you sure to delete " + potential_cb.get() + "?", icon='warning')
        if result == 'yes' :

            if potential_cb.get() in {"Box", "Ring", "Harmonic Oscillator", "Barrier", "1D lattice"} :
                showinfo(title='Information', message="You cannot delete the following default potentials : \n" + \
                         "Box, Ring, Harmonic Oscillator, Barrier, 1D lattice.")
            else :

                df_d = pa.read_csv('VP/potential_data_name.csv', index_col=0)
                d = df_d[potential_cb.get()][1]
                del df_d[potential_cb.get()]
                df_d.to_csv('VP/potential_data_name.csv')
                os.remove("VP/" + str(d) + "D/" + potential_cb.get() + ".bz2")
                os.remove("VP/" + str(d) + "D/show/" + potential_cb.get() + ".bz2")
                os.remove("VP/" + str(d) + "D/" + potential_cb.get() + "_vp.csv")

                files = glob.glob("VP/" + str(d) + "D/dynamic/" + potential_cb.get() + "/*")
                for f in files:
                    os.remove(f)
                os.rmdir("VP/" + str(d) + "D/dynamic/" + potential_cb.get())

                result2 = askquestion("Delete", "Do you want to also delete media?\n\nIf no, media will be still available in media/" 
                                      + potential_cb.get() + " (deleted 'date').", icon='warning')
                if result2 == 'yes':

                    files = glob.glob("media/" + potential_cb.get() + "/dynamic/*")
                    for f in files:
                        os.remove(f)

                    files = glob.glob("media/" + potential_cb.get() + "/Eigenstates/*")
                    for f in files:
                        os.remove(f)

                    os.rmdir("media/" + potential_cb.get() + "/dynamic")
                    os.rmdir("media/" + potential_cb.get() + "/Eigenstates")
                    os.rmdir("media/" + potential_cb.get())



                else :
                    os.rename("media/" + potential_cb.get(), "media/" + potential_cb.get() + " (deleted " + strftime("%Y-%m-%d %H-%M-%S") + ")")

                potential_cb_update("", potential_cb, False, 0)
                potential_cb.current(0)
                update_vp("")

    Eig_label = ttk.Label(frame, text = "Eigenstates :", image = img, compound="right", background="white")
    Eig_label.image = img
    Eig_label.grid(column=0, row=3, sticky=tk.W, ipadx=5, ipady=5)
    spin_box = ttk.Spinbox(frame, from_=0, to=49, command = lambda : update_vp(""))
    spin_box.grid(column=1, row = 3, pady=5)
    spin_box.insert(0, "0")

    Tip_Eig = Hovertip(Eig_label,'Displays the first fifty eigenstates (upper graphic) in position representation based on its eigenenergy.' + \
                       '\nThe eigeenergies ​​are arranged in ascending order, from the lowest energy to the twentieth lowest (lower graphic).')

    
    Wave_type_label = ttk.Label(frame, text = "Wave type :", image = img, compound="right", background="white")
    Wave_type_label.image = img
    Wave_type_label.grid(column=0, row=4, sticky=tk.W, ipadx=5, ipady=5)
    Wave_part_cb = ttk.Combobox(frame, values = ["Real & Imag", "Density of Prob.", "All"])
    Wave_part_cb.grid(column=1, row=4, sticky= tk.W, padx = 5, pady = 5)
    Wave_part_cb.current(2)
    Wave_part_cb.bind("<<ComboboxSelected>>", update_vp)
    Wave_part_cb['state'] = 'readonly'

    Tip_wave_type = Hovertip(Wave_type_label, "Only for 1D systems Select your wave type between :\nReal & Imag, Density of Prob., All.")

    potential_checkb_var = tk.BooleanVar(value = "True")
    potential_checkb = ttk.Checkbutton(frame, text = "Show Potential", \
                                        variable = potential_checkb_var, command=lambda:update_vp(""))
    potential_checkb.grid(column=1, row=5, sticky= tk.W, padx = 5)

    legend_checkb_var = tk.BooleanVar(value = "True")
    legend_checkb = ttk.Checkbutton(frame, text = "Show legend", \
                                        variable = legend_checkb_var, command=lambda:update_vp(""))
    legend_checkb.grid(column=1, row=6, sticky= tk.W, padx = 5)

    cmap_label = ttk.Label(frame, text = "Colormap :", image=img, compound="right", background="white")
    cmap_label.image = img
    cmap_label.grid(column=0, row=7, sticky=tk.W, ipadx=5, ipady=5)
    cmap_cb = ttk.Combobox(frame, values = ["magma", "afmhot", "jet", "viridis", "plasma", "gist_gray"])
    cmap_cb.grid(column=1, row=7, sticky= tk.W, padx = 5, pady = 5)
    cmap_cb.current(1)
    cmap_cb.bind("<<ComboboxSelected>>", update_vp)
    cmap_cb['state'] = 'disabled'

    Tip_cmap = Hovertip(cmap_label, "Only for 2D systems. Allow the user to change the color of the plot between :\n\
magma, afmhot, jet, viridis, plasma and gist_gray")

    interpo_label = ttk.Label(frame, text = "Interpolation :", image=img, compound="right",background="white")
    interpo_label.image = img
    interpo_label.grid(column=0, row=8, sticky=tk.W, ipadx=5, ipady=5)
    interpo_cb = ttk.Combobox(frame, values = ["gaussian", "bicubic", "quadric", "hermite", "None"])
    interpo_cb.grid(column=1, row=8, sticky= tk.W, padx = 5, pady = 5)
    interpo_cb.current(0)
    interpo_cb.bind("<<ComboboxSelected>>", update_vp)
    interpo_cb['state'] = 'disabled'

    Tip_interpo = Hovertip(interpo_label, "Only for 2D systems. Allow the user to change the interpolation between pixels\
of the plot, using :\ngaussian, bicubic, quadric, hermite or None.")

    save_eig_button = ttk.Button(frame, text = "Save All", image=img, compound="right", command = lambda : multi_threading_save())
    save_eig_button.image = img
    save_eig_button.grid(column=1, row=9, sticky=tk.W, pady = 5, padx = 5)
    Tip_save_eig = Hovertip(save_eig_button, "Save all eigenstates corresponding to the selected potential as .png files in the media/'potential name'/Eigenstates directory.\n\
Saving all eigenstates may take several seconds, during which the application may experience a temporary freeze. Please, wait until the process is complete.")

    def multi_threading_save():
        tt = threading.Thread(target=save)
        tt.setDaemon(True)
        tt.start()


    def save():
        """Save eigenstates in media"""

        try :
            popup = tk.Toplevel()
        except :
            print("An error occurs. Please restart the app to save the files.")
            return 0

        progress_bar = ttk.Progressbar(popup, orient='horizontal', mode='determinate', length=280, maximum=40)
        progress_bar.grid(column=0, row=0, columnspan=2, padx=10, pady=10)

        value_label = ttk.Label(popup, text="Current progess : 0/40")
        value_label.grid(column=0, row=1, columnspan=2, padx = 10, pady = 10)

        potential_cb.state(['disabled'])
        save_eig_button.state(['disabled'])
        delete_potential.state(['disabled'])
        new_potential.state(['disabled'])
        spin_box.state(['disabled'])
        Wave_part_cb.state(['disabled'])
        potential_checkb.state(['disabled'])
        legend_checkb.state(['disabled'])
        cmap_cb.state(['disabled'])
        interpo_cb.state(['disabled'])

        for i in range(0,40):
            update_vp(i)
            try :
                progress_bar["value"] += 1
                value_label["text"] = "Current progess : " + str(i+1) + "/40"
            except :
                potential_cb.state(['!disabled'])
                save_eig_button.state(['!disabled'])
                delete_potential.state(['!disabled'])
                new_potential.state(['!disabled'])
                spin_box.state(['!disabled'])

                potential_cb.state(['readonly'])
                return 0

        potential_cb.state(['!disabled'])
        save_eig_button.state(['!disabled'])
        delete_potential.state(['!disabled'])
        new_potential.state(['!disabled'])
        spin_box.state(['!disabled'])
        potential_cb.state(['readonly'])

        progress_bar.stop()
        popup.destroy()

    frame2 = tk.Frame(frame, bg="white")
    frame2.grid(column=2, row = 0, rowspan=100)
    figure = Figure(figsize=(12, 8), dpi=100)
    figure_canvas = FigureCanvasTkAgg(figure, frame2)
    figure_canvas.get_tk_widget().pack(side = "top", fill = "both", expand  =True)


    ax = figure.add_axes([0.07,0.45,0.8,0.5])
    ax2 = figure.add_axes([0.07,0.08,0.8,0.3])
    ax.grid()
    ax2.grid()

    ax3 = figure.add_axes([0,0.15,0.8,0.8])
    ax4 = figure.add_axes([0.75,0.7,0.2,0.25])
    ax5 = figure.add_axes([0.75,0.13,0.2,0.5])

    ax3.set_xticks([-10, -5, 0, 5, 10])

    ax3.set_axis_off()
    ax4.set_axis_off()
    ax5.set_axis_off()  
        
    return frame

def create_param_DYN(rooot):

    """
    Notebook : dynamic part
    """

    frame = tk.Frame(rooot, bg="white")

    frame.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=1)
    frame.columnconfigure(2, weight=7)

    image = Image.open(os.path.dirname(__file__) + "/info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    dynamic_list_ = os.listdir("VP/1D/dynamic/Box")
    dynamic_list = []
    for a in dynamic_list_ :
        if a[-9:] == "_info.csv" :
            dynamic_list.append(a[0:-9])

    ttk.Label(frame, text = "Dynamic : ", background="white").grid(column=0, row=1, sticky=tk.E, ipadx=5, ipady=5)
    dynamic_choice = tk.StringVar()
    dynamic_choice_cb = ttk.Combobox(frame, textvariable = dynamic_choice, values = dynamic_list)
    dynamic_choice_cb.grid(column=1, row=1, sticky= tk.W, padx = 5, pady = 5)
    dynamic_choice_cb['state'] = 'readonly'

    potential_label = ttk.Label(frame, text = "Potential : ", background="white").grid(column=0, row=0, sticky=tk.E, ipadx=5, ipady=5)
    potential_list = [K for K in pa.read_csv('VP/potential_data_name.csv', index_col=0).columns.values]


    potential_choice_dyn = tk.StringVar()
    potential_cb_dyn = ttk.Combobox(frame, textvariable = potential_choice_dyn, values = potential_list, postcommand=lambda : potential_cb_update("", potential_cb_dyn,False,0))
    potential_cb_dyn.bind("<<ComboboxSelected>>", lambda event, a = potential_cb_dyn, b = False, c = dynamic_choice_cb : potential_cb_update(event, a,b,c))
    potential_cb_dyn.grid(column=1, row=0, sticky= tk.W, padx = 5, pady = 5)
    potential_cb_dyn.current(0)
    
    potential_cb_dyn['state'] = 'readonly'

    new_dynamic_button = ttk.Button(frame, text = "New Dynamic", image=img, compound="right", command = lambda:start_new_dynamic())
    new_dynamic_button.grid(column=1, row=2, sticky=tk.W, columnspan=2, padx = 10)

    Tip_new_dyn = Hovertip(new_dynamic_button, "Create a new dynamic simulation.")

    def delete_dyn():
        result = askquestion("Delete", "Are you sure to delete " + dynamic_choice_cb.get() + "?", icon='warning')
        if result == 'yes' :
            if dynamic_choice_cb.get() == "":
                return 0
            
            df_d = pa.read_csv('VP/potential_data_name.csv', index_col=0)
            d = df_d[potential_cb_dyn.get()][1]
            os.remove("VP/" + str(d) + "D/dynamic/" + potential_cb_dyn.get() + "/" + dynamic_choice_cb.get() + ".csv")
            os.remove("VP/" + str(d) + "D/dynamic/" + potential_cb_dyn.get() + "/" + dynamic_choice_cb.get() + "_info.csv")
            os.remove("VP/" + str(d) + "D/dynamic/" + potential_cb_dyn.get() + "/" + dynamic_choice_cb.get() + "_momentum.csv")
            potential_cb_update("", potential_cb_dyn, False, dynamic_choice_cb)

    delete_dynamic = ttk.Button(frame, text = "Delete Dynamic", image=img, compound="right", command = lambda : delete_dyn())
    delete_dynamic.grid(column=1, row = 3, sticky=tk.W, padx=10)

    Tip_delete = Hovertip(delete_dynamic,"Delete a dynamic simulaton by choosing it using the above combobox.")

    def start_new_dynamic():
        new_dynamic_button.state(['disabled'])
        dcg.create_new_dynamic()
        new_dynamic_button.state(['!disabled'])
        potential_cb_update("", potential_cb_dyn, False, dynamic_choice_cb)

    legend_checkb_var_dyn = tk.BooleanVar(value = "True")
    legend_checkb_dyn = ttk.Checkbutton(frame, text = "Show legend", variable = legend_checkb_var_dyn)
    legend_checkb_dyn.grid(column=1, row=4, sticky= tk.W, padx = 5, pady = 5)

    cmap_label = ttk.Label(frame, text = "Colormap :", image=img, compound="right", background="white")
    cmap_label.image = img
    cmap_label.grid(column=0, row=5, sticky=tk.W, ipadx=5, ipady=5)
    cmap_cb = ttk.Combobox(frame, values = ["magma", "afmhot", "jet", "viridis", "plasma", "gist_gray"])
    cmap_cb.grid(column=1, row=5, sticky= tk.W, padx = 5, pady = 5)
    cmap_cb.current(1)

    Tip_cmap = Hovertip(cmap_label, "Only for 2D systems. Allow the user to change the color of the plot between :\n\
magma, afmhot, jet, viridis, plasma and gist_gray")

    interpo_label = ttk.Label(frame, text = "Interpolation :", image=img, compound="right", background="white")
    interpo_label.image = img
    interpo_label.grid(column=0, row=6, sticky=tk.W, ipadx=5, ipady=5)
    interpo_cb = ttk.Combobox(frame, values = ["gaussian", "bicubic", "quadric", "hermite", "None"])
    interpo_cb.grid(column=1, row=6, sticky= tk.W, padx = 5, pady = 5)
    interpo_cb.current(0)

    Tip_interpo = Hovertip(interpo_label, "Only for 2D systems. Allow the user to change the interpolation between pixels\
of the plot, using :\ngaussian, bicubic, quadric, hermite or None.")

    start_button = ttk.Button(frame, text = "Start", image = img, compound="right", command = lambda:start_dynamique(""))
    start_button.image = img
    start_button.grid(column=1, row=7, columnspan=2, sticky=tk.W, padx=10)

    Tip_start = Hovertip(start_button, "Initiate dynamic evolution.\n\nIn 1D (resp. in 2D), the upper graphic (resp. the left graphic) \
displays the state evolution in position representation, \n\
providing frame count, norm, mean (resp. the middle right graphic), and standard deviation indicators. \n\n\
The lower graphic (resp. the upper right graphic) illustrates the state evolution in momentum representation, \n\
featuring mean (resp. the bottom right graphic), standard deviation, and Heisenberg inequality indicators.\n\n\
You can pause the animation by clicking on it. \n\n\
The animation may take several seconds, depending on the values of t and the discretization level.")

    save_button = ttk.Button(frame, text = "Save", image=img, compound="right", command = lambda:start_dynamique("save"))
    save_button.image = img
    save_button.grid(column=1, row=8, sticky=tk.W, columnspan=2, padx = 10)

    Tip_save_dyn = Hovertip(save_button, "Save the animation as a mp4 in media/'potential name'/Dynamic directory.\n\
Saving the animation may take several seconds/minutes (depending on the values of t and the discretization level),\n\
during which the application may experience a temporary freeze. Please, wait until the process is complete.")

    frame2 = tk.Frame(frame, bg="white")
    frame2.grid(column=2, row = 0, rowspan=100)
    figure = Figure(figsize=(12, 8), dpi=100)
    figure_canvas = FigureCanvasTkAgg(figure, frame2)

    ax = figure.add_axes([0.07,0.45,0.8,0.5])
    ax2 = figure.add_axes([0.07,0.08,0.8,0.3])
    ax.grid()
    ax2.grid()

    ax3 = figure.add_axes([0,0.12,0.82,0.82])
    ax4 = figure.add_axes([0.71,0.68,0.28,0.28])
    ax5 = figure.add_axes([0.75,0.38,0.2,0.25])
    ax6 = figure.add_axes([0.75,0.06,0.2,0.25])

    ax3.set_xticks([-10, -5, 0, 5, 10])

    ax3.set_axis_off()
    ax4.set_axis_off()
    ax5.set_axis_off() 
    ax6.set_axis_off()  

    figure_canvas.get_tk_widget().pack()

    def start_dynamique(event):

        """
        Function that initialize and start / save the animation
        """

        ax.cla()
        ax2.cla()
        ax3.cla()
        ax4.cla()
        ax5.cla()
        ax6.cla()

        df_pot = pa.read_csv('VP/potential_data_name.csv', index_col=0)
        try : 
            d = int(df_pot[potential_cb_dyn.get()][1])
        except :
            showinfo(title='Information', message="This file doesn't exist anymore.")
            return 0

        df_ = "VP/" + str(d) + "D/dynamic/" + potential_cb_dyn.get() + "/"

        try : 
            psi = (pa.read_csv(df_ + dynamic_choice_cb.get() + ".csv", index_col=0)).to_numpy(dtype = np.complex128)
        except :
            showinfo(title='Information', message="File doesn't exist.")
            return 0
        
        psi_p = (pa.read_csv(df_ + dynamic_choice_cb.get() + "_momentum.csv", index_col=0)).to_numpy(dtype = np.complex128)

        if d == 1 :
            
            x_mean, x_inc, Px_mean, Px_inc, psi_norm = (pa.read_csv(df_ + dynamic_choice_cb.get() + "_info.csv", index_col=0)).to_numpy()

            L = 10
            N = len(psi[0])
            dx = 2*L/N
            x = np.arange(-L, L, dx)
            X = np.copy(x)
            p = np.linspace(-np.pi/dx, np.pi/dx, len(x))
            
            line_real, = ax.plot([], [], color = "blue", label = "Real part")
            line_imag, = ax.plot([], [], color = "red", label = "Imaginary part")
            line_proba, = ax.plot([], [], linewidth = 4, color = "black", label = "Density of prob.")

            V = eval(str(df_pot[potential_cb_dyn.get()][0]))
            if np.max(np.abs(V))==0 :
                ax.plot(x, 0*x, "--", color = "green",label = "Potential")
            else:
                ax.plot(x, V * np.max(np.abs(np.array(psi))**2)/np.max(np.abs(V)),\
                        '--', color = "green", label = "Potential")

            if legend_checkb_var_dyn.get() == True:
                ax.legend(loc='lower right')

            ax3.set_axis_off()
            ax4.set_axis_off()
            ax5.set_axis_off()
            ax6.set_axis_off() 

            ax.set_xlim(-L, L)
            ax.set_ylim(np.min(np.array(psi).real), np.max(np.max(np.array(psi).real)))
            ax.grid()
            ax.plot([-L,L], [0,0], linestyle = "dotted", color = "black")
            ax.set_title("Wave function")
            box = ax.text(0.07,0.05, "", bbox={'facecolor':'w', 'alpha':0.8, 'pad':5},
                    transform=ax.transAxes, ha="center")
            

            line_real_p, = ax2.plot([], [], color = "blue")
            line_imag_p, = ax2.plot([], [], color = "red")
            line_proba_p, = ax2.plot([], [], linewidth = 4, color = "black")

            ax2.set_xlim(p[0], p[-1])
            ax2.set_ylim(np.min(np.array(psi_p).real), np.max(np.array(psi_p).real))
            ax2.grid()
            ax2.plot([-np.pi/dx, np.pi/dx], [0,0], linestyle = "dotted", color = "black")
            ax2.set_title("Wave function of momentum")
            box2 = ax2.text(0.07,-0.65, "", bbox={'facecolor':'w', 'alpha':0.8, 'pad':5},
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
                    box2.set_text("<P> = " + "%.2f" % Px_mean[i] + "\nΔP = " + "%.2f" % Px_inc[i] + \
                                "\nΔXΔP = " + "%.2f" % (x_inc[i] * Px_inc[i]))

                    return line_real, line_imag, line_proba, box, box2, line_real_p, line_imag_p, line_proba_p

            ani = animation.FuncAnimation(figure, animate, frames=len(psi), fargs = (x, psi, psi_norm, psi_p, p),
                                            interval= 25, blit=True)
            
            if event == 'save':
                writervideo = animation.FFMpegWriter(fps=25)
                try :
                    ani.save('media/' + potential_choice_dyn.get() + '/dynamic/' + strftime("%Y-%m-%d %H-%M-%S") + '.mp4',writer=writervideo, dpi = 100)
                except :
                    showinfo("Information", "Saving animation is impossible because you do not have ffmpeg.")

            figure.canvas.mpl_connect('button_press_event', onClick)
            figure_canvas.draw()

        else :

            x_mean, y_mean, x_inc, y_inc, Px_mean, Py_mean, Px_inc, Py_inc, psi_norm = (pa.read_csv(df_ + dynamic_choice_cb.get() + "_info.csv", index_col=0)).to_numpy()

            L = 10
            N = int(np.sqrt(len(psi[0])))
            dx = 2*L/N

            ax.set_axis_off()
            ax2.set_axis_off()

            ax4.set_xticks([])
            ax4.set_yticks([])

            ax5.grid()
            ax6.grid()

            im = ax3.imshow((psi[0].reshape((N,N))*np.conjugate(psi[0].reshape((N,N)))).real, interpolation=interpo_cb.get(), vmin = 0, cmap = cmap_cb.get())
            im_p = ax4.imshow((psi_p[0].reshape((N,N))*np.conjugate(psi_p[0].reshape((N,N)))).real, interpolation=interpo_cb.get(), vmin = 0, cmap = cmap_cb.get())
            mean_plot, = ax5.plot([], [], "b.", markersize = 2)
            mean_plot_pt, = ax5.plot([], [], "b.", markersize = 20)
            mean_p_plot, = ax6.plot([], [], "r.", markersize = 2)
            mean_p_plot_pt, = ax6.plot([], [], "r.", markersize = 20)

            title = ax3.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                    transform=ax.transAxes, ha="center")
            
            ax5.set_xlim(-L,L)
            ax5.set_ylim(-L,L)

            ax3.set_title('Probability density of the wave function')
            ax4.set_title('PD of momentum')
            ax5.set_title('Mean position')
            ax6.set_title('Mean momentum')

            ax6.set_xlim(-2*pi/dx, 2*pi/dx)
            ax6.set_ylim(-2*pi/dx, 2*pi/dx)

            anim_running2D = True

            def onClick2D(event):
                """Function that pause the animation when clicking on it"""
                nonlocal anim_running2D
                if anim_running2D:
                    anim.event_source.stop()
                    anim_running2D = False
                else:
                    anim.event_source.start()
                    anim_running2D = True

            def animate(i):

                psii = (psi[i].reshape((N,N)) * np.conjugate(psi[i].reshape((N,N)))).real
                psii_p = (psi_p[i].reshape((N,N)) * np.conjugate(psi_p[i].reshape((N,N)))).real

                im.set_array(psii)
                im_p.set_array(psii_p)

                mean_plot.set_data([x_mean[:i], y_mean[:i]])
                mean_plot_pt.set_data([x_mean[i], y_mean[i]])

                mean_p_plot.set_data([Px_mean[1:i], Py_mean[1:i]])
                mean_p_plot_pt.set_data([Px_mean[i], Py_mean[i]])

                title.set_text("Frame = " + str(i) + " ; norm = " + str(psi_norm[i]))

                return [im, im_p, mean_plot, title, mean_p_plot, mean_p_plot_pt, mean_plot_pt]

            anim = animation.FuncAnimation(figure, animate, frames=len(psi),
                                            interval= 25, blit=True)        #fargs = (x, psi, psi_norm, psi_p, p)
            
            if event == 'save':
                writervideo = animation.FFMpegWriter(fps=25)
                try :
                    anim.save('media/' + potential_choice_dyn.get() + '/dynamic/' + strftime("%Y-%m-%d %H-%M-%S") + '.mp4',writer=writervideo, dpi = 120)
                except :
                    showinfo("Information", "Saving animation is impossible because you do not have ffmpeg.")

            figure.canvas.mpl_connect('button_press_event', onClick2D)
            figure_canvas.draw()

        # Plot ΔPΔX in function of time

        # plt.rcParams["figure.figsize"] = (12,10)
        # plt.plot(t, 1/2 * np.ones(len(t)), "r-.", linewidth=3)
        # plt.plot(t, np.array(x_inc)[1:]*np.array(p_inc)[1:], color = "blue")
        # plt.grid()
        # plt.xlim(t[0], t[-1])
        # plt.ylim(bottom = 0)
        # plt.xlabel("t")
        # plt.ylabel("ΔPΔX")
        # plt.savefig("heisenberg_incertitude.png", dpi = 200)

    return frame

def main():

    """
    Main
    """
    root = tk.Tk()

    notebook_param = ttk.Notebook(root)
    notebook_param.pack()

    root.title('QuantumSimLite')
    root.geometry("1300x800+10+10")
    icon = ImageTk.PhotoImage(Image.open(os.path.dirname(__file__) + "/maximelogo.ico", mode = "r"))
    root.iconphoto(True, icon)

    frame1 = create_param_VP(notebook_param)
    frame1.pack()
    frame2 = create_param_DYN(notebook_param)
    frame2.pack()

    notebook_param.add(frame1, text = "Eigenstates")
    notebook_param.add(frame2, text = "Dynamic")

    root.mainloop()

main()

#######################################################
#######   QuantumSimLite made by @maximev131   ########
####    Reach me out at maxime.vinteler@yahoo.fr   ####
#######################################################