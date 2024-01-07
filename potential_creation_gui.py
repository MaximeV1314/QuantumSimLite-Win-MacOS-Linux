import tkinter as tk
from tkinter import ttk
import pandas as pa
import OneD_VP
from PIL import Image, ImageTk
import os
from idlelib.tooltip import Hovertip
from tkinter.messagebox import showinfo

def create_new_potential():

    """
    Frame that enables the creation of a new potential. Call OneD_VP.py function to calculate eigenvalues and eigenvectors
    """

    root2 = tk.Toplevel()
    root2.geometry("340x180")
    root2.resizable(False, False)
    root2.title("New Potential")

    root2.columnconfigure(0, weight=1)
    root2.columnconfigure(1, weight=3)

    image = Image.open(os.path.dirname(__file__) + "\info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    name_new_potential = tk.StringVar()
    ttk.Label(root2, text = "Name :").grid(column=0, row=0, sticky=tk.E, padx = 5, pady = 5)
    name_entry = ttk.Entry(root2, textvariable=name_new_potential)
    name_entry.grid(column=1, row=0, padx = 5, pady = 5, sticky=tk.W)

    new_potential_form = tk.StringVar()
    new_pf_label = ttk.Label(root2, text = "V(X) = ", compound="right", image = img)
    new_pf_label.image = img
    new_pf_label.grid(column=0, row=1, sticky=tk.E, padx = 5, pady = 5)
    potential_entry = ttk.Entry(root2, textvariable=new_potential_form)
    potential_entry.grid(column=1, row=1, padx = 5, pady = 5, sticky=tk.W)

    Tip_new_pf = Hovertip(new_pf_label, "\nPotential formula. \n\nMathematical functions available : \n cos, sin, tan, exp, log, sinh, cosh, tanh, \n \
arcsin, arccos, arctan, arccosh, arcsinh, arctanh, rint, \n\
floor, ceil, pi, sqrt.\n\nUse ** for power operations (ex : x**2 for x squared).\n\
Note that you have to write 2*x for 'two times x' and 2x will not work.\nExample : V(x) = 2*cos(4*x**2)\n")
    
    BC_label = ttk.Label(root2, text = "Boundary\nconditions", image=img, compound="right")
    BC_label.image = img
    BC_label.grid(column=0, row=2, sticky=tk.E, ipadx=5, ipady=5)
    BC_variable = tk.StringVar(value = 0)
    BCD_rb = ttk.Radiobutton(root2, text = "Dirichlet", value = 0, variable=BC_variable)
    BCD_rb.grid(column=1, row=2, padx = 70)
    BCP_rb = ttk.Radiobutton(root2, text = "Periodic", value = 1, variable=BC_variable)
    BCP_rb.grid(column=1, row=2, sticky=tk.W)

    Tip_BC_label = Hovertip(BC_label, "\nIf 2L is the lenght of the box :\n\n\
    Periodic: The wave function ψ at position -L and time t equals ψ \n\
                at position L and time t for all t.\n\n\
    Dirichlet: The wave function ψ at position -L and time t, as well as ψ\n\
                at position L and time t, are both set to zero for all t.\n")

    slider_N_value = tk.DoubleVar()
    slider_N_label = ttk.Label(root2, text = "Discretization\n      level", image=img, compound="right")
    slider_N_label.image = img
    slider_N_label.grid(column=0, row=3, sticky=tk.E, ipadx=5, ipady=5)
    slider_N_value_label = ttk.Label(root2, text = "200")
    slider_N_value_label.grid(column=1, row=3, sticky=tk.W, ipadx=5, ipady=5, padx = 20)
    
    slider_N = ttk.Scale(root2, from_=100, to=600, variable=slider_N_value, 
                          command = lambda event : slider_N_value_label.configure(text = "%.0f" % slider_N_value.get()))
    slider_N.grid(column=1, row=3, padx=5, pady=5)
    slider_N.set(200)

    Tip_slider_N = Hovertip(slider_N_label, "\nRefers to the number of points used to discretize the interval [-L, L], \n\
where 2L represents the lenght of the box (in this case, L = 10).\n\n\
Note that calculation time increases with the number of points\n(particularly for dynamic).\n")

    create_button = ttk.Button(root2, text = "Create", command = lambda : finalisation())
    create_button.grid(row = 4, column=1, padx=5, sticky=tk.W)

    def finalisation():

        for letter in name_new_potential.get():
            if letter in '/:*?"<>|\\':
                showinfo("information", 'Potential name cannot contain the following symbols :\n / : * ? " < > | \\')
                return 0

        df = pa.read_csv('VP/1D/potential_data_name.csv', index_col=0)

        for pot_save in [K for K in pa.read_csv('VP/1D/potential_data_name.csv', index_col=0).columns.values]:
            if pot_save == name_new_potential.get() :
                showinfo(title='Information', message='There is already a potential with that name.')
                return 0

        try :
            OneD_VP.main(new_potential_form.get(), name_new_potential.get(), BC_variable.get(), int(slider_N_value.get()))
        except :
            showinfo("warning", message = "An error has occurred in your potential formula." + \
         " Please refer to the provided information on correctly formatting a potential formula located next to the entry.")
            return 0
        
        df[name_new_potential.get()] = [new_potential_form.get()]
        df.to_csv('VP/1D/potential_data_name.csv')
        os.makedirs("media/" + name_new_potential.get())
        os.makedirs("media/" + name_new_potential.get() + "/dynamic")
        os.makedirs("media/" + name_new_potential.get() + "/Eigenstates")

        root2.quit()
        root2.destroy()


    cancel_button = ttk.Button(root2, text = "Cancel", command = lambda: def_cancel_button())
    cancel_button.grid(row = 4, column=0)

    def def_cancel_button():
        root2.quit()
        root2.destroy()

    root2.protocol("WM_DELETE_WINDOW", def_cancel_button)
    root2.mainloop()
