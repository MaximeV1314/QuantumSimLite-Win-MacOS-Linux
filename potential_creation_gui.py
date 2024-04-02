import tkinter as tk
from tkinter import ttk
import pandas as pa
import OneD_VP
from PIL import Image, ImageTk
import os
from idlelib.tooltip import Hovertip
from tkinter.messagebox import showinfo
import threading
import sys

def create_new_potential():

    """
    Frame that enables the creation of a new potential. Call OneD_VP.py function to calculate eigenvalues and eigenvectors
    """

    root2 = tk.Toplevel()
    root2.geometry("340x260")
    root2.resizable(False, False)
    root2.title("New Potential")

    root2.columnconfigure(0, weight=1)
    root2.columnconfigure(1, weight=3)

    image = Image.open(os.path.dirname(__file__) + "/info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    name_new_potential = tk.StringVar()
    ttk.Label(root2, text = "Name :").grid(column=0, row=0, sticky=tk.E, padx = 5, pady = 5)
    name_entry = ttk.Entry(root2, textvariable=name_new_potential)
    name_entry.grid(column=1, row=0, padx = 5, pady = 5, sticky=tk.W)

    dimension_label = ttk.Label(root2, text = "Dimension : ", compound="right", image = img)
    dimension_label.image = img
    dimension_label.grid(column=0, row=1, sticky=tk.E, ipadx=5, ipady=5)
    dimension_variable = tk.StringVar(value = 1)
    dimension_1D = ttk.Radiobutton(root2, text = "1D", value = 1, variable=dimension_variable, command=lambda : RB())
    dimension_1D.grid(column=1, row=1, sticky=tk.W)
    dimension_2D = ttk.Radiobutton(root2, text = "2D", value = 2, variable=dimension_variable, command=lambda : RB())
    dimension_2D.grid(column=1, row=1, sticky=tk.W, padx = 50)

    Tip_dimension = Hovertip(dimension_label, "Space dimension of the system ")

    def RB():
        if int(dimension_variable.get()) == 1:
            new_pf_label['text'] = "V(X) = "
            slider_N["from_"] = 100
            slider_N["to"] = 500
            slider_N.set(200)
            BCD_rb['state'] = "!disabled"
            BCP_rb['state'] = "!disabled"

        else :
            new_pf_label['text'] = "V(X,Y) = "
            slider_N["from_"] = 40
            slider_N["to"] = 100
            slider_N.set(60)
            BCD_rb['value'] = 0
            BCP_rb['value'] = 1
            BCD_rb['state'] = "disabled"
            BCP_rb['state'] = "disabled"


    new_potential_form = tk.StringVar()
    new_pf_label = ttk.Label(root2, text = "V(X) = ", compound="right", image = img)
    new_pf_label.image = img
    new_pf_label.grid(column=0, row=2, sticky=tk.E, padx = 5, pady = 5)
    potential_entry = ttk.Entry(root2, textvariable=new_potential_form)
    potential_entry.grid(column=1, row=2, padx = 5, pady = 5, sticky=tk.W)

    Tip_new_pf = Hovertip(new_pf_label, "\nPotential formula. \n\nMathematical functions available : \n cos, sin, tan, exp, log, sinh, cosh, tanh, \n \
arcsin, arccos, arctan, arccosh, arcsinh, arctanh, rint, \n\
floor, ceil, pi, sqrt.\n\nUse ** for power operations (ex : x**2 for x squared).\n\
Note that you have to write 2*x for 'two times x' and 2x will not work.\nExample : V(x) = 2*cos(4*x**2)\n")
    
    BC_label = ttk.Label(root2, text = "Boundary\nconditions", image=img, compound="right")
    BC_label.image = img
    BC_label.grid(column=0, row=3, sticky=tk.E, ipadx=5, ipady=5)
    BC_variable = tk.StringVar(value = 0)
    BCD_rb = ttk.Radiobutton(root2, text = "Dirichlet", value = 0, variable=BC_variable)
    BCD_rb.grid(column=1, row=3, padx = 70)
    BCP_rb = ttk.Radiobutton(root2, text = "Periodic", value = 1, variable=BC_variable)
    BCP_rb.grid(column=1, row=3, sticky=tk.W)

    Tip_BC_label = Hovertip(BC_label, "\nIf 2L is the lenght of the box :\n\n\
    Periodic: The wave function ψ at position -L and time t equals ψ \n\
                at position L and time t for all t.\n\n\
    Dirichlet: The wave function ψ at position -L and time t, as well as ψ\n\
                at position L and time t, are both set to zero for all t.\n\nNote that Periodic conditions are not available iin 2D.")

    slider_N_value = tk.DoubleVar()
    slider_N_label = ttk.Label(root2, text = "Discretization\n      level", image=img, compound="right")
    slider_N_label.image = img
    slider_N_label.grid(column=0, row=4, sticky=tk.E, ipadx=5, ipady=5)
    slider_N_value_label = ttk.Label(root2, text = "200")
    slider_N_value_label.grid(column=1, row=4, sticky=tk.W, ipadx=5, ipady=5, padx = 20)
    
    slider_N = ttk.Scale(root2, from_=50, to=500, variable=slider_N_value, 
                          command = lambda event : slider_N_value_label.configure(text = "%.0f" % slider_N_value.get()))
    slider_N.grid(column=1, row=4, padx=5, pady=5)
    slider_N.set(200)

    Tip_slider_N = Hovertip(slider_N_label, "\nRefers to the number of points used to discretize the interval [-L, L], \n\
where 2L represents the lenght of the box (in this case, L = 10).\n\n\
Note that calculation time increases with the number of points\n(particularly for dynamic).\n")

    def finalisation():

        for letter in name_new_potential.get():
            if letter in '/:*?"<>|\\ ':
                showinfo("information", 'Potential name cannot contain the following symbols :\n / : * ? " < > | \\ and spaces.')
                return 0

        df = pa.read_csv('VP/potential_data_name.csv', index_col=0)

        for pot_save in [K for K in df.columns.values]:
            if pot_save == name_new_potential.get() :
                showinfo(title='Information', message='There is already a potential with that name.')
                return 0
            
        create_button['state'] = "disabled"
        cancel_button['state'] = "disabled"

        popup = tk.Toplevel()
        popup.title('QuantumSimLite')
        def on_closing():
            popup.destroy()
            attention['text'] = "Please, wait to avoid errors..."

        popup.protocol("WM_DELETE_WINDOW", on_closing)

        progress_bar = ttk.Progressbar(popup, orient='horizontal', mode='indeterminate', length=280)
        progress_bar.grid(column=0, row=0, columnspan=2, padx=10, pady=20)
        progress_bar.start()
        wait_label = ttk.Label(popup, text="Please wait ...")
        wait_label.grid(column=0, row=1, columnspan=2, padx=10, pady=5)

        back = OneD_VP.main(new_potential_form.get(), name_new_potential.get(), BC_variable.get(), \
                         int(slider_N_value.get()), int(dimension_variable.get()))
        
        if back == 0 :
            create_button['state'] = "!disabled"
            cancel_button['state'] = "!disabled"
            popup.destroy()
        
            showinfo("warning", message = "An error has occurred in your potential formula." + \
          " Please refer to the provided information on correctly formatting a potential formula located next to the entry.")

            return 0
        
        try :
            progress_bar.stop()        
        except :
            d = int(dimension_variable.get())
            os.remove("VP/" + str(d) + "D/" + name_new_potential.get() + ".bz2")
            os.remove("VP/" + str(d) + "D/show/" + name_new_potential.get() + ".bz2")
            os.remove("VP/" + str(d) + "D/" + name_new_potential.get() + "_vp.csv")
            
            try :
                create_button['state'] = "!disabled"
            except :
                return 0
            
            cancel_button['state'] = "!disabled"
            attention['text'] = ""

            return 0
        
        df[name_new_potential.get()] = [new_potential_form.get(), dimension_variable.get()]
        df.to_csv('VP/potential_data_name.csv')
        os.makedirs("media/" + name_new_potential.get())
        os.makedirs("media/" + name_new_potential.get() + "/dynamic")
        os.makedirs("media/" + name_new_potential.get() + "/Eigenstates")
        os.makedirs("VP/" + dimension_variable.get() + "D/dynamic/" + name_new_potential.get())

        create_button['state'] = "!disabled"
        cancel_button['state'] = "!disabled"
    
        popup.destroy()
        root2.destroy()
        root2.quit()


    def start_finalisation():
        threading.Thread(target=finalisation).start()

    create_button = ttk.Button(root2, text = "Create", image=img, compound="right", command = lambda : start_finalisation())
    create_button.grid(row = 5, column=1, padx=5, sticky=tk.W)

    Tip_create = Hovertip(create_button, "Initiate diagonalization. Please, don't close windows during calculation to avoid errors.\n\
If you have to, please close the main window.")

    cancel_button = ttk.Button(root2, text = "Quit", command = lambda: def_cancel_button())
    cancel_button.grid(row = 5, column=0)

    def def_cancel_button():
        root2.destroy()
        root2.quit()

    attention = ttk.Label(root2, text = "")
    attention.grid(column=1, row=6, ipadx=5, ipady=5)

    root2.protocol("WM_DELETE_WINDOW", def_cancel_button)
    root2.mainloop()

#######################################################
#######   QuantumSimLite made by @maximev131   ########
####    Reach me out at maxime.vinteler@yahoo.fr   ####
#######################################################
    
