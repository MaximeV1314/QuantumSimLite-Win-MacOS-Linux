import tkinter as tk
from tkinter import ttk
import pandas as pa
import OneD_VP
from PIL import Image, ImageTk
import os
from idlelib.tooltip import Hovertip
from tkinter.messagebox import showinfo
import threading

def create_new_dynamic():

    """
    Frame that enables the creation of a new potential. Call OneD_VP.py function to calculate eigenvalues and eigenvectors
    """

    frame = tk.Toplevel(bg="white")
    frame.geometry("400x250")
    frame.resizable(False, False)
    frame.title("New Potential")

    frame.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=3)
    frame.columnconfigure(2, weight=1)
    frame.columnconfigure(3, weight=3)

    image = Image.open(os.path.dirname(__file__) + "/info_button.png", mode = 'r')
    resized_image= image.resize((15,15))
    img = ImageTk.PhotoImage(resized_image)

    name_new_dynamic = tk.StringVar()
    ttk.Label(frame, text = "Name :", background="white").grid(column=0, row=0, sticky=tk.E, padx = 5, pady = 5)
    name_entry = ttk.Entry(frame, textvariable=name_new_dynamic)
    name_entry.grid(column=1, row=0, padx = 5, pady = 5, sticky=tk.W)

    potential_list = [K for K in pa.read_csv('VP/potential_data_name.csv', index_col=0).columns.values]

    def dimensional(event):
        df_pot = pa.read_csv('VP/potential_data_name.csv', index_col=0)
        if int(df_pot[potential_choice_dyn.get()][1]) == 1 :
            slider_y0['state'] = 'disabled'
            slider_ky0['state'] = 'disabled'
        else :
            slider_y0['state'] = '!disabled'
            slider_ky0['state'] = '!disabled'

    potential_label = ttk.Label(frame, text = "Potential : ", image=img, compound="right", background="white")
    potential_label.image = img
    potential_label.grid(column=0, row=1, sticky=tk.E, ipadx=5, ipady=5)
    potential_choice_dyn = tk.StringVar()
    potential_cb_dyn = ttk.Combobox(frame, textvariable = potential_choice_dyn, values = potential_list)
    potential_cb_dyn.grid(column=1, row=1, sticky= tk.W, padx = 5, pady = 5)
    potential_cb_dyn.bind("<<ComboboxSelected>>", dimensional)
    potential_cb_dyn.current(0)
    potential_cb_dyn['state'] = 'readonly'

    Tip_potental = Hovertip(potential_label, "Select one of the potentials for the dynamic simulation.")


    iw_label = ttk.Label(frame, text = "Initial wave : ",  image=img, compound="right", background="white")
    iw_label.image = img
    iw_label.grid(column=0, row=2, sticky=tk.E, ipadx=5, ipady=5)

    Tip_iw_label = Hovertip(iw_label, "Initial wave in position representation.\n\nYou can choose between \
                            a gaussian, lorentzian or an uniform density")

    initial_wave_list = ["Gaussian", "Lorentzian", "Uniform density"]

    def def_iw(event):
        if initial_wave_choice.get() == "Uniform density" :
            slider_eps["state"] = "disabled"
        else :
            slider_eps["state"] = "!disabled"

    initial_wave_choice = tk.StringVar()
    initial_wave_cb = ttk.Combobox(frame, textvariable = initial_wave_choice, values = initial_wave_list)   #, postcommand = lambda : def_iw
    initial_wave_cb.grid(column=1, row=2, sticky= tk.W, padx = 5, pady = 5)
    initial_wave_cb.bind("<<ComboboxSelected>>", def_iw)
    initial_wave_cb.current(0)
    initial_wave_cb['state'] = 'readonly'


    slider_eps_value = tk.DoubleVar(value=2)
    slider_eps_label = ttk.Label(frame, text = "σ", image=img, compound="right", background="white")
    slider_eps_label.image = img
    slider_eps_label.grid(column=2, row=2, ipadx=5, ipady=5, padx = 20)
    slider_eps_value_label = ttk.Label(frame, text = "2.0", background="white")
    slider_eps_value_label.grid(column=2, row=2, sticky=tk.E, ipady=5, padx=0)

    Tip_slider_eps = Hovertip(slider_eps_label, "Width at half maximum of the Gaussian/Lorentzian.")
    
    slider_eps = ttk.Scale(frame, from_=0.3, to=4, variable=slider_eps_value, 
                          command = lambda event : slider_eps_value_label.configure(text = "%.1f" % slider_eps_value.get()))
    slider_eps.grid(column=3, row=2, padx=5, pady=5)

    
    slider_x0_value = tk.DoubleVar()
    slider_x0_label = ttk.Label(frame, text = "X0", image=img, compound="right", background="white")
    slider_x0_label.grid(column=0, row=3, ipadx=5, ipady=5)
    slider_x0_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_x0_value_label.grid(column=0, row=3, sticky=tk.E, ipady=5)
    
    slider_x0 = ttk.Scale(frame, from_=-10, to=10, variable=slider_x0_value, 
                          command = lambda event : slider_x0_value_label.configure(text = "%.1f" % slider_x0_value.get()))
    slider_x0.grid(column=1, row=3, padx=5, pady=5)

    Tip_y = Hovertip(slider_x0_label, "Initial position for horizontal axis.")

    slider_y0_value = tk.DoubleVar()
    slider_y0_label = ttk.Label(frame, text = "Y0", image=img, compound="right", background="white")
    slider_y0_label.grid(column=2, row=3, ipadx=5, ipady=5)
    slider_y0_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_y0_value_label.grid(column=2, row=3, sticky=tk.E, ipady=5)
    
    slider_y0 = ttk.Scale(frame, from_=-10, to=10, variable=slider_y0_value, 
                          command = lambda event : slider_y0_value_label.configure(text = "%.1f" % slider_y0_value.get()))
    slider_y0.grid(column=3, row=3, padx=5, pady=5)
    slider_y0['state'] = 'disabled'
    
    Tip_y = Hovertip(slider_y0_label, "Initial position for vertical axis. This parameter is not available in 1D.")

    slider_kx0_value = tk.DoubleVar()
    slider_kx0_label = ttk.Label(frame, text = "Px0", image=img, compound="right", background="white")
    slider_kx0_label.grid(column=0, row=4, ipadx=5, ipady=5)
    slider_kx0_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_kx0_value_label.grid(column=0, row=4, sticky=tk.E, ipady=5)
    
    slider_kx0 = ttk.Scale(frame, from_=-15, to=15, variable=slider_kx0_value, 
                          command = lambda event : slider_kx0_value_label.configure(text = "%.1f" % slider_kx0_value.get()))
    slider_kx0.grid(column=1, row=4, padx=5, pady=5)

    Tip_k = Hovertip(slider_kx0_label, "Initial momentum for horizontal axis.")

    slider_ky0_value = tk.DoubleVar()
    slider_ky0_label = ttk.Label(frame, text = "Py0", image=img, compound="right", background="white")
    slider_ky0_label.grid(column=2, row=4, ipadx=5, ipady=5)
    slider_ky0_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_ky0_value_label.grid(column=2, row=4, sticky=tk.E, ipady=5)
    
    slider_ky0 = ttk.Scale(frame, from_=-15, to=15, variable=slider_ky0_value, 
                          command = lambda event : slider_ky0_value_label.configure(text = "%.1f" % slider_ky0_value.get()))
    slider_ky0.grid(column=3, row=4, padx=5, pady=5)
    slider_ky0['state'] = 'disabled'

    Tip_ky = Hovertip(slider_ky0_label, "Initial momentum for vertical axis. This parameter is not available in 1D.")

    slider_time_value = tk.DoubleVar()
    slider_time_label = ttk.Label(frame, text = "Time", image = img, compound= "right", background="white")
    slider_time_label.image = img
    slider_time_label.grid(column=0, row=5, ipadx=5, ipady=5)
    slider_time_value_label = ttk.Label(frame, text = "0.0", background="white")
    slider_time_value_label.grid(column=0, row=5, sticky=tk.E, ipady=5)
    
    slider_time = ttk.Scale(frame, from_=10, to=60, variable=slider_time_value,
                          command = lambda event : slider_time_value_label.configure(text = str(slider_time_value.get())[:2]))
    slider_time.grid(column=1, row=5, padx=5, pady=5)
    slider_time.set(15)

    Tip_time = Hovertip(slider_time_label, "One unit of time = 25 frames and is approximately 1 second.")

    start_button = ttk.Button(frame, text = "Start", image = img, compound="right", command = lambda : start_finalisation())
    start_button.image = img
    start_button.grid(column=1, row=6, columnspan=2, sticky=tk.W, padx=10)

    def start_finalisation():
        threading.Thread(target=finalisation).start()

    def finalisation():

        for letter in name_new_dynamic.get():
            if letter in '/:*?"<>|\\ ' :
                showinfo("information", 'Potential name cannot contain the following symbols :\n / : * ? " < > | \\ and spaces')
                return 0
            
        df_pot = pa.read_csv('VP/potential_data_name.csv', index_col=0)
        d = int(df_pot[potential_cb_dyn.get()][1])
    
        dynamic_list_ = os.listdir("VP/" + str(d) + "D/dynamic/" + potential_cb_dyn.get())
        if dynamic_list_ == []:
            dynamic_list = [""]
        else :
            dynamic_list = []
            for a in dynamic_list_ :
                if a[-9:] == "_info.csv" :
                    dynamic_list.append(a[0:-9])

        for dyn_save in dynamic_list:
            if dyn_save == name_new_dynamic.get() :
                showinfo(title='Information', message='There is already a dynamic with that name.')
                return 0
            
        start_button['state'] = "disabled"
        cancel_button['state'] = "disabled"

        OneD_VP.dynamique(potential_choice_dyn.get(), initial_wave_choice.get(), \
           slider_x0_value.get(), slider_kx0_value.get(), slider_time_value.get(), 0.04, slider_eps_value.get(), \
               slider_y0_value.get(), slider_ky0_value.get(), name_new_dynamic.get())

        frame.destroy()
        frame.quit()

    Tip_start = Hovertip(start_button, "Initiate dynamic calculaton.")

    cancel_button = ttk.Button(frame, text = "Quit", command = lambda: def_cancel_button())
    cancel_button.grid(row = 6, column=0)


    def def_cancel_button():
        frame.destroy()
        frame.quit()

    frame.protocol("WM_DELETE_WINDOW", def_cancel_button)
    frame.mainloop()

#create_new_dynamic()

#######################################################
#######   QuantumSimLite made by @maximev131   ########
####    Reach me out at maxime.vinteler@yahoo.fr   ####
#######################################################
    
