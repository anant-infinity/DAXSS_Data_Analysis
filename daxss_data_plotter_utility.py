from tkinter import filedialog
from tkinter import *
from tkinter import ttk
from PIL import Image, ImageTk
from pandas import DataFrame
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import netCDF4 as nc

LARGE_FONT= ("Verdana", 12)
NORM_FONT= ("Verdana", 10)
SMALL_FONT= ("Verdana", 8)

def popupmsg(title, msg):
    popup = Tk()
    popup.wm_title(title)
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()

def printMetaDataX():
    child_window = Toplevel(root)
    child_window.geometry("900x450")
    child_window.title("X-axis Variable Meta Data")
    child_window.configure(background="#c4e0e5")
    scrollbar_y = Scrollbar(child_window, orient="vertical")
    scrollbar_y.pack(side=RIGHT, fill=Y)
    text = Text(child_window, height=50, width=110, font=("Times", 12, "bold"), yscrollcommand=scrollbar_y.set,
                   bg="#c4e0e5", bd=0)

    # create label
    Text_label = Label(child_window, text="Variable Meta Data", bg="#c4e0e5")
    Text_label.config(font=("Corier", 14))

    Text_label.pack(padx=10, pady=10)
    text.pack()

    var_array = []
    for var in daxsslevel1.variables.values():
        var_array.append(var)
    plot_variable_x = variable_x.get()
    for i in range(0, len(var_array)):
        if(var_array[i].name == plot_variable_x):
            print(var_array[i])
            text.insert(END, var_array[i])

    scrollbar_y.config(command=text.yview)

def printMetaDataY():
    child_window = Toplevel(root)
    child_window.geometry("900x450")
    child_window.title("Y-axis Variable Meta Data")
    child_window.configure(background="#c4e0e5")
    scrollbar_y = Scrollbar(child_window, orient="vertical")
    scrollbar_y.pack(side=RIGHT, fill=Y)
    text = Text(child_window, height=50, width=110, font=("Times", 12, "bold"), yscrollcommand=scrollbar_y.set,
                   bg="#c4e0e5", bd=0)

    # create label
    Text_label = Label(child_window, text="Variable Meta Data", bg="#c4e0e5")
    Text_label.config(font=("Corier", 14))

    Text_label.pack(padx=10, pady=10)
    text.pack()

    var_array = []
    for var in daxsslevel1.variables.values():
        var_array.append(var)
    plot_variable_x = variable_y.get()
    for i in range(0, len(var_array)):
        if(var_array[i].name == plot_variable_x):
            print(var_array[i])
            text.insert(END, var_array[i])

    scrollbar_y.config(command=text.yview)


def loadRawDataUser():
    global daxsslevel1
    global level_1_var_list
    #root = Tk()
    #root.withdraw()
    daxss_level1_file_path = filedialog.askopenfilename()
    daxsslevel1 = nc.Dataset(daxss_level1_file_path)
    level_1_var_list = list(daxsslevel1.variables)

    # Code for Scrolling
    main_frame = Frame(root)
    main_frame.pack(fill=BOTH, expand=1)

    main_canvas = Canvas(main_frame)
    main_canvas.pack(side=LEFT, fill=BOTH, expand=1)

    my_scrollbar = Scrollbar(main_frame, orient=VERTICAL, command=main_canvas.yview)
    my_scrollbar.pack(side=RIGHT, fill=Y)

    main_canvas.configure(yscrollcommand=my_scrollbar.set)
    main_canvas.bind('<Configure>', lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all")))

    global second_frame
    second_frame = Frame(main_canvas)
    main_canvas.create_window((0, 0), window=second_frame, anchor='nw')

    label_width = 20
    label_height = 10
    border_width = 2
    label_relief = "ridge"

    # Select x-axis plotting variable
    # X-axis Variable
    lx = Label(second_frame, fg='black', text="Select X Axis Variable", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    lx.grid(row=0, column=2,sticky='nesw')
    global variable_x
    variable_x = StringVar(root)
    variable_x.set("----")  # default value
    w_x = OptionMenu(second_frame, variable_x, *level_1_var_list)
    w_x.config(width=20)
    w_x.config(height=1)
    w_x.config(bg='lightblue')
    w_x.config(font=("Helvetica", 10))
    w_x.grid(row=0, column=3,sticky='nesw')

    # Show X-Meta Data Button
    button = Button(second_frame, width=20, height=1, text="Show X-Variable Metadata",
                    command=printMetaDataX, font=("Helvetica", 10),
                    bg="lightblue")
    button.grid(row=1, column=3, columnspan=1,sticky='nesw')


    global x_lim_low
    x_lim_low = StringVar()
    l4 = Label(second_frame, fg='black', text="X Limit Lower", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l4.grid(row=2, column=2)
    e4 = Entry(second_frame, fg='blue', textvariable=x_lim_low, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e4.grid(row=2, column=3)

    global x_lim_upper
    x_lim_upper = StringVar()
    l5 = Label(second_frame, fg='black', text="X Limit Upper", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l5.grid(row=3, column=2)
    e6 = Entry(second_frame, fg='blue', textvariable=x_lim_upper, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e6.grid(row=3, column=3)

    global x_scale
    x_scale = StringVar()
    l3 = Label(second_frame, fg='black', text="X Scale Type", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l3.grid(row=4, column=2)
    e3 = Entry(second_frame, fg='blue', textvariable=x_scale, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e3.grid(row=4, column=3)

    global x_lable
    x_lable = StringVar()
    l8 = Label(second_frame, fg='black', text="X Axis Label", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l8.grid(row=5, column=2)
    e9 = Entry(second_frame, fg='blue', textvariable=x_lable, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e9.grid(row=5, column=3)

    global x_dim_1
    x_dim_1 = StringVar()
    l81 = Label(second_frame, fg='black', text="X Dimension-1", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l81.grid(row=6, column=2)
    e91 = Entry(second_frame, fg='blue', textvariable=x_dim_1, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e91.grid(row=6, column=3)

    global x_dim_2
    x_dim_2 = StringVar()
    l811 = Label(second_frame, fg='black', text="X Dimension-2", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l811.grid(row=7, column=2)
    e911 = Entry(second_frame, fg='blue', textvariable=x_dim_2, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e911.grid(row=7, column=3)


    # Select y-axis plotting variable
    # Y-axis Variable
    ly = Label(second_frame, fg='black', text="Select Y Axis Variable", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    ly.grid(row=0, column=5,sticky='nesw')
    global variable_y
    variable_y = StringVar(root)
    variable_y.set("---")  # default value
    w_y = OptionMenu(second_frame, variable_y, *level_1_var_list)
    w_y.config(width=20)
    w_y.config(height=1)
    w_y.config(bg='lightblue')
    w_y.config(font=("Helvetica", 10))
    w_y.grid(row=0, column=6,sticky='nesw')

    # Show Y-Meta Data Button
    button = Button(second_frame, width=20, height=1, text="Show Y-Variable Metadata",
                    command=printMetaDataY, font=("Helvetica", 10),
                    bg="lightblue")
    button.grid(row=1, column=6, columnspan=1,sticky='nesw')

    #Plot Parameters Entry Labels
    global y_lim_low
    y_lim_low = StringVar()
    l1 = Label(second_frame, fg='black', text="Y Limit Lower", font=('Arial', 12, 'bold'), width=label_width,
                borderwidth=border_width, relief=label_relief, bg="lightblue")
    l1.grid(row=2, column=5)
    e1 = Entry(second_frame, fg='blue', textvariable=y_lim_low, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e1.grid(row=2, column=6)

    global y_lim_upper
    y_lim_upper = StringVar()
    l2 = Label(second_frame, fg='black', text="Y Limit Upper", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l2.grid(row=3, column=5)
    e2 = Entry(second_frame, fg='blue', textvariable=y_lim_upper, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e2.grid(row=3, column=6)

    global y_scale
    y_scale = StringVar()
    l3 = Label(second_frame, fg='black', text="Y Scale Type", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l3.grid(row=4, column=5)
    e3 = Entry(second_frame, fg='blue', textvariable=y_scale, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e3.grid(row=4, column=6)

    global y_lable
    y_lable = StringVar()
    l7 = Label(second_frame, fg='black', text="Y Axis Label", font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l7.grid(row=5, column=5)
    e8 = Entry(second_frame, fg='blue', textvariable=y_lable, font=('Arial', 12, 'bold'), width=label_width,
               borderwidth=border_width, relief=label_relief)
    e8.grid(row=5, column=6)

    global y_dim_1
    y_dim_1 = StringVar()
    l81 = Label(second_frame, fg='black', text="Y Dimension-1", font=('Arial', 12, 'bold'), width=label_width,
                borderwidth=border_width, relief=label_relief, bg="lightblue")
    l81.grid(row=6, column=5)
    e91 = Entry(second_frame, fg='blue', textvariable=y_dim_1, font=('Arial', 12, 'bold'), width=label_width,
                borderwidth=border_width, relief=label_relief)
    e91.grid(row=6, column=6)

    global y_dim_2
    y_dim_2 = StringVar()
    l811 = Label(second_frame, fg='black', text="Y Dimension-2", font=('Arial', 12, 'bold'), width=label_width,
                 borderwidth=border_width, relief=label_relief, bg="lightblue")
    l811.grid(row=7, column=5)
    e911 = Entry(second_frame, fg='blue', textvariable=y_dim_2, font=('Arial', 12, 'bold'), width=label_width,
                 borderwidth=border_width, relief=label_relief)
    e911.grid(row=7, column=6)

    global plot_type
    plot_type = StringVar()
    l7 = Label(second_frame, fg='black', text="Select Plot Type", font=('Arial', 12, 'bold'), width=15,
               borderwidth=border_width, relief=label_relief, bg="lightblue")
    l7.grid(row=5, column=0, sticky='nesw')
    e8 = Entry(second_frame, fg='blue', textvariable=plot_type, font=('Arial', 12, 'bold'), width=15,
               borderwidth=border_width, relief=label_relief)
    e8.grid(row=5, column=1, sticky='nesw')

    textlabel = Label(second_frame,
                      text="Instructions for Plotting:\n"
                           "1. Select plot type ('line' or 'scatter')\n"
                           "2. Select X-axis and Y-axis plotting variable \n"
                           "3. Select X and Y-axis lower and upper limits \n"
                           "4. Select X and Y-axis scale type ('linear','log',\n"
                           "'symlog','logit','function','functionlog' \n"
                           "5. Press Plot\n"
                           "Plot will be generated",
                      font=("Helvetica", 12), bg="lightblue",borderwidth=border_width, relief=label_relief, justify="left", fg="black")

    textlabel.grid(row=0, column=0, rowspan=5, columnspan=2)

    # Plot Button
    button = Button(second_frame, width=7, height=2, text="PLOT", command=generatePlot, font=("Helvetica", 12),
                    bg='#58F')
    button.grid(row=0, column=8, rowspan=6, sticky='nesw')

    mainloop()

def generatePlot():
    # Plotting Spectrum
    # X-axis plot parameters
    plot_variable_x = variable_x.get()

    if (x_dim_2.get() == 'NA'):
        if(x_dim_1.get() == ':'):
            x_plot = daxsslevel1[plot_variable_x][0, :]
        else:
            x_plot = daxsslevel1[plot_variable_x][0, int(x_dim_1.get())]
    else:
        if (x_dim_2.get() == ':'):
            x_plot = daxsslevel1[plot_variable_x][0, int(x_dim_1.get()), :]
        else:
            x_plot = daxsslevel1[plot_variable_x][0, int(x_dim_1.get()), int(x_dim_2.get())]

    plt.xlim([float(x_lim_low.get()), float(x_lim_upper.get())])
    plt.xscale(str(x_scale.get()))
    plt.xlabel(str(x_lable.get()))

    #Y axis plot parameters
    plot_variable_y = variable_y.get()
    if (y_dim_2.get() == 'NA'):
        if(y_dim_1.get() == ':'):
            y_plot = daxsslevel1[plot_variable_y][0, :]
        else:
            y_plot = daxsslevel1[plot_variable_y][0, int(y_dim_1.get())]
    else:
        if (y_dim_2.get() == ':'):
            y_plot = daxsslevel1[plot_variable_y][0, int(y_dim_1.get()), :]
        else:
            y_plot = daxsslevel1[plot_variable_y][0, int(y_dim_1.get()), int(y_dim_2.get())]


    plt.ylim([float(y_lim_low.get()), float(y_lim_upper.get())])
    plt.yscale(str(y_scale.get()))
    plt.ylabel(str(y_lable.get()))

    if(plot_type.get()=="line"):
        plt.plot(x_plot, y_plot, color='b')
    elif (plot_type.get() == "scatter"):
        plt.scatter(x_plot, y_plot, color='b')
    #plt.suptitle('DAXSS Plot')
    plt.legend()
    plt.show()

root = Tk()
root.title("DAXSS Data Plotter")
menu = Menu(root)
root.config(menu=menu)
menu.add_command(label="Select DAXSS Level-1 netCDF file", command=loadRawDataUser)
#filemenu.add_command(label="Download netCDF file", command=testPlot)

bg = ImageTk.PhotoImage(Image.open("images/background3.jpg").resize((1200,200), Image.ANTIALIAS))
daxss_icon = ImageTk.PhotoImage(Image.open("images/daxss_logo.PNG").resize((200,200), Image.ANTIALIAS))
inspire_icon = ImageTk.PhotoImage(Image.open("images/Inspire-logo.jpg").resize((200,200), Image.ANTIALIAS))


#create a label
window_label = Label(root, image = bg)
window_label.place(x=0, y=0)

#create canvas
my_canvas = Canvas(root, width = 1200, height = 200)
my_canvas.pack(fill = "both", expand = True)

#label to add DAXSS icon
isro_label = Label(root, image = daxss_icon, bd = 0, width = 200, height=200)
isro_label_window = my_canvas.create_window(0, 0, anchor="nw", window = isro_label)

#label to add INSPIRE icon
astre_label = Label(root, image = inspire_icon, bd = 0, width = 200, height=200)
astre_label_window = my_canvas.create_window(1000, 0, anchor="nw", window = astre_label)

#set image in canvas
my_canvas.create_image(0,0, image=bg, anchor = "nw")

#add a label project name
my_canvas.create_text(600, 100, text = "DAXSS Data Plotter Utility", font = ("Helvetica",30), fill = "#b6d0cf")

root.mainloop()