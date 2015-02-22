#Bauer, Matt
#1000-631-613
#2015-02-8
#Assignment_03

from tkinter import *
from tkinter import simpledialog
from tkinter import filedialog
import time

class cl_widgets:
    def __init__(self, ob_root_window, ob_world=[]):
        self.ob_root_window = ob_root_window
        self.ob_world = ob_world
        #self.toolbar = cl_toolbar(self)
        self.loader_panel = LoaderPanel(self)
        self.rotator_panel = RotatorPanel(self)
        self.scale_panel = ScalePanel(self)
        self.translate_panel = TranslationPanel(self)
        self.ob_canvas_frame = cl_canvas_frame(self)
        #self.status = cl_statusBar_frame(self)
        self.ob_world.add_canvas(self.ob_canvas_frame.canvas)


class cl_canvas_frame:
    def __init__(self, master):
        # mouse variables
        self.left_mouse_click_hold = False
        self.right_mouse_click_hold = False

        self.master = master
        self.canvas = Canvas(master.ob_root_window, width=640, height=640, bg="yellow")
        self.canvas.pack(expand=YES, fill=BOTH)

        self.canvas.bind('<Configure>', self.canvas_resized_callback)
        self.canvas.bind("<ButtonPress-1>", self.left_mouse_click_callback)
        self.canvas.bind("<ButtonRelease-1>", self.left_mouse_release_callback)
        self.canvas.bind("<ButtonPress-3>", self.right_mouse_click_callback)
        self.canvas.bind("<ButtonRelease-3>", self.right_mouse_release_callback)
        self.canvas.bind("<Motion>", self.mouse_motion_callback)
        #self.canvas.bind("<Key>", self.key_pressed_callback)    
        self.canvas.bind("<Up>", self.up_arrow_pressed_callback)
        self.canvas.bind("<Down>", self.down_arrow_pressed_callback)
        self.canvas.bind("<Right>", self.right_arrow_pressed_callback)
        self.canvas.bind("<Left>", self.left_arrow_pressed_callback)
        self.canvas.bind("<Shift-Up>", self.shift_up_arrow_pressed_callback)
        self.canvas.bind("<Shift-Down>", self.shift_down_arrow_pressed_callback)
        self.canvas.bind("<Shift-Right>", self.shift_right_arrow_pressed_callback)
        self.canvas.bind("<Shift-Left>", self.shift_left_arrow_pressed_callback)
        self.canvas.bind("f", self.f_key_pressed_callback)
        self.canvas.bind("b", self.b_key_pressed_callback)

    def key_pressed_callback(self, event):
        print('key pressed')

    def up_arrow_pressed_callback(self, event):
        print('pressed up')

    def down_arrow_pressed_callback(self, event):
        print('pressed down')

    def right_arrow_pressed_callback(self, event):
        print('pressed right')

    def left_arrow_pressed_callback(self, event):
        print('pressed left')

    def shift_up_arrow_pressed_callback(self, event):
        self.canvas.world.translate(0, .1, 0, 1)

    def shift_down_arrow_pressed_callback(self, event):
        pass

    def shift_right_arrow_pressed_callback(self, event):
        pass

    def shift_left_arrow_pressed_callback(self, event):
        pass

    def f_key_pressed_callback(self, event):
        pass

    def b_key_pressed_callback(self, event):
        pass

    def left_mouse_click_callback(self, event):
        self.left_mouse_click_hold = True
        self.right_mouse_click_hold = False

    def left_mouse_release_callback(self, event):
        self.left_mouse_click_hold = False

    def right_mouse_click_callback(self, event):
        self.right_mouse_click_hold = True
        self.left_mouse_click_hold = False

    def right_mouse_release_callback(self, event):
        self.right_mouse_click_hold = False

    def mouse_motion_callback(self, event):
        if self.right_mouse_click_hold:
            self.master.ob_world.scale_with_mouse(event)
            self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas, event)
        if self.left_mouse_click_hold:
            self.master.ob_world.rotate_with_mouse(event)
            self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas, event)

    def canvas_resized_callback(self, event):
        self.canvas.config(width=event.width - 4, height=event.height - 4)

        self.canvas.pack()
        self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas, event)


class TranslationPanel:
    def __init__(self, master):
        self.master = master
        frame = Frame(master.ob_root_window)
        frame.pack()

        self.point = StringVar()
        self.point.set("[10,10,10]")

        # start widgets
        Label(frame, text="Translation ([dx,dy,dz])").pack(side=LEFT)

        self.point_text_field = Entry(frame, width=10, textvariable=self.point)
        self.point_text_field.pack(side=LEFT)

        Label(frame, text="Steps").pack(side=LEFT)
        self.degree_spin_box = Spinbox(frame, width=3, from_=0, to=360)
        self.degree_spin_box.pack(side=LEFT)

        self.file_dialog_button = Button(frame, text="Translate", fg="blue", command=self.rotate)
        self.file_dialog_button.pack(side=LEFT)

    def rotate(self):
        theta_segment = int(self.degree_spin_box.get()) / int(self.steps_spin_box.get())
        axis = int(self.radio_button_group.get())

        if axis != 4:
            for x in range(0, int(self.steps_spin_box.get())):
                self.master.ob_world.rotate_theta(axis, theta_segment)
                self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)
                self.master.ob_root_window.update()
                time.sleep(.05)

        else:
            point_a_vector = str(self.point_a.get())[1:-1].strip().split(',')
            point_b_vector = str(self.point_b.get())[1:-1].strip().split(',')
            for x in range(0, int(self.steps_spin_box.get())):
                self.master.ob_world.rotate_around_a_line(point_a_vector.copy(), point_b_vector.copy(), theta_segment)
                self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)
                self.master.ob_root_window.update()
                time.sleep(.05)


class RotatorPanel:
    def __init__(self, master):
        self.master = master
        frame = Frame(master.ob_root_window)
        frame.pack()

        self.point_a = StringVar()
        self.point_a.set("[0.0,0.0,0.0]")

        self.point_b = StringVar()
        self.point_b.set("[1.0,1.0,1.0]")

        # start widgets
        Label(frame, text="Rotation Axis:").pack(side=LEFT)

        self.radio_button_group = IntVar()
        Radiobutton(frame, text="X", variable=self.radio_button_group, value=1,
                    command=self.disable_text_fields).pack(side=LEFT)
        Radiobutton(frame, text="Y", variable=self.radio_button_group, value=2,
                    command=self.disable_text_fields).pack(side=LEFT)
        selected_radio = Radiobutton(frame, text="Z", variable=self.radio_button_group, value=3,
                                     command=self.disable_text_fields)
        selected_radio.pack(side=LEFT)
        Radiobutton(frame, text="Line AB", variable=self.radio_button_group, value=4,
                    command=self.activate_text_fields).pack(side=LEFT)
        selected_radio.select()

        Label(frame, text="A:").pack(side=LEFT)
        self.point_a_text_field = Entry(frame, width=10, textvariable=self.point_a)
        self.point_a_text_field.pack(side=LEFT)

        Label(frame, text="B:").pack(side=LEFT)
        self.point_b_text_field = Entry(frame, width=10, textvariable=self.point_b)
        self.point_b_text_field.pack(side=LEFT)
        self.disable_text_fields()

        Label(frame, text="Degree:").pack(side=LEFT)
        self.degree_spin_box = Spinbox(frame, width=3, from_=0, to=360)
        self.degree_spin_box.pack(side=LEFT)

        Label(frame, text="Steps:").pack(side=LEFT)
        self.steps_spin_box = Spinbox(frame, width=3, from_=1, to=10)
        self.steps_spin_box.pack(side=LEFT)

        self.file_dialog_button = Button(frame, text="Rotate", fg="blue", command=self.rotate)
        self.file_dialog_button.pack(side=LEFT)

    def rotate(self):
        theta_segment = int(self.degree_spin_box.get()) / int(self.steps_spin_box.get())
        axis = int(self.radio_button_group.get())

        if axis != 4:
            for x in range(0, int(self.steps_spin_box.get())):
                self.master.ob_world.rotate_theta(axis, theta_segment)
                self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)
                self.master.ob_root_window.update()
                time.sleep(.05)

        else:
            point_a_vector = str(self.point_a.get())[1:-1].strip().split(',')
            point_b_vector = str(self.point_b.get())[1:-1].strip().split(',')
            for x in range(0, int(self.steps_spin_box.get())):
                self.master.ob_world.rotate_around_a_line(point_a_vector.copy(), point_b_vector.copy(), theta_segment)
                self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)
                self.master.ob_root_window.update()
                time.sleep(.05)

    def activate_text_fields(self):
        self.point_a_text_field.configure(state="normal")
        self.point_b_text_field.configure(state="normal")

    def disable_text_fields(self):
        self.point_a_text_field.configure(state="disabled")
        self.point_b_text_field.configure(state="disabled")


class ScalePanel:
    def __init__(self, master):
        self.master = master
        frame = Frame(master.ob_root_window)
        frame.pack()

        self.point = StringVar()
        self.point.set("[0.0,0.0,0.0]")
        self.scale_factor_string = StringVar()
        self.scale_factor_string.set("[1,1,1]")

        # start widgets
        Label(frame, text="Scale Ratio:").pack(side=LEFT)

        self.radio_button_group = IntVar()
        selected_radio = Radiobutton(frame, text="All", variable=self.radio_button_group, value=1,
                                     command=self.activate_scale_factor)
        selected_radio.pack(side=LEFT)
        selected_radio.select()
        self.scale_factor = Spinbox(frame, width=4, from_=0.25, to=10, format="%.2f", increment=.25)
        self.scale_factor.pack(side=LEFT)
        Radiobutton(frame, text="[Sx,Sy,Sz]", variable=self.radio_button_group, value=2,
                    command=self.activate_entry).pack(side=LEFT)
        self.point_text_field = Entry(frame, width=10, textvariable=self.point)
        self.point_text_field.pack(side=LEFT)
        self.scale_factor_text_field = Entry(frame, width=10, textvariable=self.scale_factor_string)
        self.scale_factor_text_field.pack(side=LEFT)
        self.activate_scale_factor()

        Label(frame, text="Steps:").pack(side=LEFT)
        self.steps_spin_box = Spinbox(frame, width=3, from_=1, to=10)
        self.steps_spin_box.pack(side=LEFT)

        self.file_dialog_button = Button(frame, text="Scale", fg="blue", command=self.scale)
        self.file_dialog_button.pack(side=LEFT)

    def scale(self):
        if self.radio_button_group.get() == 1:
            factor = (float(self.scale_factor.get()) - 1.0) / float(self.steps_spin_box.get())
            saved_vertices = self.master.ob_world.vertices.copy()
            for x in range(0, int(self.steps_spin_box.get())):
                self.master.ob_world.vertices = saved_vertices.copy()
                self.master.ob_world.scale(factor * (x+1) + 1, factor * (x+1) + 1, factor * (x+1) + 1)
                self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)
                self.master.ob_root_window.update()
                time.sleep(.05)

        elif self.radio_button_group.get() == 2:

            split_factor = str(self.scale_factor_string.get())[1:-1].strip().split(',')

            factor_x = (float(split_factor[0]) - 1.0) / float(self.steps_spin_box.get())
            factor_y = (float(split_factor[1]) - 1.0) / float(self.steps_spin_box.get())
            factor_z = (float(split_factor[2]) - 1.0) / float(self.steps_spin_box.get())
            saved_vertices = self.master.ob_world.vertices.copy()
            for x in range(0, int(self.steps_spin_box.get())):
                self.master.ob_world.vertices = saved_vertices.copy()
                self.master.ob_world.scale(factor_x * (x+1) + 1, factor_y * (x+1) + 1, factor_z * (x+1) + 1)
                self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)
                self.master.ob_root_window.update()
                time.sleep(.05)

    def activate_scale_factor(self):
        self.scale_factor.configure(state="normal")
        self.scale_factor_text_field.configure(state="disabled")
        self.point_text_field.configure(state="disabled")

    def activate_entry(self):
        self.scale_factor.configure(state="disabled")
        self.scale_factor_text_field.configure(state="normal")


class LoaderPanel:
    def __init__(self, master):
        self.master = master
        frame = Frame(master.ob_root_window)
        frame.pack()

        self.var_filename = StringVar()
        self.var_filename.set('')

        self.label = Label(frame, text="Filename: ")
        self.label.pack(side=LEFT)

        self.text_field = Entry(frame, textvariable=self.var_filename)
        self.text_field.pack(side=LEFT)

        self.file_dialog_button = Button(frame, text="Browse", fg="blue", command=self.browse_file)
        self.file_dialog_button.pack(side=LEFT)

        self.load_button = Button(frame, text="Load", command=self.load_file)
        self.load_button.pack(side=LEFT)

    def browse_file(self):
        self.var_filename.set(filedialog.askopenfilename(filetypes=[("allfiles", "*"), ("pythonfiles", "*.txt")]))

    def load_file(self):
        self.master.ob_world.reset_lists()
        file = open(self.var_filename.get())
        for line in file.read().splitlines():
            leading_char = line[:1]
            variable_list = []

            for value in line[1:].strip().split():

                if leading_char == "v":
                    variable_list.append(float(value))
                elif leading_char == "f":
                    variable_list.append(int(value) - 1)
                elif leading_char == "r":
                    variable_list.append(float(value))
                elif leading_char == "n":
                    variable_list.append(float(value))
                elif leading_char == "u":
                    variable_list.append(float(value))
                elif leading_char == "p":
                    variable_list.append(float(value))
                elif leading_char == "w":
                    variable_list.append(float(value))
                elif leading_char == "s":
                    variable_list.append(float(value))

            if leading_char == "v":
                variable_list.append(float(1.0))
                self.master.ob_world.vertices.append(variable_list)
            elif leading_char == "f":
                self.master.ob_world.faces.append(variable_list)
            elif leading_char == "r":
                self.master.ob_world.r = variable_list
            elif leading_char == "n":
                self.master.ob_world.n = variable_list
            elif leading_char == "u":
                self.master.ob_world.u = variable_list
            elif leading_char == "p":
                self.master.ob_world.p = variable_list
            elif leading_char == "w":
                self.master.ob_world.window = variable_list
                self.master.ob_world.center_of_window = [(variable_list[0] + variable_list[1]) / 2.0,
                                                         (variable_list[2] + variable_list[3]) / 2.0,
                                                         0]
            elif leading_char == "s":
                self.master.ob_world.viewport = variable_list
        file.close()
        self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)




class cl_pannel_01:
    def __init__(self, master):
        self.master = master
        frame = Frame(master.ob_root_window)
        frame.pack()

        self.var_filename = StringVar()
        self.var_filename.set('')
        self.button = Button(frame, text="Hello", fg="red", command=self.say_hi)
        self.button.pack(side=LEFT)

        self.hi_there = Button(frame, text="Ask for a string", command=self.ask_for_string)
        self.hi_there.pack(side=LEFT)

        self.hi_there = Button(frame, text="Ask for a float", command=self.ask_for_string)
        self.hi_there.pack(side=LEFT)
        self.file_dialog_button = Button(frame, text="Open File Dialog", fg="blue", command=self.browse_file)
        self.file_dialog_button.pack(side=LEFT)

    def say_hi(self):
        print("hi there, everyone!")

    def ask_for_string(self):
        s = simpledialog.askstring('My Dialog', 'Please enter a string')
        print(s)

    def ask_for_float(self):
        f = simpledialog.askfloat('My Dialog', 'Please enter a string')
        print(f)

    def browse_file(self):
        self.var_filename.set(filedialog.askopenfilename(filetypes=[("allfiles", "*"), ("pythonfiles", "*.txt")]))
        filename = self.var_filename.get()
        print(filename)


class cl_pannel_02:
    def __init__(self, master):
        self.master = master
        frame = Frame(master.ob_root_window)
        frame.pack()
        self.button = Button(frame, text="Open Dialog", fg="blue", command=self.open_dialog_callback)
        self.button.pack(side=LEFT)

        self.hi_there = Button(frame, text="button 2", command=self.button2_callback)
        self.hi_there.pack(side=LEFT)


    def open_dialog_callback(self):
        d = MyDialog(self.master.ob_root_window)
        print(d.result)
        print("mydialog_callback pressed!")

    def button2_callback(self):
        print("button2 pressed!")


class MyDialog(simpledialog.Dialog):
    def body(self, master):

        Label(master, text="Integer:").grid(row=0, sticky=W)
        Label(master, text="Float:").grid(row=1, column=0, sticky=W)
        Label(master, text="String:").grid(row=1, column=2, sticky=W)
        self.e1 = Entry(master)
        self.e1.insert(0, 0)
        self.e2 = Entry(master)
        self.e2.insert(0, 4.2)
        self.e3 = Entry(master)
        self.e3.insert(0, 'Default text')

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=1, column=3)

        self.cb = Checkbutton(master, text="Hardcopy")
        self.cb.grid(row=3, columnspan=2, sticky=W)


    def apply(self):
        try:
            first = int(self.e1.get())
            second = float(self.e2.get())
            third = self.e3.get()
            self.result = first, second, third
        except ValueError:
            tkMessageBox.showwarning(
                "Bad input",
                "Illegal values, please try again"
            )


        #class StatusBar:

        #def __init__(self, master):
        #self.master=master
        #self.label = Label(self, bd=1, relief=SUNKEN, anchor=W)
        #self.label.pack(fill=X)

        #def set(self, format, *args):
        #self.label.config(text=format % args)
        #self.label.update_idletasks()

        #def clear(self):
        #self.label.config(text="")
        #self.label.update_idletasks()       


class cl_statusBar_frame:
    def __init__(self, master):
        self.master = master
        status = StatusBar(master.ob_root_window)
        status.pack(side=BOTTOM, fill=X)
        status.set('%s', 'This is the status bar')


    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()

    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()


class cl_menu:
    def __init__(self, master):
        self.master = master
        self.menu = Menu(master.ob_root_window)
        master.ob_root_window.config(menu=self.menu)
        self.filemenu = Menu(self.menu)
        self.menu.add_cascade(label="File", menu=self.filemenu)
        self.filemenu.add_command(label="New", command=self.menu_callback)
        self.filemenu.add_command(label="Open...", command=self.menu_callback)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.menu_callback)
        self.dummymenu = Menu(self.menu)
        self.menu.add_cascade(label="Dummy", menu=self.dummymenu)
        self.dummymenu.add_command(label="Item1", command=self.menu_item1_callback)
        self.dummymenu.add_command(label="Item2", command=self.menu_item2_callback)

        self.helpmenu = Menu(self.menu)
        self.menu.add_cascade(label="Help", menu=self.helpmenu)
        self.helpmenu.add_command(label="About...", command=self.menu_help_callback)

    def menu_callback(self):
        print("called the menu callback!")

    def menu_help_callback(self):
        print("called the help menu callback!")

    def menu_item1_callback(self):
        print("called item1 callback!")

    def menu_item2_callback(self):
        print("called item2 callback!")


class cl_toolbar:
    def __init__(self, master):
        self.master = master
        self.toolbar = Frame(master.ob_root_window)
        self.button = Button(self.toolbar, text="Draw", width=16, command=self.toolbar_draw_callback)
        self.button.pack(side=LEFT, padx=2, pady=2)

        self.button = Button(self.toolbar, text="Toolbar Button 2", width=16, command=self.toolbar_callback)
        self.button.pack(side=RIGHT, padx=2, pady=2)

        self.toolbar.pack(side=TOP, fill=X)

    def toolbar_draw_callback(self):
        self.master.ob_world.create_graphic_objects(self.master.ob_canvas_frame.canvas)
        #temp_canvas=self.master.ob_canvas_frame.canvas
        #line1=temp_canvas.create_line(0,0,temp_canvas.cget("width"),temp_canvas.cget("height"))
        #line2=temp_canvas.create_line(temp_canvas.cget("width"),0,0,temp_canvas.cget("height"))
        #oval=temp_canvas.create_oval(int(0.25*int(temp_canvas.cget("width"))),
        #int(0.25*int(temp_canvas.cget("height"))),
        #int(0.75*int(temp_canvas.cget("width"))),
        #int(0.75*int(temp_canvas.cget("height"))))

        print("called the draw callback!")

    def toolbar_callback(self):
        print("called the toolbar callback!")


