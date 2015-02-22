#Bauer, Matt
#1000-631-613
#2015-02-8
#Assignment_03

from tkinter import messagebox
from Bauer_widgets_01 import *
from Bauer_graphics_01 import *


def close_window_callback(root):
    if messagebox.askokcancel("Quit", "Do you really wish to quit?"):
        root.destroy()


ob_root_window = Tk()
ob_root_window.protocol("WM_DELETE_WINDOW", lambda root_window=ob_root_window: close_window_callback(root_window))
ob_world = ClWorld()
cl_widgets(ob_root_window, ob_world)
ob_root_window.mainloop()    
    