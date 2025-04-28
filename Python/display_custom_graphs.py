from display_functions import panel_display, restricted_nr
import matplotlib
import numpy as np


example = np.array([[3, -1, -1, 0, -1], [-1, 2, -1, 0, 0], [-1, -1, 4, -1, -1],
                        [0, 0, 0, 0, 0], [0, 0, 0, -1, 1]])

name = ""
title = f"{name}'s Restricted Numerical Range"
color = None
panel_display(example, name, color)
restricted_nr(example, title, name, color)