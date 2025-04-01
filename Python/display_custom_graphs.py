from display_functions import panel_display, restricted_nr
import numpy as np


example = np.array([[3, -1, -1, 0, -1], [-1, 2, -1, 0, 0], [-1, -1, 4, -1, -1],
                        [0, 0, 0, 0, 0], [0, 0, 0, -1, 1]])

name = "Thomas"
title = "Thomas's Restricted Numerical Range"
panel_display(example, name)
restricted_nr(example, title, name)

