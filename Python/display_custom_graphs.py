from display_functions import panel_display, restricted_nr
import numpy as np


example = np.array([[3, -1, -1, -1, 0], [0, 1, -1, 0, 0], [-1, -1, 4, -1, -1],
                        [0, 0, 0, 0, 0], [0, 0, 0, -1, 1]])

panel_display(example, "Michael")
restricted_nr(example, "Michael's Restricted Numerical Range", name = "Michael")

