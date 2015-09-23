__author__ = 'xlinfr'

import Suews_wrapper_v10
import os
import Tkinter
import FileDialog

# working_path = os.path.dirname(__file__)
working_path = os.getcwd()
Suews_wrapper_v10.wrapper(working_path)