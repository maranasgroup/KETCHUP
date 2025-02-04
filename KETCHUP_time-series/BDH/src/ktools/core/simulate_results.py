"""
FORmulate Kinetic model

Formulates the kinetic models for the solver
"""

from pyomo.environ import *
from pyomo.dae import *
import numpy as np
import pandas as pd

from pyomo.dae.simulator import Simulator

