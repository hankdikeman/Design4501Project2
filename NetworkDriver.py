"""
This file implements the reactor network optimization algorithm for the OME3-5
production process using the ReactorNetwork object and component objects.
"""
import numpy as np
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber, Feed, Removal, ProductRemoval, Mixer, Splitter, FlashTank
from ConversionFunctions import MethanolReactor, FormaldehydeReactor, OMEReactor

if __name__ == "__main__":
