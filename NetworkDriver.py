"""
This file implements the reactor network optimization algorithm for the OME3-5
production process using the ReactorNetwork object and component objects.
"""
import numpy as np
from ReactorNetwork import Network
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber, Feed, Removal, ProductRemoval, Mixer, Splitter, FlashTank
from ConversionFunctions import MethanolReactor, FormaldehydeReactor, OMEReactor

if __name__ == "__main__":
    # generate reactor network
    net1 = Network()
    # feed to reactor

    # mixer for recycle and feed

    # methanol reactor

    # splitter for recycle

    # purge outlet

    # flash tank for reactor outlet

    # distillation column for methanol

    # water outlet

    # splitter for pure methanol

    # mixer for oxygen and methanol into formaldehyde

    # formaldehyde reactor

    # formaldehyde absorber

    # air removal

    # OME + methanol mixer

    # recycle + OME + methanol mixer

    # OME reactor

    # OME column

    # product outlet

    # water adsorber

    # water outlet
