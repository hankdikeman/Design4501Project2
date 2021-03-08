"""
This file implements the reactor network optimization algorithm for the OME3-5
production process using the ReactorNetwork object and component objects.
"""
import numpy as np
from ReactorNetwork import Network
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber, Feed, Removal, ProductRemoval, Mixer, Splitter, FlashTank, StreamGen
from ConversionFunctions import MethanolReactor, FormaldehydeReactor, OMEReactor

TEMP_MEOH = 573
PRESS_MEOH = 37500

if __name__ == "__main__":
    # generate reactor network
    net1 = Network()
    # feed to reactor
    i1outlets = {'i1': StreamGen(H2=60, CO2=20)}
    net1.add_component('I1', Feed(i1outlets))
    # mixer for recycle and feed
    mix1inlets = {
        'i1': StreamGen(H2=60, CO2=20),
        's11': StreamGen()
    }
    mix1outlets = {'m11': StreamGen()}
    net1.add_component('M1', Mixer(mix1inlets, mix1outlets))
    # methanol reactor
    r1inlets = {'m11': StreamGen()}
    r1outlets = {'r11': StreamGen()}
    net1.add_component('R1', Reactor(TEMP_MEOH, PRESS_MEOH,
                                     r1inlets, r1outlets, MethanolReactor))
    # splitter for recycle
    s1inlets = {'f11': StreamGen()}
    s1outlets = {
        's11': StreamGen(),
        's12': StreamGen()
    }
    net1.add_component('S1', Splitter(s1inlets, s1outlets, 0.995, 's11'))
    # purge outlet
    purge1inlet = {'s12': StreamGen()}
    net1.add_component('Purge1', Removal(purge1inlet))
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
