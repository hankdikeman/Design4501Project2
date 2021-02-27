"""
Implements an inheritance-based component substructure for reactor network,
couples components to inlets and outlets to other components. Determines outlet
compositions as a function of inlet compositions
"""
import numpy as np


# baseclass for component
class Component:
    def __init__(self, inlets, outlets):
        self.inlets = {}
        self.outlets = {}
        # add all inlet streams to inlet stream dict
        for comp, inflow in inlets:
            self.inlets[comp] = inflow
        # add all outlet streams to outlet stream dict
        for comp, outflow in outlets:
            self.outlets[comp] = outflow

    def change_stream(self):
        pass


class HeatExchanger(Component):
    def __init__(self, inlettemp, outlettemp, inlets, outlets):
        self.in_temp = inlettemp
        self.out_temp = outlettemp
        # run super constructor
        super(HeatExchanger, self).__init__(inlets, outlets)


class Compressor(Component):
    def __init__(self, in_t, out_t, in_p, out_p, inlets, outlets):
        self.in_temp = in_t
        self.out_temp = out_t
        self.in_press = in_p
        self.out_temp = out_p
        # run super constructor
        super(Compressor, self).__init__(inlets, outlets)


class Turbine(Component):
    def __init__(self, in_t, out_t, in_p, out_p, inlets, outlets):
        self.in_temp = in_t
        self.out_temp = out_t
        self.in_press = in_p
        self.out_temp = out_p
        # run super constructor
        super(Turbine, self).__init__(inlets, outlets)


# reactor baseclass
class Reactor(Component):
    def __init__(self, temperature, pressure, inlets, outlets):
        self.temperature = temperature
        self.pressure = pressure
        # run super constructor
        super(Reactor, self).__init__(inlets, outlets)
