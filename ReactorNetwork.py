"""
Network class to contain all reaction network components and iterate over
components to find the steady state solution of system
"""
import numpy as np
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber
from ConversionFunctions import MethanolReactor, FormaldehydeReactor, OMEReactor


class Network:
    LEARNING_PARAM = 0.05

    def __init__(self):
        self.component_set = {}
        self.equilibrated = False

    def add_component(self, name, component):
        # check if the components name is taken, and if so then raise an error
        if not isinstance(component, Component):
            raise ValueError('this is not a component type object')
        elif name in self.component_set.keys():
            raise ValueError('this component name is already taken')
        else:
            self.component_set[name] = component

    def iterate_network(self):
        # make dictionary to hold current iteration by stream name and vector
        current_inlets = {}
        next_inlets = {}
        # loop through components and store inlets in dictionary
        for component in self.component_set.keys():
            # merge current inlets dict
            current_inlets = {
                **self.component_set[component].get_inlets(), **current_inlets}
            next_inlets = {
                **self.component_set[component].calc_outlets(), **next_inlets}
        # use intersection of dictionary keysets to avoid key error
        # feed and outlets will not be iterated in this set-up (might need to be changed)
        iterable_streams = current_inlets.keys() & next_inlets.keys()
        # loop through inlet streams and iterate values
        for component in self.component_set.keys():
            for inlets in component.get_inlets():
                for inflow in (inlets.keys() & iterable_streams):
                    inlets[inflow] = inlets[inflow] - self.LEARNING_PARAM * \
                        (next_inlets[inflow] - inlets[inflow])

    def equilibrate_network(self):
        n_iter = 0
        # while the network is not marked as equilibrated
        while not self.equilibrated:
            print("Iteration ", n_iter)
            # call iteration function
            self.iterate_network()
            # check if all components are at steady state solution
            self.equilibrated = True
            for component in self.component_set:
                self.equilibrated &= component.check_solution()
            print(self.equilibrated)
            n_iter += 1


# unit tests for reaction network class
if __name__ == "__main__":
    print(__doc__)
    # generate network class
    net = Network()
    # generate several components of different classes, add each to network
    net.equilibrate_network()
