"""
Network class to contain all reaction network components and iterate over
components to find the steady state solution of system
"""
import numpy as np


class Network:
    LEARNING_PARAM = 0.05

    def __init__(self):
        self.component_set = {}
        self.equilibrated = False

    def add_component(self, name, component):
        # check if the components name is taken, and if so then raise an error
        if name in self.component_set.keys():
            raise ValueError("this component name is already taken")
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
        # while the network is not marked as equilibrated
        while ~self.equilibrated:
            # call iteration function
            self.iterate_network()
            # check if all components are at steady state solution
            self.equilibrated = True
            for component in self.component_set:
                self.equilibrated &= component.check_solution()
