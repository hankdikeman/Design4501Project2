"""
Network class to contain all reaction network components and iterate over
components to find the steady state solution of system
"""
import numpy as np
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber, Feed, Removal, ProductRemoval, Mixer, Splitter, FlashTank
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
            next_inlets = {
                **self.component_set[component].calc_outlets(), **next_inlets}
            current_inlets = {
                **self.component_set[component].get_inlets(), **current_inlets}
        # use intersection of dictionary keysets to avoid key error
        # feed and outlets will not be iterated in this set-up (might need to be changed)
        iterable_streams = current_inlets.keys() & next_inlets.keys()
        # loop through inlet streams and iterate values
        for component in self.component_set.keys():
            if not component.isfixed():
                for inlets in component.get_inlets():
                    for inflow in (inlets.keys() & iterable_streams):
                        inlets[inflow] = inlets[inflow] - self.LEARNING_PARAM * \
                            (next_inlets[inflow] - inlets[inflow])
        return next_inlets

    def equilibrate_network(self):
        n_iter = 0
        # while the network is not marked as equilibrated
        while not self.equilibrated:
            print("Iteration ", n_iter)
            # call iteration function
            next_inlets = self.iterate_network()
            # check if all components are at steady state solution
            self.equilibrated = True
            for component in self.component_set:
                self.equilibrated &= component.check_solution(next_inlets)
            print(self.equilibrated)
            n_iter += 1


# unit tests for reaction network class
if __name__ == "__main__":
    print(__doc__)
    # generate network class
    net = Network()
    dict1 = {'mu1': np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])}
    dict2 = {'mu2': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    dict3 = {'mu3': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    dict4 = {'mu4': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    dict5 = {'mu5': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    dict6 = {'mu6': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    absorber1 = Absorber(300, 2, 0.97, 5, 3, dict1, dict2, dict3, dict4)
    recov = {'LK': (0.99, 1), 'HK': (0.001, 2)}
    distillationColumn = DistillationColumn(300, 2, recov, dict4, dict5, dict6)
    net.add_component('D1', distillationColumn)
    net.add_component('A1', absorber1)
    # generate several components of different classes, add each to network
    net.equilibrate_network()
