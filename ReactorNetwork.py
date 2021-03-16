"""
Network class to contain all reaction network components and iterate over
components to find the steady state solution of system
"""
import numpy as np
import pprint
import collections
from Components import Component, HeatExchanger, Compressor, Turbine, Reactor, DistillationColumn, Absorber, Feed, Removal, ProductRemoval, Mixer, Splitter, FlashTank
from ConversionFunctions import MethanolReactor, FormaldehydeReactor, OMEReactor


class Network:
    LEARNING_PARAM = 0.15

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
            next_inlets = {**(self.component_set[component].calc_outlets()), **next_inlets}
            current_inlets = {**(self.component_set[component].get_inlets()), **current_inlets}
        # change all negative values to zero
        for k in current_inlets.keys():
            current_inlets[k][current_inlets[k]<0] = 0
        for k in next_inlets.keys():
            next_inlets[k][next_inlets[k]<0] = 0
        # use intersection of dictionary keysets to avoid key error
        # feed and outlets will not be iterated in this set-up (might need to be changed)
        iterable_streams = current_inlets.keys() & next_inlets.keys()
        # loop through components
        for component in self.component_set.keys():
            # check if component is fixed
            if not self.component_set[component].is_fixed():
                # get inlet streams of component
                inlets = self.component_set[component].get_inlets()
                print('component:',component,'\ninlets:')
                pprint.pprint(inlets)
                for inflow in (inlets.keys() & iterable_streams):
                    inlets[inflow] = np.array(inlets[inflow] + self.LEARNING_PARAM * \
                        (next_inlets[inflow] - inlets[inflow]))
                    pprint.pprint(inlets[inflow])
        return next_inlets

    def equilibrate_network(self):
        # check stream coupling
        self.check_stream_coupling()
        # initialize count variable
        n_iter = 0
        self.to_csv()
        # while the network is not marked as equilibrated
        while not self.equilibrated and n_iter < 1700:
            print("Iteration ", n_iter)
            # call iteration function
            next_inlets = self.iterate_network()
            # check if all components are at steady state solution
            self.equilibrated = True
            for component in self.component_set:
                if not self.component_set[component].is_fixed():
                    print('Converged:',self.component_set[component].check_solution(next_inlets))
                    self.equilibrated &= self.component_set[component].check_solution(
                        next_inlets)
            print('\nSystem equilibrated:',self.equilibrated)
            print('OME3-5 Output:',np.sum(list(self.component_set['ProductOutlet'].get_inlets().values())[0][10:13]),'moles')
            print('\n'+'='*100+'\n')
            n_iter += 1


    def to_csv(self):
        cols = np.array(['Species: ', 'H2', 'CO2', 'CO', 'H2O', 'MEOH', 'FA', 'N2',
                  'O2', 'OME1', 'OME2', 'OME3-', 'OME4', 'OME5', 'OME6'])
        with open('/Users/patrickgibbons/Desktop/git/Design4501Project2/Report2Data/'+'MassTable'+'.csv', 'a+b') as f:
            np.savetxt(f, cols, fmt = '%s', newline=',', delimiter = ',')
            f.write(b"\n")
        # loop through components and store inlets in dictionary
        inlets = {}
        for component in self.component_set.keys():
            # merge current inlets dict
            inlets = {**(self.component_set[component].get_inlets()), **inlets}

        for k, v in inlets.items():
            line = np.array([str(k)])
            line = np.append(line, v, axis=0)
            with open('/Users/patrickgibbons/Desktop/git/Design4501Project2/Report2Data/'+'MassTable'+'.csv', 'a+b') as f:
                np.savetxt(f, line, fmt = '%s', newline=',', delimiter = ',')
                f.write(b"\n")







    def check_stream_coupling(self):
        # check that every inlet has outlet and vice-versa
        inlet_names = []
        outlet_names = []
        for component in self.component_set.values():
            for outlets in component.outlets.keys():
                outlet_names = outlet_names + [outlets]
            for inlets in component.inlets.keys():
                inlet_names = inlet_names + [inlets]
        # sort lists of streams
        inlet_names.sort()
        outlet_names.sort()
        # print stream names
        print('-'*100 + '\n')
        print('inlet names:', inlet_names)
        print('outlet names:', outlet_names, '\n')
        print('-'*100 + '\n')
        # print unresolved streams
        print('UNRESOLVED STREAMS:')
        print(set(inlet_names) ^ set(outlet_names))
        if (set(inlet_names) ^ set(outlet_names)) == set():
            print('\tAll Streams Resolved Correctly! Reactor Network Ready to Go')
        # raise ValueError if some streams are unmatched
        if collections.Counter(inlet_names) != collections.Counter(outlet_names):
            raise ValueError('Not all streams have a start and end!!')


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
