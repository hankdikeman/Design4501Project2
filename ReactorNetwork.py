class Network:
    def __init__(self):
        self.component_set = {}

    def add_component(self, name, component):
        if name in self.component_set.keys():
            raise ValueError("this component name is already taken")
        else:
            self.component_set[name] = component

    def check_solution(self):
        pass
