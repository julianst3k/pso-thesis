import numpy as np


class ArchiveTreeController:

    def __init__(self, evaluated_values, particles, hypercubes=None, is_agg = False):
        self.values = evaluated_values
        self.particles = particles
        self.tree = self.create_tree(is_agg)
        self.hypercubes = hypercubes
        self.non_dominated_nodes = []

    def create_tree(self, is_aggs = False):
        tree = ArchiveTree()
        for i, val in enumerate(self.values):
            if is_aggs:
                tree.add_value(val, 0, i, self.particles[i].copy())
            else:
                tree.add_value(val, 0, i, self.particles[i, :].copy())
        return tree

    def update_tree(self, nodes):
        node_dictionary = {}
        for node in nodes:
            val = node.get_value()
            index = self.get_hypercubes().return_index(val)
            if index[0]>=10 or index[0]<0 or index[1]>=10 or index[1]<0:
                if 0 in node_dictionary:
                    node_dictionary[0].append(node)
                else:
                    node_dictionary[0] = [node]
            elif self.get_hypercubes().storage[index[0], index[1]] is None:
                if 0 in node_dictionary:
                    node_dictionary[0].append(node)
                else:
                    node_dictionary[0] = [node]
            else:
                if len(self.get_hypercubes().storage[index[0], index[1]]) in node_dictionary:
                    node_dictionary[len(self.get_hypercubes().storage[index[0], index[1]])].append(node)
                else:
                    node_dictionary[len(self.get_hypercubes().storage[index[0], index[1]])] = [node]
        keys = node_dictionary.keys()
        keys_sorted = sorted(keys)
        for key in keys_sorted:
            for node in node_dictionary[key]:
                self.tree.add_value_node(node)


    def return_non_dominant_nodes(self):
        if self.tree.root is None:
            return []
        return self.tree.return_nodes()

    def get_hypercubes(self):
        return self.hypercubes


class ArchiveTree:

    def __init__(self, top=None):
        self.root = None
        self.nodes = []
        self.dominated_nodes = []
        self.size = 0
        self.top = top

    def add_value_node(self, node):
        value = node.get_value()
        key = node.get_key()
        index = node.get_index()
        particle = node.get_particle()
        if value[0]<-1e-3 and value[1]<-1e-3:
            self.add_value(value, key, index, particle, debug=True)

    def add_value(self, value, key, index, particle, debug=False):
        if self.root is None:
            self.root = ArchiveTreeNode(value, key, index, particle)
            self.size = 1
        elif self.root.not_dominated(value):
            if self.top is not None and self.size >= self.top:
                val = self.root.conditional_add(value, key, index, particle)
            else:
                val = self.root.add(value, key, index, particle, debug=debug)
            self.size += val
        else:
            new_root, deleted = self.root.prune_dominated(value)
            self.size -= deleted
            if new_root is not None:
                self.root = new_root
                self.size += new_root.add(value, key, index, particle, debug=debug)
            else:
                self.root = ArchiveTreeNode(value, key, index, particle, debug=debug)
                self.size = 1
        if debug:
            print(value, particle)

    def print(self):
        self.root.print()

    def return_nodes(self):
        return self.root.return_nodes()


class ArchiveTreeNode:

    def __init__(self, value, key, index, particle, debug=False):
        self.key_value = value[key]
        self.key = key
        self.value = value
        self.particle = particle
        self.index = index
        self.left = None
        self.right = None

    def get_key(self):
        return self.key

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def get_index(self):
        return self.index

    def get_particle(self):
        return self.particle

    def get_value(self):
        return self.value

    def get_key_value(self):
        return self.key_value

    def set_left(self, node):
        self.left = node

    def set_right(self, node):
        self.right = node

    def is_dominated(self, other_node):
        return all([self.value[i] < val for i, val in enumerate(other_node.get_value())])

    def not_dominated(self, value):
        return any([self.value[i] < val for i, val in enumerate(value)])

    def add(self, value, key, index, particle, debug=False, root = False):
        if all([self.value[i] <= val for i, val in enumerate(value)]):
            return 0
        if self.not_dominated(value):
            if value[key] < self.get_key_value():
                if self.get_left() is None:
                    new_node = ArchiveTreeNode(value, key, index, particle, debug)
                    self.set_left(new_node)
                    return 1
                else:
                    return self.get_left().add(value, key, index, particle, debug)
            else:
                if self.get_left() is not None:
                    explore_bool = self.get_left().explore_dom(value)
                    if explore_bool:
                        return 0
                if self.get_right() is None:
                    new_node = ArchiveTreeNode(value, key, index, particle, debug)
                    self.set_right(new_node)
                    return 1
                else:
                    return self.get_right().add(value, key, index, particle, debug)
        else:
            new_node, nodes_deleted = self.prune_dominated(value)
            if new_node is not None:
                self.replace_attributes(new_node)
                return self.add(value, key, index, particle, debug)-nodes_deleted
            else:
                self.replace_attributes(ArchiveTreeNode(value, key, index, particle, debug))
                return 1-nodes_deleted

    def conditional_add(self, value, key, index, particle, root = False):
        if all([self.value[i] <= val for i, val in enumerate(value)]):
            return 0
        if self.not_dominated(value):
            if value[key] < self.get_key_value():
                if self.get_left() is None:
                    return 0
                else:
                    return self.get_left().conditional_add(value, key, index, particle)
            else:
                if self.get_left() is not None:
                    explore_bool = self.get_left().explore_dom(value)
                    if explore_bool:
                        return 0
                if self.get_right() is None:
                    return 0
                else:
                    return self.get_right().conditional_add(value, key, index, particle)
        else:
            new_node, nodes_deleted = self.prune_dominated(value)
            if new_node is not None:
                self.replace_attributes(new_node)
                return self.add(value, key, index, particle)-nodes_deleted
            else:
                self.replace_attributes(ArchiveTreeNode(value, key, index, particle))
                return 1-nodes_deleted




    def replace_attributes(self, node):
        self.set_left(node.get_left())
        self.set_right(node.get_right())
        self.set_particle(node.get_particle())
        self.set_key_value(node.get_key_value())
        self.set_index(node.get_index())
        self.set_value(node.get_value())

    def set_index(self, index):
        self.index = index

    def set_particle(self, particle):
        self.particle = particle

    def set_value(self, value):
        self.value = value

    def set_key_value(self, key_value):
        self.key_value = key_value

    def explore_dom(self, value):
        if all([self.value[i] <= val for i, val in enumerate(value)]):
            return True
        elif self.get_right() is not None:
            return self.get_right().explore_dom(value)
        return False

    def prune_dominated(self, value):
        deleted_left = deleted_right = 0
        if self.get_left() is not None:
            node, deleted_left = self.get_left().find_dominated(value)
            self.set_left(node)
        if self.get_right() is not None:
            node, deleted_right = self.get_right().find_dominated(value)
            self.set_right(node)
        if self.get_left() is not None:
            replace_node = self.get_left().replace_max()
            if replace_node == self.get_left():
                replace_node.set_left(replace_node.get_left())
                replace_node.set_right(self.get_right())
            else:
                replace_node.set_left(self.get_left())
                replace_node.set_right(self.get_right())
        elif self.get_right() is not None:
            replace_node = self.get_right().replace_min()
            if replace_node == self.get_right():
                replace_node.set_left(self.get_left())
                replace_node.set_right(replace_node.get_right())
            else:
                replace_node.set_left(self.get_left())
                replace_node.set_right(self.get_right())
        else:
            return None, deleted_left+deleted_right+1
        return replace_node, deleted_left+deleted_right+1

    def replace_min(self):

        if self.get_left() is None:
            return self
        elif self.get_left().get_left() is None:
            return_node = self.get_left()
            self.set_left(return_node.get_right())
            return return_node
        else:
            return self.get_left().replace_min()

    def replace_max(self):
        if self.get_right() is None:
            return self
        elif self.get_right().get_right() is None:
            return_node = self.get_right()
            self.set_right(return_node.get_left())
            return return_node

        else:
            return self.get_right().replace_max()

    def find_dominated(self, value):
        cantidad = 0
        if self.not_dominated(value):
            if value[1] < self.get_value()[1] and self.get_left() is not None:
                node, cantidad_new = self.get_left().find_dominated(value)
                cantidad += cantidad_new
                self.set_left(node)
            elif self.get_right() is not None:
                node, cantidad_new = self.get_right().find_dominated(value)
                cantidad += cantidad_new
                self.set_right(node)
            return self, cantidad
        else:
            return self.prune_dominated(value)

    def print(self):
        if self.get_left() is not None:
            self.get_left().print()
        print(self.get_value(), self.get_particle())
        if self.get_right() is not None:
            self.get_right().print()

    def return_nodes(self):
        node_list = []
        if self.get_left() is not None:
            node_list.extend(self.get_left().return_nodes())
        node_list.append(self)
        if self.get_right() is not None:
            node_list.extend(self.get_right().return_nodes())

        return node_list

class ArchiveTreeNodeCollector:

    def __init__(self):
        self.nodes = []

    def add_node(self, node):
        self.nodes.append(node)

    def get_nodes(self):
        return self.nodes

    def __len__(self):
        return len(self.nodes)

    def remove(self, node):
        self.nodes.remove(node)

    def get_random_node(self):
        return self.nodes[np.random.choice(len(self.nodes))]

