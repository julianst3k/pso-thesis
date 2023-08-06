import pytest
import optimization.Archive_Tree as at
import numpy as np


class TestTree:

    def test_one(self):
        array_test = np.zeros((3,2))
        array_test[0,:] = [2,3]
        array_test[1,:] = [3,2]
        array_test[2,:] = [1,4]
        archive_tree = at.ArchiveTree()
        for val in array_test:
            archive_tree.add_value(val, 0, None, None)
        assert archive_tree.root.get_value()[0] == 2

    def test_two(self):
        array_test = np.zeros((3, 2))
        array_test[0, :] = [2, 3]
        array_test[1, :] = [3, 2]
        array_test[2, :] = [2, 2.5]
        archive_tree = at.ArchiveTree()
        for val in array_test:
            archive_tree.add_value(val, 0, None, None)
        assert archive_tree.root.get_right() is None

    def test_three(self):
        array_test = np.zeros((5, 2))
        array_test[0, :] = [5, 5]
        array_test[1, :] = [4, 6]
        array_test[2, :] = [6, 4]
        array_test[3, :] = [7, 3]
        array_test[4, :] = [4, 2]
        archive_tree = at.ArchiveTree()
        for val in array_test:
            archive_tree.add_value(val, 0, None, None)
        print(archive_tree.root.get_value())
        assert archive_tree.root.get_right() is None and archive_tree.root.get_left() is None

    def test_four(self):
        array_test = np.zeros((6, 2))
        array_test[0, :] = [5, 5]
        array_test[1, :] = [4, 6]
        array_test[2, :] = [6, 4]
        array_test[3, :] = [7, 3]
        array_test[4, :] = [3, 7]
        array_test[5, :] = [3.5, 3.5]
        archive_tree = at.ArchiveTree()
        for val in array_test:
            archive_tree.add_value(val, 0, None, None)
        archive_tree.print()
        assert archive_tree.root.get_right().get_value()[0] == 7

    def test_five(self):
        array_test = np.random.rand(100, 2)
        archive_tree = at.ArchiveTree()
        for val in array_test:
            archive_tree.add_value(val, 0, None, None)
        test_list = archive_tree.return_nodes()
        non_dominated = True
        node_values = [node.get_value() for node in test_list]
        non_dominated_index = []
        print(node_values)
        for value_of_node in node_values:
            for j,val in enumerate(array_test):
                if all([val[i]==value_of_node[i] for i,_ in enumerate(val)]):
                    non_dominated_index.append(j)
        for value_of_node in node_values:
            for j,val in enumerate(array_test):
                if j in non_dominated_index:
                    continue
                else:
                    non_dominated = non_dominated and not all([val[i]<=value_of_node[i] for i, _ in enumerate(val)])

                    if all([val[i]<=value_of_node[i] for i, _ in enumerate(val)]):
                        print(val, value_of_node)
        assert non_dominated
        assert len(test_list) == archive_tree.size