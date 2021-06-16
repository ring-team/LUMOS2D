"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """

import numpy as np
from lib_data_structure.import_export_data import *


def normal_vector(point1, point2):
    """
    compute normalized normal vector of a segment
    @param point1: first segment endpoint
    @param point2: first segment endpoint
    @return: normalized normal vector
    """
    vect_n = np.array([-(point2.y_coord_ - point1.y_coord_), (point2.x_coord_ - point1.x_coord_)])
    norm = np.linalg.norm(vect_n)
    vect_n /= norm
    return vect_n


class Line:
    """
       A class to represent a polygonal line.

       @Attributes
       line_file_path_ : File path containing the description of the line (.msh or .mesh). The points have to be ordered
       l_vertices_ : liste of points (Point)
       l_normal_verctors_ : list of normal vectors of each line segment

       @Methods
       import_line(): define local edge index in the mesh
       compute_normal_vectors(): compute normalized normal vectors

    """
    def __init__(self, line_file_path):
        self.line_file_path_ = line_file_path
        self.l_vertices_ = []
        self.import_line()
        self.l_normal_verctors_ = []
        self.compute_normal_vectors()

    def import_line(self):
        """
        Fill the l_vertices_ by index order after reading line_file_path_
        """
        open_gmsh_file(self.line_file_path_)
        self.l_vertices_ = import_nodes()
        close_gmsh_file()

        self.l_vertices_.sort(key=lambda pt: pt.point_idx_)

    def compute_normal_vectors(self):
        """
        Compute normalized normal vector and fill the list of normal vectors of each line segment (l_normal_verctors_)
        """
        self.l_normal_verctors_.append(normal_vector(self.l_vertices_[0], self.l_vertices_[1]))

        for vertex in range(1, len(self.l_vertices_) - 1):
            self.l_normal_verctors_.append(normal_vector(self.l_vertices_[vertex - 1], self.l_vertices_[vertex + 1]))

        self.l_normal_verctors_.append(normal_vector(self.l_vertices_[-2], self.l_vertices_[-1]))


if __name__ == '__main__':
    my_line = Line("D:/Programmation/ringlab/LUMOS/python/data/frac1.msh")
    print(my_line.l_vertices_)
    print(my_line.l_normal_verctors_)
