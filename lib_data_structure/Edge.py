"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
import numpy as np


class Edge:
    """
       A class to represent an edge.

       @Attributes
       point1_ : first end point (Point)
       point2_ : second end point (Point)
       edge_idx_ : edge index
       edge_tag_ : edge tag (line tag in the mesh)

       @Methods
       define_local_edge_idx(val): define local edge index in the mesh
       compute_normal_vector(): compute normalized normal vector
       compute_length(): compute

    """
    def __init__(self, point1, point2, edge_idx, edge_tag):
        self.point1_ = point1
        self.point2_ = point2
        self.edge_idx_ = edge_idx
        self.edge_tag_=edge_tag

    def define_local_edge_idx(self, val):
        self.local_edge_idx_ = val

    def compute_normal_vector(self):
        """
        Compute normalized normal vector
        @return: normalized normal vector
        """
        vect = np.array(
            [-(self.point2_.y_coord_ - self.point1_.y_coord_), (self.point2_.x_coord_ - self.point1_.x_coord_)])
        norm_vect = vect / np.linalg.norm(vect)
        return norm_vect

    def get_point_idx(self):
        return self.point1_.point_idx_, self.point2_.point_idx_

    def compute_length(self):
        """
        Compute edge length
        @return: edge length
        """
        vect = np.array(
            [(self.point2_.y_coord_ - self.point1_.y_coord_), (self.point2_.x_coord_ - self.point1_.x_coord_)])
        return np.linalg.norm(vect)