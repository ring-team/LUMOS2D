"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
import numpy as np


class Triangle:
    """
       A class to represent a triangle.

       @Attributes
       l_trgl_vertices_ : list of triangle vertices (Point)
       l_trgl_vertices_idx_ : list of triangle vertices indexes (int)
       l_trgl_edges_ : list of triangle edges (Edge)
       triangle_idx_ : triangle index
       triangle_reg_: triangle tag (=region index)
       l_local_edge_index_: list of triangle edges indexes (int)
       surf_: triangle surface

       @Methods
       find_local_point_index(): define local vertices indexes in the mesh
       find_local_edge_index(): define local edges indexes in the mesh
       find_common_edge(other): find common edge between two adjacent triangles
       orientation(): compute cross product in the triangle. The sign indicates triangle orientation
       orientation_new_vrtx(vertices_idx, new_point): compute cross product in the potential new triangle.
                                                    The sign indicates triangle orientation
       choose_direction(self, mesh_trgl_orient, new_point): Find the edge to continue the walk on triangles
       compute_surface(self, vertices_idx, point): compute triangle surface
       determine_barycentric_coeff(self, point): compute the barycentric coefficients for a point in the triangle
       compute_gradient_matrix(): compute the triangle gradient
       height(base): compute triangle height length (from triangle area and base edge)
       minimal_height(): compute triangle minimal height length


    """
    def __init__(self, trgl_vertices, trgl_edges, triangle_idx, triangle_reg):
        self.l_trgl_vertices_ = trgl_vertices
        self.l_trgl_vertices_idx_ = [vertex.point_idx_ for vertex in self.l_trgl_vertices_]
        self.l_trgl_edges_ = trgl_edges
        self.triangle_idx_ = triangle_idx
        self.triangle_reg_ = triangle_reg
        self.l_local_edge_index_ = []

        self.find_local_point_index()
        self.find_local_edge_index()

        self.surf_ = self.compute_surface((0, 1), self.l_trgl_vertices_[2])

    def get_l_trgl_vertices(self):
        return self.l_trgl_vertices_

    def get_vertices_coordinates(self):
        l_coordinates = []
        for vertex in self.l_trgl_vertices_:
            l_coordinates.append((vertex.x_coord_, vertex.y_coord_, vertex.z_coord_))

        return l_coordinates

    def find_local_point_index(self):
        """
        Define local index for each vertex
        """
        for vertex in self.l_trgl_vertices_:
            local_id = self.triangle_idx_ * 3 + self.l_trgl_vertices_.index(vertex)
            vertex.define_local_point_idx(local_id)

    def find_local_edge_index(self):
        """
        Define local index for each edge
        """

        for edge in self.l_trgl_edges_:
            local_id = self.triangle_idx_ * 3 + self.l_trgl_edges_.index(edge)
            edge.define_local_edge_idx(local_id)

            self.l_local_edge_index_.append(local_id)

    def find_common_edge(self, other):
        """
        Find common edge between two adjacent triangles
        :param other: other Triangle
        :return: end points indexes
        """
        l_common_vertices = [other.l_trgl_vertices_idx_[id] for id in range(3) if
                             other.l_trgl_vertices_idx_[id] in self.l_trgl_vertices_idx_]
        return l_common_vertices[0], l_common_vertices[1]

    def orientation(self):
        """
        det = (x1-x0)*(y2-y1)-(y1-y0)*(x2-x1)
        Compute cross product in the triangle. The sign indicates triangle orientation

        """
        l_coordinates = self.get_vertices_coordinates()

        det = (l_coordinates[1][0] - l_coordinates[0][0]) * (l_coordinates[2][1] - l_coordinates[1][1]) - (
                l_coordinates[1][1] - l_coordinates[0][1]) * (l_coordinates[2][0] - l_coordinates[1][0])

        return det

    def orientation_new_vrtx(self, vertices_idx, new_point):
        """
        det = (x1-x0)*(y2-y1)-(y1-y0)*(x2-x1)
        Compute cross product in the new potential triangle. The sign indicates triangle orientation
        @param vertices_idx: common two vertices index between the triangle and the potential new one
        @param new_point: third vertex of the new potential triangle
        @return:
        """
        vtx_0 = self.l_trgl_vertices_[vertices_idx[0]]
        vtx_1 = self.l_trgl_vertices_[vertices_idx[1]]
        vtx_2 = new_point

        det = (vtx_1.x_coord_ - vtx_0.x_coord_) * (vtx_2.y_coord_ - vtx_1.y_coord_) - (
                vtx_1.y_coord_ - vtx_0.y_coord_) * (vtx_2.x_coord_ - vtx_1.x_coord_)

        return det

    def choose_direction(self, mesh_trgl_orient, new_point):
        """
        Find the edge to continue the walk on triangles
        @param mesh_trgl_orient: triangles orientation in the mesh (All trianglea have the same orientation)
        @param new_point: Goal point of the walk in triangles
        @return:
        """
        l_couple_vrtx = [(1, 2), (2, 0), (0, 1)]

        for edge in l_couple_vrtx:
            new_orient = self.orientation_new_vrtx(edge, new_point)

            if mesh_trgl_orient * new_orient < 0:
                edge_id = l_couple_vrtx.index(edge)
                local_edge_id = self.triangle_idx_ * 3 + edge_id
                return local_edge_id

    def compute_surface(self, vertices_idx, point):
        """
        Compute triangle surface
        @param vertices_idx: common two vertices index between the triangle and the potential new one
        @param point: third vertex of the new potential triangle
        @return: triangle surface
        """
        vtx_0 = self.l_trgl_vertices_[vertices_idx[0]]
        vtx_1 = self.l_trgl_vertices_[vertices_idx[1]]
        vtx_2 = point

        coord_matrix = np.array(
            [[vtx_0.x_coord_, vtx_1.x_coord_, vtx_2.x_coord_], [vtx_0.y_coord_, vtx_1.y_coord_, vtx_2.y_coord_],
             [1, 1, 1]])

        surface = np.linalg.det(coord_matrix)

        return surface

    def determine_barycentric_coordinates(self, point):
        """
        Compute barycentric coordinates
        @param point: internal Point
        @return: list of barycentric coordinates [C0, C1, C2]
        """
        l_couple_vrtx = [(1, 2), (2, 0), (0, 1)]
        l_barycentric_coordinates = []
        for vertex in l_couple_vrtx:
            surf_partiel = self.compute_surface(vertex, point)

            coeff = surf_partiel / self.surf_

            l_barycentric_coordinates.append(coeff)

        return l_barycentric_coordinates

    def compute_gradient_matrix(self):
        """
        Compute the matrix in the triangle
        M = inv( (x2-x1)   (y2-y1)
                 (x3-x1)   (y3-y1))
        @return: matrix M
        """
        grad_mat = np.array([[self.l_trgl_vertices_[1].x_coord_ - self.l_trgl_vertices_[0].x_coord_,
                              self.l_trgl_vertices_[1].y_coord_ - self.l_trgl_vertices_[0].y_coord_],
                             [self.l_trgl_vertices_[2].x_coord_ - self.l_trgl_vertices_[0].x_coord_,
                              self.l_trgl_vertices_[2].y_coord_ - self.l_trgl_vertices_[0].y_coord_]])

        inv_grad_mat = np.linalg.inv(grad_mat)

        return inv_grad_mat

    def height(self, base):
        """
        compute triangle height length (from triangle area and base edge)
        @param base: base edge
        @return: height length
        """
        base_width = base.compute_length()
        return abs(2*self.surf_/base_width)

    def minimal_height(self):
        """
        compute triangle minimal height length
        @return: minimal height length
        """
        min_height=self.height(self.l_trgl_edges_[0])
        for edge in self.l_trgl_edges_[1:]:
            height=self.height(edge)
            if height<min_height:
                min_height=height

        return min_height

if __name__ == '__main__':
    print("Hello")

    from lib_data_structure.Point import Point2d

    mon_point1 = Point2d(1, 1, 0)
    mon_point2 = Point2d(2, 1.5, 1)
    mon_point3 = Point2d(3, 0.5, 1)
    l_pts = [mon_point1, mon_point2, mon_point3]

    mon_point = Point2d(4, 1, 0.2)

    from lib_data_structure.Edge import Edge

    edge1 = Edge(mon_point1, mon_point2, 1)
    edge2 = Edge(mon_point2, mon_point3, 2)
    edge3 = Edge(mon_point3, mon_point1, 3)

    l_edges = [edge1, edge2, edge3]

    mon_triangle = Triangle(l_pts, l_edges, 1)

    print(mon_triangle.determine_barycentric_coordinates(mon_point))
