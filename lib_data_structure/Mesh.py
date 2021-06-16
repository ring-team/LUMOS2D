"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """

from lib_data_structure.import_export_data import import_all_data
import shapely.geometry as shp


class Mesh:
    """
       A class to represent a triangular mesh.

       @Attributes
       msh_file_path_ : file path to mesh file
       l_vertices_ : list of mesh vertices (list of Point)
       l_edges_ : list of mesh edges (list of Edge)
       dct_point_edge_ : dictionary that contains the point indexes as key and edge index as value
       l_triangles_ : list of mesh triangles (list of Triangle)
       dct_global_local_edge_idx_ : dictionary that contains the edge mesh index as key and edge local index(es)
                                            as value
       mesh_trgl_orientation_ : triangle orientation (positive or negative determinant)
       adjacencies_tab_ : triangle adjacencies table. For each mesh triangle, its adjacent triangle(s) is(are) stored
       l_edges_tag_ : list of edges tag (one tag for each edge)
       boundary_polygon_: shapely object that represents the boundary box (Polygon)
       l_geo_region_: inner and outer region
       tag_region_boundary_edges_: list of edges tag

       @Methods
       import_mesh(): define local edges indexes in the mesh
       build_adjacencies(): Build adjacencies table. For each mesh triangle, its adjacent triangle(s) is(are) stored.
       find_region_boundary_edges(): Change the boundary edge tags (add 100)
       l_contact(vertices_idx, new_point): Return possible pairs of adjacent region between inner and outer mesh parts
       search_following_edge_tuple(mesh_trgl_orient, new_point): Find edge that have a point in common (adjacent edge)
       build_sorted_list_boundary_point(l_edge_end_point): Order the list of points, to build the list of ordered
                                                                 boundary points (based on edge adjacencies)
       build_boundary_polygon(): Build bounding box
       find_first_triangle(point): Find the triangle where the point is in.
       find_next_triangle(triangle, end_point): (In walk on triangle algorithm) find the adjacent relevant triangle index
       walk_on_triangle(first_point, final_point): Find the triangle where the final point is in using the walk on
                                                   triangle algorithm. The first point is given to initialize the walk
       find_adjacent_triangles(triangle): Find all adjacent triangles to a triangle
       find_min_height(): Compute the minimal height of mesh triangles

    """

    def __init__(self, msh_file_path, l_geo_region=[]):
        self.msh_file_path_ = msh_file_path
        self.import_mesh()
        self.mesh_trgl_orientation_ = self.l_triangles_[0].orientation()
        self.adjacencies_tab_ = [0] * (3 * len(self.l_triangles_))
        self.l_edges_tag_ = [0] * len(self.adjacencies_tab_)
        self.build_adjacencies()
        self.boundary_polygon_ = self.build_boundary_polygon()
        self.l_geo_region_ = l_geo_region
        if self.l_geo_region_:
            self.tag_region_boundary_edges_ = self.change_tag_boundary_edges()

    def import_mesh(self):
        """
        Import all mesh information from mesh file : (lists of vertices, edges, triangles, and dictionaries that
        contain the point indexes as key and edge index as value and
        the edge mesh index as key and edge local index(es) as value
        """
        self.l_vertices_, self.l_edges_, self.dct_point_edge_, self.l_triangles_, self.dct_global_local_edge_idx_ = \
            import_all_data(self.msh_file_path_)

    def build_adjacencies(self):
        """
        Build adjacencies table. For each mesh triangle, its adjacent triangle(s) is(are) stored.
        """
        l_edges_tag = [edge.edge_tag_ for edge in self.l_edges_]

        for edge in self.dct_global_local_edge_idx_:
            if len(self.dct_global_local_edge_idx_[edge]) == 2:  # edge isn't boundary
                self.adjacencies_tab_[self.dct_global_local_edge_idx_[edge][0] - 3] = \
                    self.dct_global_local_edge_idx_[edge][1]
                self.l_edges_tag_[self.dct_global_local_edge_idx_[edge][0] - 3] = l_edges_tag[edge - 1]
                self.adjacencies_tab_[self.dct_global_local_edge_idx_[edge][1] - 3] = \
                    self.dct_global_local_edge_idx_[edge][0]
                self.l_edges_tag_[self.dct_global_local_edge_idx_[edge][1] - 3] = l_edges_tag[edge - 1]

            elif len(self.dct_global_local_edge_idx_[edge]) == 1:
                self.l_edges_tag_[self.dct_global_local_edge_idx_[edge][0] - 3] = l_edges_tag[edge - 1]


    def l_contact(self):
        """
        Return possible pairs of adjacent region between inner and outer mesh parts
        @return: list of index pair
        ex: inner mesh part : tag 1
            outer mesh part : tags 2 and 3
            return : [(1,2), (2,1), (1,3), (3,1)]
        """
        l_contact = []
        for reg1 in self.l_geo_region_[0]:
            for reg2 in self.l_geo_region_[1]:
                l_contact.append((reg1, reg2))
                l_contact.append((reg2, reg1))
        return l_contact

    def change_tag_boundary_edges(self):
        """
        Change the boundary edge tags (add 100)
        @return: list of edges tag
        """
        l_edges_tag = self.l_edges_tag_
        l_contact = self.l_contact()
        for edge_id in range(len(self.adjacencies_tab_)):

            triangle = self.l_triangles_[(edge_id // 3)]
            if l_edges_tag[edge_id] < 100 and self.adjacencies_tab_[edge_id] != 0:

                adj_triangle = self.l_triangles_[(self.adjacencies_tab_[edge_id] // 3) - 1]
                if (triangle.triangle_reg_, adj_triangle.triangle_reg_) in l_contact:
                    l_edges_tag[edge_id] += 100
                    l_edges_tag[self.adjacencies_tab_[edge_id] - 3] += 100

            else:
                continue

        return l_edges_tag

    def find_following_edge_tuple(self, l_edge_end_point, common_point):
        """
        Find edge that have a point in common (adjacent edge)
        @param l_edge_end_point: list of pairs of points
        @param common_point: point in common
        @return: the adjacent edge endpoints and the other end point (different of the common point)
        """
        for edge_tuple in l_edge_end_point:

            if common_point in edge_tuple:

                if common_point == edge_tuple[0]:
                    return edge_tuple, edge_tuple[1]

                else:
                    return edge_tuple, edge_tuple[0]

    def build_sorted_list_boundary_point(self, l_edge_end_point):
        """
        Order the list of points, to build the list of ordered boundary points (based on edge adjacencies)
        @param l_edge_end_point: list of boundary points
        @return: list of ordered boundary points
        """
        l_boundary_point = [l_edge_end_point[0][0]]
        other_point = l_edge_end_point[0][1]
        l_edge_end_point.remove(l_edge_end_point[0])
        for point in range(len(l_edge_end_point)):
            l_boundary_point.append(other_point)
            tuple, other_point = self.find_following_edge_tuple(l_edge_end_point, other_point)
            l_edge_end_point.remove(tuple)

        return l_boundary_point

    def build_boundary_polygon(self):
        """
        Build bounding box
        @return: shapely object that represents the boundary box (Polygon)
        """
        l_edge_end_point = [self.l_edges_[edge - 1].get_point_idx() for edge in self.dct_global_local_edge_idx_ if
                            len(self.dct_global_local_edge_idx_[edge]) == 1]
        l_boundary_point_idx = self.build_sorted_list_boundary_point(l_edge_end_point)
        l_boundary_point = [self.l_vertices_[point_idx - 1] for point_idx in l_boundary_point_idx]
        coords = [point.get_coordinates() for point in l_boundary_point]
        boundary_polygon = shp.Polygon(coords)

        return boundary_polygon

    def find_in_triangle(self, point):
        """
        Find the triangle where the point is in.
        @param point: Point
        @return: the Triangle where the point is in
        """
        for trgl in self.l_triangles_:

            in_triangle = point.in_triangle(trgl)
            if in_triangle:
                return trgl

    def find_next_triangle(self, triangle, next_point):
        """
        (In walk on triangle algorithm) find the adjacent relevant triangle index.
        @param triangle: previous triangle
        @param next_point: the goal of the walk on triangle
        @return: next triangle index
        """
        edge_local_id = triangle.choose_direction(self.mesh_trgl_orientation_, next_point)

        adjacent_edge = self.adjacencies_tab_[edge_local_id - 3]

        adjacent_trgl_id = adjacent_edge // 3

        return adjacent_trgl_id

    def walk_on_triangle(self, first_point, final_point):
        """
        Find the triangle where the final point is in using the walk on triangle algorithm. The first point is given to
        initialize the walk
        @param first_point: initial point (Point)
        @param final_point: final point (Point)
        @return: the triangle where the final point is in (Triangle)
        """
        triangle = self.find_in_triangle(first_point)
        in_triangle = final_point.in_triangle(triangle)

        while not in_triangle:
            next_trgl_id = self.find_next_triangle(triangle, final_point)
            triangle = self.l_triangles_[next_trgl_id - 1]
            in_triangle = final_point.in_triangle(triangle)

        return triangle

    def find_adjacent_triangles(self, triangle):
        """
        Find all adjacent triangles to a triangle
        @param triangle: Triangle object
        @return: list of adjacent triangles
        """
        triangle_id = triangle.triangle_idx_
        l_adj_triangles = []

        for edge_id in range(0, 3):
            local_edge_id = 3 * triangle_id + edge_id
            adjacent_edge = self.adjacencies_tab_[local_edge_id - 3]

            if adjacent_edge != 0:
                adjacent_trgl_id = adjacent_edge // 3
                l_adj_triangles.append(self.l_triangles_[adjacent_trgl_id - 1])

        return l_adj_triangles

    def find_min_height(self):
        """
        Compute the minimal height of mesh triangles
        @return: minimal height
        """
        min_height=self.l_triangles_[0].minimal_height()
        minheight_trgl_idx = 0
        for triangle in self.l_triangles_[1:]:
            height = triangle.minimal_height()
            if height<min_height:
                min_height=height
                minheight_trgl_idx = triangle.triangle_idx_
        return min_height, minheight_trgl_idx


if __name__ == '__main__':
    my_mesh = Mesh("/home/legentil/Bureau/Test_jupyter/Training_data/frac2/mesh_frac1.mesh")
    print(my_mesh.adjacencies_tab_)

    """from lib_data_structure.Point import Point2d

    first_point = Point2d(12, 1.5, -1.5)
    final_point = Point2d(13, 4, 1.5)

    my_mesh.walk_on_triangle(first_point, final_point)

    my_mesh2 = Mesh("/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/ext_line/region/curve1_cut.o.mesh", [[4],[1, 3]])
"""
