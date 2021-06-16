"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.


 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """


class Point:
    """
       A class to represent a point.

       @Attributes
       x_coord_ : x coordinate
       y_coord_ : y coordinate
       z_coord_ : z coordinate
       point_idx_ : point index

       @Methods
       define_local_point_idx(val): define local point index in the mesh
       """
    def __init__(self, point_idx, x_coord, y_coord, z_coord):
        self.x_coord_ = x_coord
        self.y_coord_ = y_coord
        self.z_coord_ = z_coord
        self.point_idx_ = point_idx

    def define_local_point_idx(self, val):
        self.local_point_idx_ = val


class Point2d(Point):
    """
           A class to represent a point.

           @Attributes
           x_coord_ : x coordinate
           y_coord_ : y coordinate
           point_idx_ : point index

           @Methods
           define_local_point_idx(val): define local mesh index
           in_triangle(triangle): Determine if a point is in the triangle
           """

    def __init__(self, point_idx, x_coord, y_coord):
        Point.__init__(self, point_idx, x_coord, y_coord, 0)

    def get_coordinates(self):
        return self.x_coord_, self.y_coord_

    def get_attributes(self):
        return self.point_idx_, self.x_coord_, self.y_coord_

    def in_triangle(self, trgl):
        '''
        Determine if a point is in the triangle
        :param trgl: triangle (Triangle object)

        c1 = (x2-x1)*(yp-y1)-(y2-y1)*(xp-x1)
        c2 = (x3-x2)*(yp-y2)-(y3-y2)*(xp-x2)
        c3 = (x1-x3)*(yp-y3)-(y1-y3)*(xp-x3)
        '''

        l_vertices_coordinates = trgl.get_vertices_coordinates()  ## [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]

        c1 = (l_vertices_coordinates[1][0] - l_vertices_coordinates[0][0]) * (
                    self.y_coord_ - l_vertices_coordinates[0][1]) - (
                     l_vertices_coordinates[1][1] - l_vertices_coordinates[0][1]) * (
                         self.x_coord_ - l_vertices_coordinates[0][0])
        c2 = (l_vertices_coordinates[2][0] - l_vertices_coordinates[1][0]) * (
                    self.y_coord_ - l_vertices_coordinates[1][1]) - (
                     l_vertices_coordinates[2][1] - l_vertices_coordinates[1][1]) * (
                         self.x_coord_ - l_vertices_coordinates[1][0])
        c3 = (l_vertices_coordinates[0][0] - l_vertices_coordinates[2][0]) * (
                    self.y_coord_ - l_vertices_coordinates[2][1]) - (
                     l_vertices_coordinates[0][1] - l_vertices_coordinates[2][1]) * (
                         self.x_coord_ - l_vertices_coordinates[2][0])

        if (c1 <= 0 and c2 <= 0 and c3 <= 0) or (c1 >= 0 and c2 >= 0 and c3 >= 0):
            return True
        else:
            return False


if __name__ == '__main__':
    print("Hello")
    mon_point1 = Point2d(1, 0, 0)
    mon_point2 = Point2d(2, 1, 0)
    mon_point3 = Point2d(3, 0.5, 1)
    l_pts = [mon_point1, mon_point2, mon_point3]

    mon_point = Point2d(4, 0.5, 1)

    from lib_data_structure.Edge import Edge

    edge1 = Edge(mon_point1, mon_point2, 1)
    edge2 = Edge(mon_point2, mon_point3, 2)
    edge3 = Edge(mon_point3, mon_point1, 3)

    l_edges = [edge1, edge2, edge3]

    from lib_data_structure.Triangle import Triangle

    mon_triangle = Triangle(l_pts, l_edges, 1)

    print(mon_point.in_triangle(mon_triangle))
