"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
from lib_data_structure.Mesh import Mesh
from lib_data_structure.Line import Line
import numpy as np
# from scipy.linalg import pinv2
from lib_data_structure.import_export_data import export_data
import shapely.geometry as shp
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import lsqr


def fill_equation_adja_matrix(eq_matrix, l_unknowns, eq_coeff, l_unknowns_trgl, eq_id):
    """
    Fill equation system AX=B with equation from triangle adjacencies
    @param eq_matrix: matrix A
    @param l_unknowns: vertices indexes (equation order in the matrices)
    @param eq_coeff: coeff to add to eq_matrix
    @param l_unknowns_trgl: indexes of triangle vertices
    @param eq_id: index of the new equation
    """
    eq_coeff_norm = [i / max(eq_coeff) for i in eq_coeff]
    if l_unknowns_trgl[-1] not in l_unknowns:
        l_unknowns.append(l_unknowns_trgl[-1])
    for coeff in range(0, 4):
        col_id = l_unknowns.index(l_unknowns_trgl[coeff])
        eq_matrix[eq_id, col_id] += eq_coeff_norm[coeff]


def solve_system(eq_matrix, res, l_unknowns):
    """
    Solve equation system AX=B
    @param eq_matrix: matrix A
    @param res: matrix B
    @param l_unknowns: vertices indexes (equation order in the matrices)
    @return: sol : the solution, the vector X
             l_unknowns : unchanged (vertices indexes (equation order in the matrices))
    """
    print("inversion")
    sparse_matrix = csr_matrix(eq_matrix)
    sol = lsqr(sparse_matrix, res)[0]
    print(sol)
    return sol, l_unknowns


def compute_gradient_matrix(triangle):
    """
    Compute the gradient matrix and return the expression as a function of f1, f2, f3 (scalar field value at triangle vertex)
      _
     \/f = (grad x  = inv( (x2-x1)   (y2-y1)  .  ((f2-f1)
            grad y)        (x3-x1)   (y3-y1))     (f3-f1))
    @param triangle: Triangle object
    @return: l_grad_x : coefficient of grad x
             l_grad_y : coefficient of grad x
    """

    grad_mat_inv = triangle.compute_gradient_matrix()
    l_grad_x = np.zeros(3)
    l_grad_x[0] = -(grad_mat_inv[0][0] + grad_mat_inv[0][1])
    l_grad_x[1] = grad_mat_inv[0][0]
    l_grad_x[2] = grad_mat_inv[0][1]

    l_grad_y = np.zeros(3)
    l_grad_y[0] = -(grad_mat_inv[1][0] + grad_mat_inv[1][1])
    l_grad_y[1] = grad_mat_inv[1][0]
    l_grad_y[2] = grad_mat_inv[1][1]

    return l_grad_x, l_grad_y


def update_unknown_and_equations(dct_trgl_coeff, l_unknown, triangle, coeff):
    """
    Update equation system AX=B
    @param dct_trgl_coeff: dictionary that contains Triangle object as key and the barycentric coefficients
                                 describing the line point (that is inside the Triangle)
    @param l_unknown: vertices indexes (equation order in the matrices)
    @param triangle: Triangle
    @param coeff: barycentric coefficient to add
    """
    if triangle not in dct_trgl_coeff:
        dct_trgl_coeff[triangle] = [coeff]
    else:
        dct_trgl_coeff[triangle].append(coeff)

    for vrtx in triangle.l_trgl_vertices_:
        if vrtx.point_idx_ not in l_unknown:
            l_unknown.append(vrtx.point_idx_)


def diff_adja_gradient(triangle1, triangle2):
    """
    Compute the difference between gradients of two triangles
    \/f = inv( (x2-x1)   (y2-y1)  .  (f2-f1
                   (x3-x1)   (y3-y1))     f3-f1)
    diff= (grad x = \/f1-\/f2
           grad y)
    @param triangle1: first Triangle
    @param triangle2: second Triangle adjacent to triangle1
    @return: diff_l_grad_x: coefficient of grad x as a function of f1, f2, f3
             diff_l_grad_y: coefficient of grad y as a function of f1, f2, f3
             l_unknowns: vertices indexes (equation order in the matrices)
    """
    diff_l_grad_x = np.zeros(4)
    diff_l_grad_y = np.zeros(4)

    l_grad_x1, l_grad_y1 = compute_gradient_matrix(triangle1)
    diff_l_grad_x[:3] = l_grad_x1
    diff_l_grad_y[:3] = l_grad_y1
    l_unknowns = [i for i in triangle1.l_trgl_vertices_idx_]
    l_grad_x2, l_grad_y2 = compute_gradient_matrix(triangle2)

    for coeff in range(0, 3):
        if triangle2.l_trgl_vertices_idx_[coeff] in triangle1.l_trgl_vertices_idx_:
            id_in_trgl1 = triangle1.l_trgl_vertices_idx_.index(triangle2.l_trgl_vertices_idx_[coeff])
            diff_l_grad_x[id_in_trgl1] -= l_grad_x2[coeff]
            diff_l_grad_y[id_in_trgl1] -= l_grad_y2[coeff]
        else:
            diff_l_grad_x[3] = -l_grad_x2[coeff]
            diff_l_grad_y[3] = -l_grad_y2[coeff]
            l_unknowns.append(triangle2.l_trgl_vertices_idx_[coeff])

    return diff_l_grad_x, diff_l_grad_y, l_unknowns


def fill_equation_matrix(eq_matrix, triangle, eq_id, eq_coeff, l_unknowns):
    """
    Fill equation system AX=B with equation
    @param eq_matrix: matrix A
    @param triangle: Triangle concerned by the new equation
    @param eq_id: index of new equation
    @param eq_coeff: coefficients to add in A
    @param l_unknowns: vertices indexes (equation order in the matrices)
    """
    eq_coeff_norm = [i / max(eq_coeff) for i in eq_coeff]
    for coeff in range(0, 3):
        col_id = l_unknowns.index(triangle.l_trgl_vertices_[coeff].point_idx_)
        eq_matrix[eq_id][col_id] = eq_coeff_norm[coeff]


class SDF:
    """
       A class to represent an edge.

       @Attributes
       mesh_ : mesh on which the Signed Distance Function is computed
       line_ : polygonal line (the scalar field represents the distance to this line)
       output_file_ : file to save the solution (.sol)
       sol_ : list of value at each vertex

       @Methods
       find_first_point_in_boundary_box(): Find the first point that inside the mesh
       in_the_boundary_box(point_id): Test if a point is inside the mesh
       intersected_triangle(): Find triangles intersected by the level-set and update the equation system with
                               information in these triangles. At each line vertex the scalar field is equal to 0.
       normal_eq(triangle, line_point_id): Equation to impose gradient direction
       normal_adjacent_eq(triangle1, triangle2): Determine equation coefficients to represents the continuity of the
                                                 scalar field between two adjacent triangles
       init_matrix(dct_trgl_coeff, l_unknowns): Initialization of all matrices of equation system (AX=B)
        build_matrix_system_intersected_triangles(dct_trgl_coeff, l_unknowns): Fill the coefficient matrix A and B of
                                                                the system (AX=B) with intersected triangles information
       build_matrix_system(dct_trgl_coeff, l_unknowns): Build the equation system with matrix: AX=B
       compute_SDF(): Main function to compute the signed Distance function

    """

    def __init__(self, mesh_file_path, line_file_path, output_file_path):
        self.mesh_ = Mesh(mesh_file_path)
        self.line_ = Line(line_file_path)
        self.output_file_ = output_file_path
        self.sol_ = []

    def find_first_point_in_boundary_box(self):
        """
        Find the first point that inside the mesh
        @return: point_id: index of the first line point that is in the mesh
        """
        boundary_polygon = self.mesh_.boundary_polygon_
        first_point = shp.Point(self.line_.l_vertices_[0].get_coordinates())
        point_id = 0

        while (not first_point.intersects(boundary_polygon)) and point_id < len(self.line_.l_vertices_):
            point_id += 1
            first_point = shp.Point(self.line_.l_vertices_[point_id].get_coordinates())

        return point_id

    def in_the_boundary_box(self, point_id):
        """
        Test if a point is inside the mesh
        @param point_id: line point index
        @return: True if the point is inside the mesh, False otherwise
        """
        point = shp.Point(self.line_.l_vertices_[point_id].get_coordinates())

        return point.intersects(self.mesh_.boundary_polygon_)

    def intersected_triangle(self):
        """
        Find triangles intersected by the level-set and update the equation system with information in these triangles
        At each line vertex the scalar field is equal to 0.
        @return: dct_trgl_coeff: dictionary that contains Triangle object as key and the barycentric coefficients
                                 describing the line point (that is inside the Triangle)
                 l_unknows: vertices indexes (equation order in the second matrix)
        """
        dct_trgl_coeff = {}
        l_unknowns = []

        first_point_in_boundary_box_id = self.find_first_point_in_boundary_box()
        first_triangle = self.mesh_.find_in_triangle(self.line_.l_vertices_[first_point_in_boundary_box_id])
        coeff = first_triangle.determine_barycentric_coordinates(self.line_.l_vertices_[first_point_in_boundary_box_id])
        coeff.append(first_point_in_boundary_box_id)
        update_unknown_and_equations(dct_trgl_coeff, l_unknowns, first_triangle, coeff)
        triangle = first_triangle
        last_point_in_box_idx = first_point_in_boundary_box_id
        for point_idx in range(first_point_in_boundary_box_id + 1, len(self.line_.l_vertices_)):
            if self.in_the_boundary_box(point_idx):

                if self.line_.l_vertices_[point_idx].in_triangle(triangle):
                    coeff = triangle.determine_barycentric_coordinates(self.line_.l_vertices_[point_idx])

                else:
                    triangle = self.mesh_.walk_on_triangle(self.line_.l_vertices_[last_point_in_box_idx],
                                                           self.line_.l_vertices_[point_idx])
                    coeff = triangle.determine_barycentric_coordinates(self.line_.l_vertices_[point_idx])
                last_point_in_box_idx = point_idx
                coeff.append(point_idx)
                update_unknown_and_equations(dct_trgl_coeff, l_unknowns, triangle, coeff)

        return (dct_trgl_coeff, l_unknowns)

    def normal_eq(self, triangle, line_point_id):
        """
        Equation to impose gradient direction
         _
        \/f = inv( (x2-x1)   (y2-y1)  .  (f2-f1
                   (x3-x1)   (y3-y1))     f3-f1)

        We want \/f . n = 1
        @param triangle: triangle intersected by the level-set and containing a line point
        @param line_point_id: index of the line point
        @return: eq_coeff: equation coefficients
        """
        l_grad_x, l_grad_y = compute_gradient_matrix(triangle)

        norm_vect = self.line_.l_normal_verctors_[line_point_id]
        eq_coeff = l_grad_x * norm_vect[0] + l_grad_y * norm_vect[1]

        return eq_coeff

    def normal_adjacent_eq(self, triangle1, triangle2):
        """
        Determine equation coefficients to represents the continuity of the scalar field between two adjacent triangles
         _
        \/f = inv( (x2-x1)   (y2-y1)  .  (f2-f1
                   (x3-x1)   (y3-y1))     f3-f1)

        We want (\/f1 -\/f2). n = 0
        @param triangle1: first Triangle
        @param triangle2: second Triangle adjacent to triangle1
        @return: eq_coeff_norm: equation coefficients in triangle2
                 l_unknowns: vertices order (=equation order in matrices)
        """
        l_grad_x, l_grad_y, l_unknowns = diff_adja_gradient(triangle1, triangle2)

        common_edge_vertices = triangle1.find_common_edge(triangle2)
        edge = self.mesh_.dct_point_edge_[common_edge_vertices]

        norm_vect = self.mesh_.l_edges_[edge].compute_normal_vector()
        eq_coeff = l_grad_x * norm_vect[0] + l_grad_y * norm_vect[1]
        eq_coeff_norm = [i / max(eq_coeff) for i in eq_coeff]

        return eq_coeff_norm, l_unknowns

    def init_matrix(self, dct_trgl_coeff, l_unknowns):
        """
        Initialization of all matrices of equation system (AX=B)
        @param dct_trgl_coeff: dictionary that contains Triangle object as key and the barycentric coefficients
                                 describing the line point (that is inside the Triangle)
        @param l_unknowns: vertices order (=equation order in matrices)
        @return: eq_matrix: matrix A (sparse matrix)
                 res: matrix B
                 eq_nb_intersection - 1: number of triangles that are intersected by the line
        """
        unknown_nb = len(self.mesh_.l_vertices_)
        equation_nb = len(self.line_.l_vertices_) + len(self.mesh_.l_triangles_)

        eq_matrix = np.zeros((equation_nb, unknown_nb))
        res = np.zeros(equation_nb)

        eq_matrix_intersection, res_intersection = self.build_matrix_system_intersected_triangles(dct_trgl_coeff,
                                                                                                  l_unknowns)
        eq_nb_intersection, unknown_nb_intersect = eq_matrix_intersection.shape

        eq_matrix[:eq_nb_intersection, :unknown_nb_intersect] = eq_matrix_intersection
        res[:eq_nb_intersection] = res_intersection

        return lil_matrix(eq_matrix), res, eq_nb_intersection - 1

    def build_matrix_system_intersected_triangles(self, dct_trgl_coeff, l_unknowns):
        """
        Fill the coefficient matrix A and B of the system (AX=B) with intersected triangles information
        @param dct_trgl_coeff: dictionary that contains Triangle object as key and the barycentric coefficients
                                 describing the line point (that is inside the Triangle)
        @param l_unknowns: vertices order (=equation order in matrices)
        @return: eq_matrix: matrix A filled
                 res: matrix B

        """
        unknown_nb = len(l_unknowns)
        equation_nb = len(self.line_.l_vertices_) + len(dct_trgl_coeff)

        eq_matrix = np.zeros((equation_nb, unknown_nb))
        res = np.zeros(equation_nb)
        eq_id = 0
        for triangle in dct_trgl_coeff:
            for barycentric_coeff in dct_trgl_coeff[triangle]:
                fill_equation_matrix(eq_matrix, triangle, eq_id, barycentric_coeff, l_unknowns)
                eq_id += 1
            norm_eq_coeff = self.normal_eq(triangle, barycentric_coeff[3])
            fill_equation_matrix(eq_matrix, triangle, eq_id, norm_eq_coeff, l_unknowns)
            res[eq_id] = 1 / max(norm_eq_coeff)
            eq_id += 1
        return eq_matrix, res

    def build_matrix_system(self, dct_trgl_coeff, l_unknowns):
        """
        Build the equation system with matrix: AX=B
        @param dct_trgl_coeff: dictionary that contains Triangle object as key and the barycentric coefficients
                                 describing the line point (that is inside the Triangle)
        @param l_unknowns: vertices order (=equation order in matrices)
        @return: eq_matrix: matrix A filled
                 res: matrix B
        """
        eq_matrix, res, eq_id = self.init_matrix(dct_trgl_coeff, l_unknowns)

        l_remaining_triangles = [self.mesh_.l_triangles_[id] for id in range(len(self.mesh_.l_triangles_)) if
                                 self.mesh_.l_triangles_[id] not in dct_trgl_coeff]
        l_processed_triangles = [trgl for trgl in dct_trgl_coeff]

        while l_remaining_triangles:
            l_triangles_traites_new = []
            for triangle in l_processed_triangles:
                l_adj_trgl = self.mesh_.find_adjacent_triangles(triangle)
                for adjacent in l_adj_trgl:
                    if adjacent in l_remaining_triangles:
                        eq_id += 1
                        eq_coeff, l_unknowns_trgl = self.normal_adjacent_eq(triangle, adjacent)
                        fill_equation_adja_matrix(eq_matrix, l_unknowns, eq_coeff, l_unknowns_trgl, eq_id)
                        l_triangles_traites_new.append(adjacent)
                        l_remaining_triangles.remove(adjacent)

            l_processed_triangles = l_triangles_traites_new

        return eq_matrix, res

    def compute_SDF(self):
        """
        main function to compute the signed Distance function

        """
        dct_trgl_coeff, l_unknowns = self.intersected_triangle()
        M, res = self.build_matrix_system(dct_trgl_coeff, l_unknowns)
        l_sol, l_unknowns = solve_system(M, res, l_unknowns)
        print("writing")
        export_data(self.output_file_, l_sol, l_unknowns)


if __name__ == '__main__':
    # my_mesh = "D:/Donnees_Linux/mesh_merge1234_2d.mesh"
    # my_line = "D:/Donnees_Linux/water_oil_line.msh"
    # output_file_path = "D:/Donnees_Linux/mesh_merge1234_2d.sol"
    my_mesh = "/home/legentil/Bureau/Test_jupyter/Training_data/frac2/mesh_frac1.mesh"
    my_line = "/home/legentil/Bureau/Test_jupyter/Training_data/frac1.msh"
    output_file_path = "/home/legentil/Bureau/Test_jupyter/Training_data/frac.sol"
    # my_mesh = "D:/Data/_2D/Simple_Models/A1/surf1/Code_Test/ext_line/curve1_cut.mesh"
    # my_line = "D:/ligne.msh"
    # output_file_path = "D:/test_ext.sol"
    my_SDF = SDF(my_mesh, my_line, output_file_path)

    my_SDF.compute_SDF()
