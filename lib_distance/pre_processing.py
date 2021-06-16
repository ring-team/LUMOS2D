"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

This program is a Trade Secret of the ASGA and it is not to be:
 * reproduced, published, or disclosed to other,
 * distributed or displayed,
 * used for purposes or on Sites other than described in the
   RING-GOCAD Advancement Agreement Phase III,
   without the prior written authorization of the ASGA.

Licensee agrees to attach or embed this Notice on all copies of the program,
including partial copies or modified versions thereof.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
from lib_data_structure.Mesh import *
from lib_cavity.write_local_mesh_file import *
import numpy as np
from lib_data_structure.import_export_data import *

def read_sol(sol_path):
    """
    Extract the distances from a sol file
    @param sol_path: path to the sol file to read
    @return: list with the distance values of each point.
    They have the same index as in the mesh object
    """
    dist = []
    sol_file = open(sol_path, "r")
    lines = sol_file.readlines()
    for line in lines[9:]:  # Skipping the head
        if "End" in line:  # Checking for end of file
            break
        dist.append(float(line))
    return dist


def find_cut_triangles(path_to_sol, path_to_mesh, level_set, l_reg):
    """
    Function that define a list of triangle that are on the fault.
    It uses the sign of the vertices of triangles.
    @param path_to_sol: path to the sol file
    @param path_to_mesh: path to the mesh file
    @param level_set: iso value of the interface that is inserted
    @param l_reg: list of regions where the interface is inserted
    @return: list of triangles that will be impacted by the level-set discretization
    """

    mesh = Mesh(path_to_mesh)
    list_tri = [tri for tri in mesh.l_triangles_ if tri.triangle_reg_ in l_reg]
    dist = read_sol(path_to_sol)
    l_sol = [value - level_set for value in dist]

    frac = []
    eps = 0.0005

    for triangle in list_tri:

        nplus = 0  # Point in the triangle that are positive
        nminus = 0  # Point on the triangle that negative
        for vertex in triangle.l_trgl_vertices_idx_:
            if l_sol[vertex-1] >= eps:
                nplus += 1
            else:
                nminus += 1
        if (nminus == 1 or nplus == 1) and boundary_triangle(mesh, triangle):  # If the triangle have 2 signs it is on the fault
            frac.append(triangle)

    return frac


def boundary_triangle(mesh, triangle):
    """
    Detect boundary triangles (triangles that have an edge representing a limit of a region)
    @param mesh: Mesh object
    @param triangle: Triangle to test
    @return: True if the triangle have a boundary edge, False otherwise
    """
    l_adja = mesh.find_adjacent_triangles(triangle)
    for tri in l_adja:
        if triangle.triangle_reg_!=tri.triangle_reg_:
            return True
    return False


def compute_gradient(l_sol, triangle):
    """
    Compute the gradient of the scalar field in a triangle
    @param l_sol: list with the distance values of each point.
    @param triangle: Triangle where the gradient is computed
    @return: gradient norm
    """
    grad_triangle = triangle.compute_gradient_matrix()
    vertex_1 = triangle.l_trgl_vertices_idx_[0]
    vertex_2 = triangle.l_trgl_vertices_idx_[1]
    vertex_3 = triangle.l_trgl_vertices_idx_[2]
    mat_scalar_field = np.array([[l_sol[vertex_2-1]-l_sol[vertex_1-1]],[l_sol[vertex_3-1]-l_sol[vertex_1-1]]])
    grad_scalar_field = np.dot(grad_triangle, mat_scalar_field)
    return np.linalg.norm(grad_scalar_field)


def snap_point(dist_init, l_cut_triangles, epsilon, final_sol_path):
    """
    Modify the scalar field according to the rules explained in Legentil et al (2021)
    @param dist_init: list of scalar field values from .sol file
    @param l_cut_triangles: list of triangles that will be impacted by the level-set discretization
    @param epsilon: tolerance assmilated to a distance
    @param final_sol_path:  final .sol file path with the modified scalar field property
    @return:
    """
    l_sol = [value - level_set for value in dist_init]

    for triangle in l_cut_triangles:
        L_vertices = triangle.l_trgl_vertices_idx_
        L_vertices.append(triangle.l_trgl_vertices_idx_[0])
        for vertex in range (0, 2):
            if l_sol[L_vertices[vertex]-1]*l_sol[L_vertices[vertex+1]-1]<0:
                norm_grad_scalar_field=compute_gradient(l_sol, triangle)
                dist = min(abs(l_sol[L_vertices[vertex]-1]), abs(l_sol[L_vertices[vertex+1]-1]))
                if dist==abs(l_sol[L_vertices[vertex]-1]):
                    if (dist / norm_grad_scalar_field) < epsilon:
                        l_sol[L_vertices[vertex]-1] = 0

                else:
                    if (dist / norm_grad_scalar_field) < epsilon:
                        l_sol[L_vertices[vertex+1]-1] = 0


    l_sol_final = [value + level_set for value in l_sol]
    l_id=[n for n in range(1,len(l_sol_final)+1)]
    export_data(final_sol_path, l_sol_final, l_id)



if __name__ == '__main__':
    mesh = "D:/Article/simu_05_01_21/just_mesh/article_model_f.mesh"
    sol_path = "D:/Article/simu_05_01_21/just_mesh/article_model_f.sol"
    level_set = -590
    final_sol_path = sol_path[:-4] + '_' + str(level_set) + ".sol"

    l_sol = read_sol(sol_path)
    l_cut_triangle = find_cut_triangles(sol_path, mesh, level_set, [6])
    print([triangle.triangle_idx_ for triangle in l_cut_triangle])
    
    snap_point(l_sol, l_cut_triangle, 0.5, final_sol_path)