"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """

from lib_data_structure.Point import Point2d
from lib_data_structure.Edge import Edge
from lib_data_structure.Triangle import Triangle
from read_write_vtk import write_vtk_with_sol
import gmsh

"""
Scripts to import data from .msh or .mesh file and to export the signed distance in a .sol file

"""


def open_gmsh_file(file_path):
    """
    Open mesh file with gmsh
    @param file_path: absolute path
    """
    gmsh.initialize([])
    gmsh.open(file_path)


def close_gmsh_file():
    """
    Close mesh file with gmsh
    """
    gmsh.finalize()


def import_nodes():
    """
    Import information about nodes
    @return: list of nodes (Point2d objects)
    """
    array_idx, array_coord, array_tag = gmsh.model.mesh.getNodes(dim=-1, tag=-1)
    l_vertices = [0] * len(array_idx)
    for vertex_id in range(len(array_idx)):
        x_coord_id = (int(vertex_id)) * 3
        new_point = Point2d(int(array_idx[vertex_id]), float(array_coord[x_coord_id]),
                            float(array_coord[x_coord_id + 1]))

        l_vertices[int(array_idx[vertex_id] - 1)] = new_point
    return l_vertices


def import_edges(l_vertices):
    """
    Import information about edges
    @param l_vertices: list of nodes (Point2d objects)
    @return: list of edges (Edge objects) and a dictionary which contains the point indexes as key and edge index as value
    """
    l_tag = gmsh.model.getEntities(dim=1)

    l_edges = []
    dct_point_edge = {}
    nb_edges = 0
    for tag in l_tag:
        array_dim, array_edge_idx, array_vertex_idx = gmsh.model.mesh.getElements(dim=1, tag=tag[1])

        for edge in range(0, len(array_edge_idx[0])):
            first_point_idx = int(array_vertex_idx[0][edge * 2]) - 1
            end_point_idx = int(array_vertex_idx[0][edge * 2 + 1]) - 1
            new_edge = Edge(l_vertices[first_point_idx], l_vertices[end_point_idx], nb_edges + 1, tag[1])

            dct_point_edge[(l_vertices[first_point_idx].point_idx_, l_vertices[end_point_idx].point_idx_)] = nb_edges
            dct_point_edge[(l_vertices[end_point_idx].point_idx_, l_vertices[first_point_idx].point_idx_)] = nb_edges

            l_edges.append(new_edge)
            nb_edges += 1
    return l_edges, dct_point_edge


def find_triangle_vertices(l_vertices, array_vertex_idx, trgl_id):
    """
    Build list of vertices for a triangle
    @param l_vertices: list of mesh vertices
    @param array_vertex_idx: list of all triangle vertices (gmsh output)
    @param trgl_id: triangle index
    @return: list of the three triangle vertices
    """
    trgl_vertices = []

    for vertex in range(0, 3):
        point_idx = int(array_vertex_idx[0][trgl_id * 3 + vertex]) - 1
        global_point = l_vertices[point_idx]
        trgl_vertex = Point2d(global_point.get_attributes()[0], global_point.get_attributes()[1],
                              global_point.get_attributes()[2])

        trgl_vertices.append(trgl_vertex)

    return trgl_vertices


def add_edge(dct_point_edge, l_edges, end_points, dct_global_local_edge_idx):
    """
    If the edge has not been declared before, create a new Edge object
    @param dct_point_edge: dictionary which contains the point indexes as key and edge index as value
    @param l_edges: list of mesh edges
    @param end_points: edge end points indexes
    @param dct_global_local_edge_idx: dictionary which contains the edge mesh index as key and edge local index(es) as value
    @return: edge global index
    """
    new_edge = Edge(end_points[0], end_points[1], len(l_edges) + 1, 0)

    l_edges.append(new_edge)

    end_points_idx = (end_points[0].point_idx_, end_points[1].point_idx_)

    dct_point_edge[(end_points_idx[0], end_points_idx[1])] = new_edge.edge_idx_ - 1
    dct_point_edge[(end_points_idx[1], end_points_idx[0])] = new_edge.edge_idx_ - 1

    dct_global_local_edge_idx[new_edge.edge_idx_] = []

    return new_edge.edge_idx_


def find_triangle_edges(trgl_vertices, dct_point_edge, l_edges, dct_global_local_edge_idx):
    """
    Build list of edges for a triangle
    @param trgl_vertices: list of the three triangle vertices
    @param dct_point_edge: dictionary which contains the point indexes as key and edge index as value
    @param l_edges: list of mesh edges
    @param dct_global_local_edge_idx: dictionary which contains the edge mesh index as key and edge local index(es) as value
    @return: list of the three triangle edges
    """
    trgl_edges = []

    for edge in range(0, 3):
        end_points = trgl_vertices.copy()
        end_points.pop(edge)

        end_points_idx = (end_points[0].point_idx_, end_points[1].point_idx_)

        if end_points_idx in dct_point_edge:
            edge_id = dct_point_edge[end_points_idx]
            tag = l_edges[edge_id].edge_tag_
            new_edge = Edge(end_points[0], end_points[1], edge_id + 1, tag)
        else:
            edge_id = add_edge(dct_point_edge, l_edges, end_points, dct_global_local_edge_idx) - 1
            new_edge = Edge(end_points[0], end_points[1], edge_id + 1, 0)

        trgl_edges.append(new_edge)

    return trgl_edges


def init_dct_global_local_edge_idx(l_edges):
    """
    Initialize dictionary which contains the edge mesh index as key and edge local index(es) as value
    @param l_edges: list of mesh edges
    @return: Initialized dictionary
    """
    dct_global_local_edge_idx = {}
    for edge in l_edges:
        dct_global_local_edge_idx[edge.edge_idx_] = []

    return dct_global_local_edge_idx


def update_dct_global_local_edge_idx(dct_global_local_edge_idx, new_triangle):
    """
    Update dictionary which contains the edge mesh index as key and edge local index(es) as value
    Add local index for each triangle edge
    @param dct_global_local_edge_idx: dictionary which contains the edge mesh index as key and edge local index(es) as value
    @param new_triangle: new triangle (Triangle object)
    """
    l_edges_trgl = new_triangle.l_trgl_edges_
    l_local_edge_index = new_triangle.l_local_edge_index_
    for edge in l_edges_trgl:
        id = l_edges_trgl.index(edge)
        dct_global_local_edge_idx[edge.edge_idx_].append(l_local_edge_index[id])


def import_triangles(l_vertices, dct_point_edge, l_edges):
    """
    Import information about triangles
    @param l_vertices: list of mesh vertices
    @param dct_point_edge: dictionary which contains the point indexes as key and edge index as value
    @param l_edges: list of mesh edges
    @return: list of mesh triangles (Triangle objects) and
            dictionary which contains the edge mesh index as key and edge local index(es) as value
    """
    l_reg = gmsh.model.getEntities(dim=2)
    l_triangles = []
    dct_global_local_edge_idx = init_dct_global_local_edge_idx(l_edges)
    nb_triangle = 0

    for reg in l_reg:

        array_dim, array_trgl_idx, array_vertex_idx = gmsh.model.mesh.getElements(dim=2, tag=reg[1])

        for trgl in range(0, len(array_trgl_idx[0])):
            trgl_vertices = find_triangle_vertices(l_vertices, array_vertex_idx, trgl)

            trgl_edges = find_triangle_edges(trgl_vertices, dct_point_edge, l_edges, dct_global_local_edge_idx)

            new_triangle = Triangle(trgl_vertices, trgl_edges, nb_triangle + 1, reg[1])
            update_dct_global_local_edge_idx(dct_global_local_edge_idx, new_triangle)
            l_triangles.append(new_triangle)
            nb_triangle += 1

    return l_triangles, dct_global_local_edge_idx


def import_all_data(file_path):
    """
    Import all mesh information : vertices, edges, triangles
    @param file_path: absolute mesh file path
    @return: list of mesh vertices,
            list of mesh edges,
            dictionary which contains the point indexes as key and edge index as value,
            list of mesh triangles,
            dictionary which contains the edge mesh index as key and edge local index(es) as value
    """
    open_gmsh_file(file_path)

    l_vertices = import_nodes()
    l_edges, dct_point_edge = import_edges(l_vertices)

    l_triangles, dct_global_local_edge_idx = import_triangles(l_vertices, dct_point_edge, l_edges)
    close_gmsh_file()
    return l_vertices, l_edges, dct_point_edge, l_triangles, dct_global_local_edge_idx


def open_file_to_write(filename):
    """
    Open .sol file
    @param filename: absolute .sol file path
    """
    my_file = open(filename, "w")
    return my_file


def export_data(sol_file_path, l_sol, l_unknowns):
    """
    Write .sol file which contains the signed distance function
    @param sol_file_path: absolute .sol file path
    @param l_sol: list of solutions (a single value for each vertex)
    @param l_unknowns: list of vertices index which represents the order of the solutions
    """
    sol_file = open_file_to_write(sol_file_path)
    sol_file.write("MeshVersionFormatted\n2\n\nDimension\n2\n\nSolAtVertices\n")
    sol_file.write(str(len(l_unknowns)))
    sol_file.write("\n1 1\n")

    for vertex in range(1, len(l_unknowns) + 1):
        vertex_pos = l_unknowns.index(vertex)
        sol_file.write(str(l_sol[vertex_pos]) + "\n")

    sol_file.write("End")

    sol_file.close()


def export_data_vtk(vtk_filepath_to_read, vtk_filepath_to_write, l_sol, l_unknowns):
    """

    @param vtk_filepath_to_read:
    @param vtk_filepath_to_write:
    @param l_sol:
    @param l_unknowns:
    """
    l_ordered_sol = []
    for vertex in range(1, len(l_unknowns) + 1):
        vertex_pos = l_unknowns.index(vertex)
        l_ordered_sol.append(l_sol[vertex_pos])
    write_vtk_with_sol(vtk_filepath_to_read, vtk_filepath_to_write, l_ordered_sol)


if __name__ == '__main__':
    l_vertices_ex, l_edges_ex, dct_point_edge_ex, l_triangles_ex, dct_global_local_edge_idx_ex = import_all_data(
        "/home/legentil/Bureau/Test_jupyter/Training_data/input_mesh.mesh")
    print(l_triangles_ex)
