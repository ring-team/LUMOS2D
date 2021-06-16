"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
from lib_data_structure.Mesh import Mesh
import gmsh


def open_gmsh_file(file_path):
    """
    Open mesh file with gmsh
    @param file_path: absolute path
    """
    gmsh.initialize([])
    gmsh.open(file_path)


def get_triangle_tags(global_mesh_filepath):
    """
    Read the region tags that are used in a mesh
    @param global_mesh_filepath: path to the mesh file
    @return: list of region tags
    """
    gmsh.initialize([])
    gmsh.open(global_mesh_filepath)
    l_tuple_tags = gmsh.model.getEntities(dim=2)
    l_tags = [tuple[1] for tuple in l_tuple_tags]
    gmsh.finalize()

    return l_tags


def global_mesh(global_mesh_file, l_regions=[]):
    """
    Create a Mesh object that represents the initial global mesh
    @param global_mesh_file: path to the initial mesh file
    @param l_regions: list that contains tags of materials of the initial mesh
    @return: return the Mesh object
    """
    global_mesh = Mesh(global_mesh_file, l_regions)
    return global_mesh


def find_local_vertices(global_mesh, local_vertices, trgl_id):
    """
    Add vertices of a triangle in the list of vertices to keep in the new mesh file
    @param global_mesh: Mesh object that represents the initial global mesh
    @param local_vertices: list of vertex indexes to keep in the new mesh (indexes from the initial mesh)
    @param trgl_id: triangle index
    """
    triangle = global_mesh.l_triangles_[trgl_id - 1]
    for vertex_id in triangle.l_trgl_vertices_idx_:
        if vertex_id not in local_vertices:
            local_vertices.append(vertex_id)


def find_local_edges(global_mesh, local_edges, trgl_id):
    """
    Add edges of a triangle in the list of edges to keep in the new mesh file
    @param global_mesh: Mesh object that represents the initial global mesh
    @param local_edges: list of edge indexes to keep in the new mesh (indexes from the initial mesh)
    @param trgl_id: triangle index
    """
    l_couple_vrtx = [(1, 2), (2, 0), (0, 1)]
    triangle = global_mesh.l_triangles_[trgl_id - 1]

    for tuple in l_couple_vrtx:
        endpoints = triangle.l_trgl_vertices_idx_[tuple[0]], triangle.l_trgl_vertices_idx_[tuple[1]]
        endpoints = sorted(endpoints)
        # if endpoints not in local_edges:
        local_edges.append(endpoints)


def tag_reg_boundary_edges(local_edges, l_trgl_id, global_mesh):
    """
    Clean the list of edge indexes to keep in the new mesh and return the list of edge tags in the same order
    @param local_edges: list of edge indexes to keep in the new mesh (indexes from the initial mesh)
    @param l_trgl_id: list of triangle indexes of the new mesh (indexes from the initial mesh)
    @param global_mesh: Mesh object that represents the initial global mesh
    @return: l_tag: list of edge tag in the same order than l_clean_local_edges
             l_clean_local_edges: list of edge indexes to keep in the new mesh without duplicates
    """
    l_tag_global_mesh = global_mesh.tag_region_boundary_edges_
    l_tag = []
    l_clean_local_edges = []
    for edge in range(len(local_edges)):
        if local_edges[edge] not in l_clean_local_edges:
            l_clean_local_edges.append(local_edges[edge])
            triangle_id_l = edge // 3
            triangle = l_trgl_id[triangle_id_l]
            edge_global_id = (triangle - 1) * 3 + edge % 3

            l_tag.append(l_tag_global_mesh[edge_global_id])
    return l_tag, l_clean_local_edges


def find_local_vertices_and_edges(global_mesh, l_trgl_id):
    """
    Create the lists of indexes of vertices and edges to keep in the new mesh
    @param global_mesh: Mesh object that represents the initial global mesh
    @param l_trgl_id: list of triangle indexes of the new mesh (indexes from the initial mesh)
    @return: local_vertices: list of vertex indexes to keep in the new mesh (indexes from the initial mesh)
             l_clean_local_edges: list of edge indexes to keep in the new mesh without duplicates
             l_tag: list of edge tag in the same order than l_clean_local_edges
    """
    local_vertices = []
    local_edges = []

    for trgl_id in l_trgl_id:
        find_local_vertices(global_mesh, local_vertices, trgl_id)
        find_local_edges(global_mesh, local_edges, trgl_id)

    l_tag, l_clean_local_edges = tag_reg_boundary_edges(local_edges, l_trgl_id, global_mesh)
    return local_vertices, l_clean_local_edges, l_tag


def write_header(file_initial_mesh, file_new_mesh, local_l_vertices):
    """
    Write the header of the 2d .mesh file
    @param file_initial_mesh: initial mesh file opened in reading mode
    @param file_new_mesh: new mesh file opened in writing mode
    @param local_l_vertices: list of vertex indexes to keep in the new mesh (indexes from the initial mesh)
    @return: vertices_number: number of vertices in the new mesh
    """
    line = file_initial_mesh.readline()

    while "Dimension" not in line:
        file_new_mesh.write(line)
        line = file_initial_mesh.readline()

    file_new_mesh.write("Dimension")
    file_new_mesh.write(" 2\n\n")  # impose 2D
    line = file_initial_mesh.readline()

    ## "Vertices"
    while "Vertices" not in line:
        line = file_initial_mesh.readline()
        file_new_mesh.write(line)

    ## Number of vertices
    line = file_initial_mesh.readline()
    vertices_number = len(local_l_vertices)
    file_new_mesh.write(str(vertices_number) + " \n")

    return vertices_number


def write_vertices(local_vertices, file_init_mesh, file_new_mesh):
    """
    Write vertices in the new mesh file and return the list of sorted vertices (by indexes)
    @param local_vertices: list of vertex indexes (index from the .mesh file)
    @param file_init_mesh: initial mesh file opened in reading mode
    @param file_new_mesh: new mesh file opened in writing mode
    @return: list of vertex sorted by their indexes
    """
    sort_local_vertices = sorted(local_vertices)
    vertex_id_read = 0

    for vertex_id in sort_local_vertices:
        id_diff = vertex_id - vertex_id_read
        for j in range(id_diff):
            line = file_init_mesh.readline()
        split_line = line.split()

        if len(split_line) > 3:
            line = " ".join(split_line[:3]) + "\n"
        else:
            line = " ".join(split_line) + "\n"
        file_new_mesh.write(line)

        vertex_id_read = vertex_id

    return sort_local_vertices


def write_edges(local_edges, sort_local_vertices, l_edge_tag, file_new_mesh):
    """
    Write edges in the new mesh file
    @param local_edges: list of edge indexes to keep in the new mesh
    @param sort_local_vertices: list of sorted vertex indexes (in ascending order) to keep in the new mesh
    @param l_edge_tag: list of edge tag in the same order than local_edges
    @param file_new_mesh: new mesh file opened in writing mode
    """
    file_new_mesh.write("Edges\n")
    file_new_mesh.write(str(len(local_edges)) + "\n")

    id = 0
    for tuple in local_edges:
        first_end_point_id = sort_local_vertices.index(tuple[0]) + 1
        second_end_point_id = sort_local_vertices.index(tuple[1]) + 1
        edge_tag = l_edge_tag[id]
        file_new_mesh.write(
            str(first_end_point_id) + " " + str(second_end_point_id) + " " + str(edge_tag) + "\n")
        id += 1


def write_triangles(l_triangles, sort_local_vertices, file_initial_mesh, file_new_mesh):
    """
    Write triangles in the new mesh file
    @param l_triangles: list of triangle indexes of the new mesh
    @param sort_local_vertices: list of sorted vertex indexes (in ascending order) to keep in the new mesh
    @param file_initial_mesh: initial mesh file opened in reading mode
    @param file_new_mesh: new mesh file opened in writing mode
    """
    file_initial_mesh.readline()
    file_new_mesh.write("Triangles\n")
    file_new_mesh.write(str(len(l_triangles)) + "\n")

    sort_local_triangles = sorted(l_triangles)
    triangle_id_read = 0

    for triangle_id in sort_local_triangles:
        id_diff = triangle_id - triangle_id_read
        for j in range(id_diff):
            line = file_initial_mesh.readline()

        split_line = line.split(" ")

        for vrtx_id in split_line:
            if vrtx_id == "" or "\n" in vrtx_id:
                pass
            else:
                split_line[split_line.index(vrtx_id)] = str(sort_local_vertices.index(int(vrtx_id)) + 1)

        line = " ".join(split_line)
        file_new_mesh.write(line)
        triangle_id_read = triangle_id


def write_local_mesh(filepath_initial_mesh, filepath_new_mesh, l_triangles, l_triangles_gmsh, l_regions):
    """
    Write the mesh file of the new mesh (from a list of triangle of a initial mesh)
    @param filepath_initial_mesh: path to initial mesh file
    @param filepath_new_mesh: path to new mesh file
    @param l_triangles: list of triangle indexes of the new mesh
    @param l_triangles_gmsh: list of triangle indexes of the new mesh (gmsh indexes)
    @param l_regions: list that contains tags of materials of the initial mesh
    """
    file_init_mesh = open(filepath_initial_mesh, "r")
    file_new_mesh = open(filepath_new_mesh, "w")

    mesh = global_mesh(filepath_initial_mesh, l_regions)
    if not l_triangles:
        l_triangles = [i for i in range(1, len(mesh.l_triangles_) + 1)]

    local_vertices, local_edges, l_tag = find_local_vertices_and_edges(mesh, l_triangles_gmsh)

    write_header(file_init_mesh, file_new_mesh, local_vertices)

    # Vertices
    sort_local_vertices = write_vertices(local_vertices, file_init_mesh, file_new_mesh)

    # Corners
    # RequiredVertices
    # RequiredEdges

    # Edges
    line = file_init_mesh.readline()
    while "Edges" not in line:
        line = file_init_mesh.readline()
    write_edges(local_edges, sort_local_vertices, l_tag, file_new_mesh)

    # Triangles
    line = file_init_mesh.readline()
    while "Triangles" not in line:
        line = file_init_mesh.readline()
    write_triangles(l_triangles, sort_local_vertices, file_init_mesh, file_new_mesh)

    file_new_mesh.write("End")
    file_init_mesh.close()
    file_new_mesh.close()


def build_l_triangles(global_mesh_filepath, l_region_id):
    """
    Create lists of all triangles that are in the inner region and in outer region
    @param global_mesh_filepath: Path to initial mesh file
    @param l_region_id: list that contains tags of materials in inner and outer regions
    @return: l_inner_triangles: list of triangle indexes of the inner region
             l_inner_triangles_gmsh_order: list of triangle indexes of the inner region with the gmsh indexes
             l_triangles_other: list of triangle indexes of the outer region
             l_triangles_other_gmsh_order: list of triangle indexes of the outer region with the gmsh indexes
    """
    open_gmsh_file(global_mesh_filepath)
    array_dim_global, array_trgl_global_idx, array_vertex_global_idx = gmsh.model.mesh.getElements(dim=2, tag=-1)
    l_sort_trgl_global_idx = sorted(list(array_trgl_global_idx[0]))
    l_trgl_global_idx = list(array_trgl_global_idx[0])
    l_inner_triangles = []
    l_inner_triangles_gmsh_order = []
    l_outer_triangles = l_sort_trgl_global_idx.copy()
    l_outer_triangles_gmsh = l_trgl_global_idx.copy()

    for region_id in l_region_id:
        array_dim_region, array_trgl_region_idx, array_vertex_region_idx = gmsh.model.mesh.getElements(dim=2,
                                                                                                       tag=region_id)

        l_trgl_region_idx = list(array_trgl_region_idx[0])

        for trgl_id in l_trgl_region_idx:
            l_outer_triangles.remove(trgl_id)
            l_outer_triangles_gmsh.remove(trgl_id)

        l_triangles_reg = [l_sort_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_trgl_region_idx]
        l_inner_triangles += l_triangles_reg

        l_triangles_reg_gmsh = [l_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_trgl_region_idx]
        l_inner_triangles_gmsh_order += l_triangles_reg_gmsh

    l_triangles_other = [l_sort_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_outer_triangles]
    l_triangles_other_gmsh_order = [l_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_outer_triangles_gmsh]

    gmsh.finalize()

    return l_inner_triangles, l_inner_triangles_gmsh_order, l_triangles_other, l_triangles_other_gmsh_order


def build_2_local_meshes(filepath_initial_mesh, filepath_inner_region, filepath_outer_region, l_regions):
    """
    Create two sub-meshes that represent the inner (modified region) and outer region
    @param filepath_initial_mesh: path to initial mesh file
    @param filepath_inner_region: Path to mesh file of the inner region
    @param filepath_outer_region: Path to mesh file of the outer region
    @param l_regions: list that contains tags of materials in inner and outer regions
    """
    l_inner_region = l_regions[0]
    l_inner_trgl, l_triangles_gmsh, l_outer_trgl, l_outer_triangles_gmsh_order = build_l_triangles(
        filepath_initial_mesh, l_inner_region)
    write_local_mesh(filepath_initial_mesh, filepath_inner_region, l_inner_trgl, l_triangles_gmsh, l_regions)
    write_local_mesh(filepath_initial_mesh, filepath_outer_region, l_outer_trgl, l_outer_triangles_gmsh_order,
                     l_regions)


def merge(mesh_reg1, mesh_reg2, output_file):
    """
    Merge two sub-parts (iner and outer parts) of the mesh
    @param mesh_reg1:first mesh sub-part file path
    @param mesh_reg2:second mesh sub-part file path
    @param output_file: path to global mesh that contains the two sub-parts
    """
    open_gmsh_file(mesh_reg1)
    gmsh.merge(mesh_reg2)

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.write(output_file)


if __name__ == '__main__':
    filepath_to_read = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/ext_line/region/curve1_cut.o.mesh"
    filepath_to_write_cut = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/ext_line/region/mesh_regcut.mesh"
    filepath_to_write_other = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/ext_line/region/mesh_regother.mesh"
    filepath_final = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/ext_line/region/mesh_merge.mesh"
    l_regions = [[1, 3], [4]]
    build_2_local_meshes(filepath_to_read, filepath_to_write_cut, filepath_to_write_other, l_regions)

    merge("C:/Users/legentil1.UL/Downloads/test.msh", "C:/Users/legentil1.UL/Downloads/test2.msh",
          "C:/Users/legentil1.UL/Downloads/test_merge.msh")
