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
from lib_data_structure.Mesh import Mesh
import gmsh


def open_gmsh_file(file_path):
    """
    Open mesh file with gmsh
    :param file_path: absolute path
    """
    gmsh.initialize([])
    gmsh.open(file_path)


def get_triangle_tags(global_mesh_filepath):
    gmsh.initialize([])
    gmsh.open(global_mesh_filepath)
    l_tuple_tags = gmsh.model.getEntities(dim=2)
    l_tags = [tuple[1] for tuple in l_tuple_tags]
    gmsh.finalize()

    return l_tags


def write_header(file_to_read, file_to_write, local_l_vertices):
    line = file_to_read.readline()

    while "Dimension" not in line:
        file_to_write.write(line)
        line = file_to_read.readline()

    file_to_write.write("Dimension")
    file_to_write.write(" 2\n\n")  # impose 2D
    line = file_to_read.readline()

    ## "Vertices"
    while "Vertices" not in line:
        line = file_to_read.readline()
        file_to_write.write(line)

    ## Number of vertices
    line = file_to_read.readline()
    vertices_number = len(local_l_vertices)
    file_to_write.write(str(vertices_number) + " \n")

    return vertices_number


def global_mesh(global_mesh_file, l_regions = []):
    global_mesh = Mesh(global_mesh_file, l_regions)
    return global_mesh


def find_local_vertices_and_edges(global_mesh, l_trgl_id):
    local_vertices = []
    local_edges = []

    for trgl_id in l_trgl_id:
        find_local_vertices(global_mesh, local_vertices, trgl_id)
        find_local_edges(global_mesh, local_edges, trgl_id)

    l_tag, l_clean_local_edges = tag_reg_boundary_edges(local_edges, l_trgl_id, global_mesh)
    return local_vertices, l_clean_local_edges, l_tag


def find_local_vertices(global_mesh, local_vertices, trgl_id):
    triangle = global_mesh.l_triangles_[trgl_id - 1]
    for vertex_id in triangle.l_trgl_vertices_idx_:
        if vertex_id not in local_vertices:
            local_vertices.append(vertex_id)


def find_local_edges(global_mesh, local_edges, trgl_id):
    l_couple_vrtx = [(1, 2), (2, 0), (0, 1)]
    triangle = global_mesh.l_triangles_[trgl_id - 1]

    for tuple in l_couple_vrtx:
        endpoints = triangle.l_trgl_vertices_idx_[tuple[0]], triangle.l_trgl_vertices_idx_[tuple[1]]
        endpoints = sorted(endpoints)
        #if endpoints not in local_edges:
        local_edges.append(endpoints)

def tag_reg_boundary_edges(local_edges, l_trgl_id, global_mesh):
    l_tag_global_mesh = global_mesh.tag_region_boundary_edges_
    #print(local_edges)
    l_tag = []
    l_clean_local_edges = []
    for edge in range(len(local_edges)):
        if local_edges[edge] not in l_clean_local_edges:
            l_clean_local_edges.append(local_edges[edge])
            triangle_id_l = edge//3
            triangle = l_trgl_id[triangle_id_l]
            edge_global_id = (triangle-1)*3 + edge % 3

            l_tag.append(l_tag_global_mesh[edge_global_id])
    #print(l_clean_local_edges)
    return l_tag, l_clean_local_edges



def write_vertices(local_vertices, file_to_read, file_to_write):
    sort_local_vertices = sorted(local_vertices)
    vertex_id_read = 0

    for vertex_id in sort_local_vertices:
        id_diff = vertex_id - vertex_id_read
        for j in range(id_diff):
            line = file_to_read.readline()
        split_line = line.split()
        #print(line)
        if len(split_line) > 3:
            line = " ".join(split_line[:3]) + "\n"
        else:
            line = " ".join(split_line) + "\n"
        file_to_write.write(line)

        vertex_id_read = vertex_id

    return sort_local_vertices


def write_edges(local_edges, sort_local_vertices, l_edge_tag, file_to_write):
    file_to_write.write("Edges\n")
    file_to_write.write(str(len(local_edges)) + "\n")

    id = 0
    for tuple in local_edges:
        first_end_point_id = sort_local_vertices.index(tuple[0]) + 1
        second_end_point_id = sort_local_vertices.index(tuple[1]) + 1
        edge_tag = l_edge_tag[id]
        file_to_write.write(str(first_end_point_id) + " " + str(second_end_point_id) + " " + str(edge_tag) + "\n")
        id+=1


def write_triangles(l_triangles, sort_local_vertices, file_to_read, file_to_write):
    file_to_read.readline()
    file_to_write.write("Triangles\n")
    file_to_write.write(str(len(l_triangles)) + "\n")

    sort_local_triangles = sorted(l_triangles)
    triangle_id_read = 0

    for triangle_id in sort_local_triangles:
        id_diff = triangle_id - triangle_id_read
        for j in range(id_diff):
            line = file_to_read.readline()

        split_line = line.split(" ")  ##### A modifier (il n'y peut être pas toujours un espace en début de ligne)

        for vrtx_id in split_line:
            if vrtx_id == "" or "\n" in vrtx_id:
                pass
            else:
                split_line[split_line.index(vrtx_id)] = str(sort_local_vertices.index(int(vrtx_id)) + 1)
        # print(line)
        line = " ".join(split_line)
        file_to_write.write(line)
        triangle_id_read = triangle_id

    # print(file_to_read.readline())


def write_local_mesh(filepath_to_read, filepath_to_write, l_triangles, l_triangles_gmsh, l_regions):
    file_to_read = open(filepath_to_read, "r")
    file_to_write = open(filepath_to_write, "w")

    mesh = global_mesh(filepath_to_read, l_regions)
    if not l_triangles:
        l_triangles = [i for i in range(1, len(mesh.l_triangles_) + 1)]

    local_vertices, local_edges, l_tag = find_local_vertices_and_edges(mesh, l_triangles_gmsh)

    write_header(file_to_read, file_to_write, local_vertices)

    # Vertices
    sort_local_vertices = write_vertices(local_vertices, file_to_read, file_to_write)

    # Corners
    # RequiredVertices
    # RequiredEdges

    # Edges
    line = file_to_read.readline()
    while "Edges" not in line:
        line = file_to_read.readline()
    write_edges(local_edges, sort_local_vertices, l_tag, file_to_write)
    #print(len(l_tag))

    # Triangles
    line = file_to_read.readline()
    while "Triangles" not in line:
        line = file_to_read.readline()
    write_triangles(l_triangles, sort_local_vertices, file_to_read, file_to_write)

    file_to_write.write("End")
    file_to_read.close()
    file_to_write.close()


def build_l_triangles(global_mesh_filepath, l_region_id):

    open_gmsh_file(global_mesh_filepath)
    array_dim_global, array_trgl_global_idx, array_vertex_global_idx = gmsh.model.mesh.getElements(dim=2, tag=-1)
    l_sort_trgl_global_idx = sorted(list(array_trgl_global_idx[0]))
    l_trgl_global_idx = list(array_trgl_global_idx[0])
    l_triangles_cut = []
    l_triangles_cut_gmsh_order = []
    l_triangles_rest = l_sort_trgl_global_idx.copy()
    l_triangles_rest_gmsh = l_trgl_global_idx.copy()

    for region_id in l_region_id:
        array_dim_region, array_trgl_region_idx, array_vertex_region_idx = gmsh.model.mesh.getElements(dim=2,
                                                                                                   tag=region_id)

        l_trgl_region_idx = list(array_trgl_region_idx[0])

        for trgl_id in l_trgl_region_idx:
            l_triangles_rest.remove(trgl_id)
            l_triangles_rest_gmsh.remove(trgl_id)


        l_triangles_reg = [l_sort_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_trgl_region_idx]
        l_triangles_cut += l_triangles_reg

        l_triangles_reg_gmsh = [l_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_trgl_region_idx]
        l_triangles_cut_gmsh_order += l_triangles_reg_gmsh

    l_triangles_other = [l_sort_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_triangles_rest]
    l_triangles_other_gmsh_order = [l_trgl_global_idx.index(reg_trgl) + 1 for reg_trgl in l_triangles_rest_gmsh]

    gmsh.finalize()

    return l_triangles_cut, l_triangles_cut_gmsh_order, l_triangles_other, l_triangles_other_gmsh_order

def build_2_local_meshes(filepath_to_read, filepath_to_write_cut_region, filepath_to_write_other, l_regions):
    l_cut_region = l_regions[0]
    l_trgl_cut, l_triangles_gmsh, l_trgl_other, l_triangles_other_gmsh_order = build_l_triangles(filepath_to_read, l_cut_region)
    write_local_mesh(filepath_to_read, filepath_to_write_cut_region, l_trgl_cut, l_triangles_gmsh, l_regions)
    write_local_mesh(filepath_to_read, filepath_to_write_other, l_trgl_other, l_triangles_other_gmsh_order, l_regions)

def merge(mesh_reg1, mesh_reg2, output_file):
    open_gmsh_file(mesh_reg1)
    gmsh.merge(mesh_reg2)

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.write(output_file)
    # gmsh.open("D:\Simple_test\simple_mesh_rewrite_merge2.mesh")t suivre l
    # print(gmsh.model.getDimension())


if __name__ == '__main__':
    filepath_to_read = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/surf1_init_ref.o.mesh"
    filepath_to_write_cut = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_regcut.mesh"
    filepath_to_write_other = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_regother.mesh"
    filepath_final = "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_merge.mesh"

    # local_1_vertices = [3, 4, 7, 8, 9, 17, 21, 22, 24]
    # local_2_vertices = [1, 2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29]

    region1_trgl = [34, 35, 36, 37, 38, 39, 40, 41]
    region2_trgl = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                    28, 29, 30, 31, 32, 33]

    #my_mesh = global_mesh("D:\Simple_test\simple_mesh.mesh")

    #write_local_mesh(filepath_to_read, filepath_to_write, region1_trgl)

    #l_trgl_cut, l_triangles_gmsh, l_trgl_other, l_triangles_other_gmsh_order = build_l_triangles(filepath_to_read, [1,3])
    #print(l_trgl_cut)
    #write_local_mesh(filepath_to_read, filepath_to_write, l_trgl_cut, l_triangles_gmsh)
    #build_2_local_meshes(filepath_to_read, filepath_to_write_cut, filepath_to_write_other, [[2, 6], [1, 3, 4]])

    open_gmsh_file("/home/legentil/Bureau/raffinement/mesh_800/mesh_merge1234.mesh")

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.write("/home/legentil/Bureau/raffinement/mesh_800/mesh_merge1234_o.mesh")

    #merge("/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_regcut_req.o.mesh", "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_regother.mesh", filepath_final)
