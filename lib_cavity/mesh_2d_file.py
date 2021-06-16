"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """


def RepresentsInt(number):
    """
    Test if a number represents an number
    @param number: number to test (string)
    @return: True if it is a number, False otherwise
    """
    try:
        int(number)
        return True
    except ValueError:
        return False


def write_header(file_initial_mesh, file_new_mesh):
    """
    Write header of a 2D .mesh file and return the number of vertices of the initial mesh
    @param file_initial_mesh: initial mesh file opened in reading mode
    @param file_new_mesh: new mesh file opened in writing mode
    @return: number of vertices
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
    vertices_number = int(line)
    file_new_mesh.write(line)

    return vertices_number


def delete_vertices_tag(file_init_mesh, file_new_mesh, vertices_number):
    """
    Write mesh vertices from gmsh output .mesh file. There is a file format problem. This code delete a vertice tag from
     gmsh format.
    @param file_init_mesh: initial mesh file opened in reading mode (gmsh output)
    @param file_new_mesh: new mesh file opened in writing mode
    @param vertices_number: number of vertices
    """
    for vertex_id in range(0, vertices_number):
        line = file_init_mesh.readline()
        split_line = line.split()

        split_line[-1] = "\n"

        line = " ".join(split_line)
        file_new_mesh.write(line)


def write_corners(file_init_mesh, file_new_mesh):
    """
    Write lines about corners
    @param file_init_mesh: initial mesh file opened in reading mode (gmsh output)
    @param file_new_mesh: new mesh file opened in writing mode
    """
    line = file_init_mesh.readline()
    while "Edges" not in line:
        file_new_mesh.write(line)
        line = file_init_mesh.readline()

    file_new_mesh.write(line)

def read_required_edges(file_init_mesh):
    """
    Read RequiredEdges in initial mesh file
    @param file_init_mesh:  initial mesh file opened in reading mode
    @return: line: last read line
             l_id_req_edge: list of required edges
    """
    l_id_req_edge=[]
    line= file_init_mesh.readline()
    number_required = int(line)
    for i in range(number_required):
        line = file_init_mesh.readline()
        l_id_req_edge.append(int(line))
    line = file_init_mesh.readline()
    return line, l_id_req_edge

def remove_duplicates(req_before, req_after):
    """
    Concatenate two list of requiredEdges by removing duplicates
    @param req_before: first list of RequiredEdges (in the initial mesh file)
    @param req_after: second list of RequiredEdges (in the future mesh file)
    @return: the merged list without duplicates
    """
    for i in req_after:
        if i not in req_before:
            req_before.append(i)
    return req_before

def required_edges(file_init_mesh, file_new_mesh):
    """
    Write RequiredEdges
    @param file_init_mesh: initial mesh file opened in reading mode
    @param file_new_mesh: new mesh file opened in writing mode
    @return: list of requiredEdges
    """
    line = file_init_mesh.readline()
    file_new_mesh.write(line)
    line = file_init_mesh.readline()

    edge_id = 1
    l_req_edges = []
    l_id_req_edge = []

    while "Triangles" not in line:

        split_line = line.split(" ")

        if RepresentsInt(split_line[-1]) and int(split_line[-1]) >= 100:
            l_req_edges.append(edge_id)
            print(line)
            split_line[-1] = str(int(split_line[-1]) - 100) + "\n"

        line = " ".join(split_line)

        file_new_mesh.write(line)
        line = file_init_mesh.readline()
        edge_id += 1

        if "RequiredEdges" in line:
            line, l_id_req_edge = read_required_edges(file_init_mesh)

    l_required_edge = remove_duplicates(l_id_req_edge, l_req_edges)

    file_new_mesh.write("RequiredEdges\n")
    file_new_mesh.write(str(len(l_required_edge)) + "\n")

    for i in l_required_edge:
        file_new_mesh.write(str(i) + "\n")

    file_new_mesh.write("\nTriangles\n")

    return l_required_edge

def write_end_file(file_init_mesh, file_new_file):
    """
    Write end of the .mesh file
    @param file_init_mesh: initial mesh file opened in reading mode
    @param file_new_file: new mesh file opened in writing mode
    """
    line = file_init_mesh.readline()

    while "End" not in line:
        file_new_file.write(line)
        line = file_init_mesh.readline()

    file_new_file.write("End")


def write_2D_mesh_req_file(init_mesh_file_path, new_mesh_file_path):
    """
    Write 2d .mesh file with required edges
    @param init_mesh_file_path: initial mesh file to open in reading mode
    @param new_mesh_file_path: new mesh file to open in writing mode
    """
    file_init_mesh = open(init_mesh_file_path, "r")
    file_new_file = open(new_mesh_file_path, "w")

    vertices_number = write_header(file_init_mesh, file_new_file)

    write_corners(file_init_mesh, file_new_file)
    required_edges(file_init_mesh, file_new_file)

    write_end_file(file_init_mesh, file_new_file)

    file_init_mesh.close()
    file_new_file.close()


def write_2D_mesh_file_gmsh(init_mesh_file_path, new_mesh_file_path):
    """
    Write 2d .mesh file and convert the output of gmsh
    @param init_mesh_file_path: initial mesh file to open in reading mode (gmsh output)
    @param new_mesh_file_path: new mesh file to open in writing mode
    """
    file_init_mesh = open(init_mesh_file_path, "r")
    file_new_file = open(new_mesh_file_path, "w")

    vertices_number = write_header(file_init_mesh, file_new_file)

    delete_vertices_tag(file_init_mesh, file_new_file, vertices_number)
    write_corners(file_init_mesh, file_new_file)
    required_edges(file_init_mesh, file_new_file)

    write_end_file(file_init_mesh, file_new_file)

    file_init_mesh.close()
    file_new_file.close()




if __name__ == '__main__':
    # write_2D_mesh_file_gmsh("/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_merge.mesh",
    # "/home/legentil/Documents/Data/_2D/Simple_Models/A1/surf1/Code_Test/water_oil/refine/couche_2/mesh_merge_2d.mesh")

    # write_2D_mesh_file(
    #     "D:/Programmation/ringlab/LUMOS/python/data/reg_cut_req_test.mesh",
    #     "D:/Programmation/ringlab/LUMOS/python/data/reg_cut_req_test2.mesh")

    write_2D_mesh_file_gmsh("D:/Programmation/ringlab/LUMOS/python/data/mesh_merge.mesh",
                            "D:/Programmation/ringlab/LUMOS/python/data/mesh_merge_test.mesh")
