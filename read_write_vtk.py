"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
"""@file read_write_vtk.py scripts for Paraview visualisation 
*Add scalar field to a .vtk file """

def build_l_sol(sol_file_path):
    """
    Create the list of scalar field values
    @param sol_file_path: File where the scalar field is saved (.sol)
    @return: list of scalar field values
    """
    sol_file = open(sol_file_path, "r")
    for i in range(9):
        sol_file.readline()
    l_sol = []
    line = sol_file.readline()
    while "End" not in line:
        l_sol.append(float(line))
        line = sol_file.readline()

    return l_sol

def write_vtk_begin(filepath_to_read, filepath_to_write):
    """
    Write all mandatory data in the new .vtk file
    @param filepath_to_read: mesh file path without the scalar field property (.vtk)
    @param filepath_to_write: final mesh file path with the scalar field property (.vtk)
    @return: A file .vtk with all mandatory data (without the scalar field values)
    """
    file_to_read = open(filepath_to_read, "r")
    file_to_write = open(filepath_to_write, "w")

    line = file_to_read.readline()

    while "POINT_DATA" not in line:
        file_to_write.write(line)
        line = file_to_read.readline()

    file_to_write.write(line)
    line = file_to_read.readline()
    split_line = line.split()
    split_line[-1] = str(int(split_line[-1]) + 1) + "\n"
    line = " ".join(split_line)
    file_to_write.write(line)

    line = file_to_read.readline()
    while line:
        file_to_write.write(line)
        line = file_to_read.readline()

    return file_to_write

def add_solution(l_sol, file_to_write, sol_name, dim):
    """
    Add scalar field property to the .vtk file
    @param l_sol: list of scalar field values
    @param file_to_write: final mesh file path with the scalar field property (.vtk)
    @param sol_name: property name (string)
    @param dim: scalar field dimension (1)
    """
    file_to_write.write(sol_name + " " + str(dim) + " " + str(len(l_sol)) + " double\n")

    for i in l_sol:
        file_to_write.write(str(i) + "\n")

def write_vtk_with_sol(filepath_to_read, filepath_to_write, l_sol):
    """
    Write the scalar field property in .vtk file
    @param filepath_to_read:  mesh file path without the scalar field property (.vtk)
    @param filepath_to_write: final mesh file path with the scalar field property (.vtk)
    @param l_sol: list of scalar field values
    """
    file_to_write = write_vtk_begin(filepath_to_read, filepath_to_write)
    add_solution(l_sol, file_to_write, "Distance", 1)

def write_vtk_sol(filepath_vtk_to_read, filepath_sol_to_read, filepath_vtk_to_write):
    """
    Read the scalar field file (.sol) and write the scalar field property in .vtk file
    @param filepath_vtk_to_read: mesh file path without the scalar field property (.vtk)
    @param filepath_sol_to_read: file where the scalar field is saved (.sol)
    @param filepath_vtk_to_write: final mesh file path with the scalar field property (.vtk)
    """
    l_sol = build_l_sol(filepath_sol_to_read)
    write_vtk_with_sol(filepath_vtk_to_read, filepath_vtk_to_write, l_sol)
    print("Writing VTK file -- DONE")

if __name__ == '__main__':
    file_to_read = "/home/legentil/Test_jupyter/surf1_init_ref.o.vtk"
    file_to_write = "/home/legentil/Test_jupyter/surf1_init_ref_sol.o.vtk"

    l_sol = build_l_sol("/home/legentil/Bureau/Test_jupyter/surf1_init_ref.o.sol")
    print(l_sol)
    write_vtk_sol(file_to_read, "/home/legentil/Bureau/Test_jupyter/surf1_init_ref.o.sol", file_to_write)
