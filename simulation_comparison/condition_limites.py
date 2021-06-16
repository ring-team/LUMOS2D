"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """


def open_read_file(filepath):
    """
    Open file in reading mode
    @param filepath: file path to read
    @return: readable file
    """
    my_file = open(filepath, "r")
    return my_file


def open_write_file(filepath):
    """
    Open file in writing mode
    @param filepath: file path to write
    @return: file in writing mode
    """
    my_file = open(filepath, "w")
    return my_file


def write_boundary_condition(file_read, file_write, xmin, xmax, ymin, ymax):
    """
    Change the mesh file to integrate boundary condition for seismic simulation in a rectangle 2D section defined by 4 corners:
    the corners hve coordinates xmin, xmax, ymin, ymax
    (for free surface and 3 for absorbing condition)
    @param file_read: mesh file path to modify
    @param file_write: modified mesh file path that contains the boundary conditions
    @param xmin: minimal x coordinate
    @param xmax: maximal x coordinate
    @param ymin: minimal y coordinate
    @param ymax: minimal y coordinate
    """
    file_to_read = open_read_file(file_read)
    file_to_write = open_write_file(file_write)
    line = file_to_read.readline()
    file_to_write.write(line)
    line = file_to_read.readline()
    file_to_write.write(line)
    line = file_to_read.readline()
    while line:
        split_line = line.split()
        if (split_line[1] == xmin and split_line[2] == ymin) or (split_line[1] == xmin and split_line[2] == ymax) \
                or (split_line[1] == xmax and split_line[2] == ymin) or (
                split_line[1] == xmax and split_line[2] == ymax):
            split_line[3] = "0"
            print("ok")

        elif split_line[1] == xmin:
            split_line[3] = "3"

        elif split_line[1] == xmax:
            split_line[3] = "3"

        elif split_line[2] == ymax:
            split_line[3] = "1"

        elif split_line[2] == ymin:
            split_line[3] = "3"

        line = " ".join(split_line) + "\n"
        file_to_write.write(line)
        line = file_to_read.readline()


if __name__ == '__main__':
    write_boundary_condition("/home/legentil/Bureau/raffinement/ls700/mesh_refined.1_cl.node",
                             "/home/legentil/Bureau/raffinement/ls700/mesh_refined.1.node", "-5291.109375000000000",
                             "10949.264648438000222", "-1837.562988281199978", "3247.130371093699978")
