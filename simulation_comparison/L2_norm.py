"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np


def open_sismo(path, receiver_number):
    """
    Open seismogram file at one receiver
    @param path: entire file path of the folder containing all seismograms
    @param receiver_number: receiver number (int)
    @return: opened readable file
    """
    number_4d = str("%04d" % receiver_number)
    filepath = path + "/sismo_U." + number_4d + ".dat"
    my_file = open(filepath, "r")
    return my_file


def open_write_file(filepath):
    """
    Open file to write information
    @param filepath: entire file path
    @return: file to write
    """
    my_file = open(filepath, "w")
    return my_file


def norm_computation(sismo_file):
    """
    Compute the displacement norm at each time step
    @param sismo_file: opended readable file
    @return: list of the displacement norm and list of u_x (displacement along x axis)
    """
    l_u_norm = []
    l_u_x = []
    line = sismo_file.readline()
    while line:
        split_line = line.split()
        u_x, u_z = float(split_line[0]), float(split_line[1])
        l_u_x.append(u_x)
        l_u_norm.append(sqrt(u_x ** 2 + u_z ** 2))
        line = sismo_file.readline()
    return l_u_norm, l_u_x


def approximation_trapeze(l_function, delta_x):
    """
    Integration trapezoidal approximation
    @param l_function: function values
    @param delta_x: x step
    @return: Integration approximation
    """
    val_integrale = 0
    for depth in range(0, len(l_function) - 1):
        val_integrale += 0.5 * delta_x * (l_function[depth] + l_function[depth + 1])
    return val_integrale


def norm_l2(sismo_ref, sismo, nb_rec, time_step):
    """
    Compute the L2-norm between two seismograms
    @param sismo_ref: reference seismogram file path (folder containing the seismograms)
    @param sismo: seismogram file path to compare (folder containing the seismograms)
    @param nb_rec: number of receivers (=number of seismogram files)
    @param time_step: time step
    @return: l_L2_norm: list of L2-norm value at each receiver,
             max_id: receiver index where the L2-norm is maximum,
             l_inte_ref: list of Integration value of the displacement norm of the reference at each receiver.
    """
    l_L2_norm = []
    l_ref_inte = []
    max = 0
    max_id = 0
    for receiver_id in range(1, nb_rec + 1):
        my_file_ref = open_sismo(sismo_ref, receiver_id)
        my_file = open_sismo(sismo, receiver_id)

        l_u_norm_ref, l_u_x_ref = norm_computation(my_file_ref)
        l_u_norm_ref_square = [l_u_norm_ref[i] * l_u_norm_ref[i] for i in range(0, len(l_u_norm_ref))]
        l_u_norm, l_u_x = norm_computation(my_file)
        l_u_diff = [(l_u_norm_ref[i] - l_u_norm[i]) * (l_u_norm_ref[i] - l_u_norm[i]) for i in range(0, len(l_u_norm))]

        val_integrale_ref = approximation_trapeze(l_u_norm_ref, time_step)
        val_integrale_ref_square = approximation_trapeze(l_u_norm_ref_square, time_step)
        val_diff_integrale = approximation_trapeze(l_u_diff, time_step)

        val = 100 * sqrt(val_diff_integrale / val_integrale_ref_square)
        l_L2_norm.append(val)
        l_ref_inte.append(val_integrale_ref)
        if val > max:
            max = val
            max_id = receiver_id

    print(max)
    print(np.mean(l_L2_norm))
    return l_L2_norm, max_id, l_ref_inte


if __name__ == '__main__':
    sismo_ref = "/home/legentil/Bureau/raffinement/ls700/FILM/FILM"
    sismo = "/home/legentil/Bureau/raffinement/ls500/FILM/FILM"

    my_file = open_sismo("D:/Article/Simulation/Test_article/FILM", 75)

    l_u_norm, l_u_x = norm_computation(my_file)
    l_x = np.arange(0.0, 4, 0.001)
    plt.plot(l_x, l_u_x, color="blue", label='reference model')
    plt.xlabel("t(s)")
    plt.rc('font', size=18)
    plt.ylabel('amplitude')
    plt.show()
