"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """

import pandas as pd
import seaborn as sns
from simulation_comparison.L2_norm import *


def scatter_plot(l_x, l_delta_z, l_L2):
    """
    Draw the scatter plot L2-norm waveform difference for each receiver as a function
    of contact depth difference Delta_z (results are in percentages).
    The color represents the L2-norm value at one receiver for one Delta_z.
    The Delta_z value are negative when the contact of a model is below the contact of the reference model.
    @param l_x: x position of receivers
    @param l_delta_z: list of tested error depths
    @param l_L2: list of l2-norm value for each model
    """
    compil = {"x": l_x, "delta_z": l_delta_z, "l2_norm": l_L2}
    dtframe = pd.DataFrame(compil)
    error = dtframe.pivot("delta_z", "x", "l2_norm")
    plt.rc('font', size=18)
    sns.heatmap(error)
    plt.ylabel(r'$\Delta_z$')

    plt.savefig("/home/legentil/Documents/Article/first_group_simu/image_complete/scatter_x730.png")


def sum_plot(l_delta_z, l_sum):
    """
    Plot the sum of L2-norm as a function of Delta_z.
    @param l_delta_z: list of tested error depths
    @param l_sum: list of the sum of l2-norm values along receivers
    """
    plt.xlabel("$\delta_z$ (m)")
    plt.plot(l_delta_z, l_sum, label="sum difference", marker="o")
    plt.rc('font', size=18)
    plt.savefig("/home/legentil/Documents/Article/first_group_simu/image_complete/sum_x730.png")


def mean_plot(l_delta_z, l_mean):
    """
    Plot the mean of L2-norm as a function of Delta_z.
    @param l_delta_z: list of tested error depths
    @param l_mean: list of the mean of l2-norm values along receivers
    """
    plt.xlabel("$\delta_z$ (m)")
    plt.plot(l_delta_z, l_mean, label="sum difference", marker="o")
    plt.rc('font', size=18)
    plt.savefig("/home/legentil/Documents/Article/first_group_simu/image_complete/mean_x730.png")


def sismos(l_z):
    """
    Compute the l2-norm for each models
    Plot one seismogram for each model at the receiver where the l2-norm value is maximal.
    @param l_z: list of tested depths
    @return: l_dz: list of tested depth error for eache receiver,
             l_dz_sum: list of tested depth error,
             l_l2_norm_tot: list of l2-norm values for each receiver at each tested depth,
             l_sum: list of the sum of l2-norm values along receivers,
             l_mean:  list of the mean of l2-norm values along receivers.
    """
    l_x = np.arange(0.0, 4, 0.001)
    l_l2_norm_tot = []
    z_ref = -730
    l_dz = []
    l_sum = []
    l_mean = []
    l_dz_sum = [z - z_ref for z in l_z]
    root_path = "/home/legentil/Documents/Article/bdy_ok_group"
    dirname_ref = root_path + "/ls" + str(z_ref) + "/FILM"
    for z in l_z:
        if z != z_ref:
            plt.figure()
            l_dz += [z - z_ref] * 150
            dir_name = root_path + "/ls" + str(z) + "/FILM"
            print(dir_name)
            l_l2_norm, max_id, l_u_norm_ref = norm_l2(dirname_ref, dir_name, 150, 0.001)
            l_l2_norm_tot += l_l2_norm
            l_sum.append(np.sum(l_l2_norm))
            l_mean.append(np.mean(l_l2_norm))
            my_file_ref = open_sismo(dirname_ref, max_id)
            m_file = open_sismo(dir_name, max_id)
            l_u_norm_ref, l_u_x_ref = norm_computation(my_file_ref)
            l_u_norm, l_u_x = norm_computation(m_file)
            plt.plot(l_x, l_u_x_ref, linestyle='dashed', color="black", label='reference model')
            plt.plot(l_x, l_u_x, color="deepskyblue", label='delta =' + str(z - z_ref) + 'm')
            plt.xlabel("t (s)")
            plt.ylabel("$u_y$ (m)")
            axes = plt.gca()
            axes.set_xlim(1.5, 4)
            axes.set_ylim(-9e-15, 9e-15)
            plt.legend()
            plt.savefig("/home/legentil/Documents/Article/bdy_ok_group/image_complete/" + str(z) + "_x" + str(
                z_ref) + ".png", dpi=600)
        else:
            l_sum.append(0)
            l_mean.append(0)

    return l_dz, l_dz_sum, l_l2_norm_tot, l_sum, l_mean


def quality_vs_computation_time():
    """
    Plot the simulation run time as a function of the mesh quality
    """
    l_quality = [2, 2.35, 2.57, 2.77, 2.95, 4.76, 4.85, 5.1, 6.01, 6.1, 6.37, 6.89, 7.36, 7.99, 8.54, 11.31, 17.95,
                 22.55, 22.76, 24.84, 25.51]
    l_time = [46376, 42651, 38380, 35271, 33718, 21759, 21300, 20847, 19620, 18231, 16973, 17883, 20298, 18469, 17321,
              13674, 11508, 11015, 10809, 11514, 11348]
    l_depth = [-2655, -2685, -2835, -2765, -2735, -2705, -2695, -2715, -2745, -2675, -2665, -2825, -2785, -2725, -2645,
               -2815, -2795, -2845, -2775, -2805, -2755]
    dico = {l_depth[i]: [l_time[i], l_quality[i]] for i in range(len(l_depth))}
    tri1 = sorted(dico.items(), key=lambda t: t[0])
    l_depth_1 = [model[0] for model in tri1]
    l_quality_1 = [model[1][1] for model in tri1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pl1 = ax.plot(l_time, l_quality, marker="o", color="black", label='Computation time')
    plt.ylabel("Quality", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Computation time (s)", fontsize=18)
    axes = ax.twiny()
    axes.set_xlabel("Contact depth (m)", fontsize=18)

    axes.plot(l_depth_1, l_quality_1, label='Contact depth', linestyle='dashed', marker="+")

    fig.legend()

    fig = plt.figure()
    plt.plot(l_depth_1, l_quality_1, label='Contact depth', linestyle='dashed', marker="+")
    plt.ylabel("Quality", fontsize=18)
    plt.xlabel("Contact depth (m)", fontsize=18)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(l_quality, l_time, label='Contact depth', marker=".", color="black")
    plt.ylabel("Computation time (s)", fontsize=18)
    plt.xlabel("Quality", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.show()


if __name__ == '__main__':
    l_z = [-730, -720, -710, -700, -690, -680, -670, -660, -650, -640, -630, -620, -610, -600, -590, -580, -570, -560,
           -550, -540, -530]
    l_x_dist = [i for i in range(0, (150) * 100, 100)] * (len(l_z) - 1)

    l_dz, l_dz_sum, l_l2_norm_tot, l_sum, l_mean = sismos(l_z)
    plt.figure(figsize=(18, 12))
    scatter_plot(l_x_dist, l_dz, l_l2_norm_tot)
    plt.figure(figsize=(18, 12))
    sum_plot(l_dz_sum, l_sum)
    plt.figure(figsize=(18, 12))
    mean_plot(l_dz_sum, l_mean)
    plt.show()

    quality_vs_computation_time()
