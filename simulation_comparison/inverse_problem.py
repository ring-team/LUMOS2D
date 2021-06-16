"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
from simulation_comparison.L2_norm import *
from math import *
import matplotlib.cbook as cbook
from mpl_toolkits.axes_grid1 import make_axes_locatable


def likelihood(dirname_ref, dirname, nb_rec, time_step):
    """
    Compute the likelihood value for a model
    @param dirname_ref: folder containing the reference model simulation results
    @param dirname: folder containing the simulation results
    @param nb_rec: number of receivers
    @param time_step: time step
    @return: likelihood value
    """
    l_L2_norm, max_id, l_u_norm_ref_square = norm_l2(dirname_ref, dirname, nb_rec, time_step)
    delta_d_i = [(i / 100) ** 2 for i in l_L2_norm]
    sum_x_delta_d_i = sum(delta_d_i)
    exposant = -0.5 * (sum_x_delta_d_i)
    sigma = sqrt(sum(l_u_norm_ref_square))
    likelihood_val = (exp(exposant)) / (sigma * sqrt(2 * pi))
    return likelihood_val


def densite_priori(z):
    """
    Return the prior probability density
    In our case: 1/200 if 402.5m <= z <= 602.5m; 0 otherwise
    @param z: depth
    @return: prior probability
    """
    if 402.5 <= z <= 602.5:
        return 1 / 200
    else:
        return 0


def proba_posteriori(l_z, z_ref):
    """
    Return the posterior probability density
    @param l_z: list of tested depths
    @param z_ref: reference depth
    @return: list of probability of all depths
    """
    l_dz = []
    l_likelihood = []
    l_prior = []

    ### denominateur
    root_path = "/home/legentil/Documents/Article/first_group_simu"
    dirname_ref = root_path + "/ls" + str(z_ref) + "/FILM"
    for z in l_z:
        l_dz += [z - z_ref] * 150
        if z >= -650:
            root_path = "/home/legentil/Documents/Article/bdy_ok_group"
        dir_name = root_path + "/ls" + str(z) + "/FILM"
        likelihood_val = likelihood(dirname_ref, dir_name, 150, 0.001)
        l_likelihood.append(likelihood_val)
        z_real = z + 1132.5
        l_prior.append(densite_priori(z_real))

    l_denominateur = [l_likelihood[i] * l_prior[i] for i in range(len(l_likelihood))]
    integrale = approximation_trapeze(l_denominateur, 10)

    l_posterior = [(l_prior[i] * l_likelihood[i]) / integrale for i in range(len(l_prior))]

    fig, axs = plt.subplots(nrows=2, gridspec_kw={'height_ratios': [1, 1]})
    fig.set_figheight(10)
    fig.set_figwidth(15)
    l_posterior_final, l_final_z, l_prior_final = real_z(-5084, 0, l_z, l_posterior)

    with cbook.get_sample_data("/home/legentil/Images/mode_fl.PNG") as datafile:
        image = plt.imread(datafile)

    # extent = (-5291, 10949, -1837, 3247)
    extent = (0, 16240, -5084, 0)

    axs[0].imshow(image, extent=extent)
    axs[0].set_ylabel("z (m)")
    axs[0].set_xlabel("x (m)")

    divider = make_axes_locatable(axs[0])
    axright = divider.append_axes("right", size=1.2, pad=0.8, sharey=axs[0])

    axright.plot(l_posterior_final, l_final_z, label="posterior")
    axright.plot(l_prior_final, l_final_z, linestyle='dashed', label="prior")
    axright.set_ylabel("z (m)")
    axright.set_xlabel("Probability")

    axs[1].plot(l_posterior_final, l_final_z)
    axs[1].set_ylim((200 - 3247, 820 - 3247))
    axs[1].plot(l_prior_final, l_final_z, linestyle='dashed')
    axs[1].set_ylabel("z (m)")
    axs[1].set_xlabel("Probability")
    fig.legend()
    plt.show()
    fig.savefig("/home/legentil/Documents/Article/bdy_ok_group/posterior" + "_y" + str(z_ref) + ".png", dpi=600)
    return l_posterior


def real_z(zmin, zmax, l_z, l_posterior):
    """
    Translation of depth for visualization
    @param zmin: real minimal model depth
    @param zmax: real maximal model depth
    @param l_z: list f tested depths
    @param l_posterior: posterior probability density
    @return: l_res_final: translated posterior probability density,
             l_final_z: translated list of depths,
             l_prior_final: translated prior probability density,
    """
    l_real_z = [i + 1132.5 for i in l_z]
    before_zmin = [i for i in range(int(zmin), int(min(l_real_z)), 10)]
    res_before = [0] * len(before_zmin)
    after_zmax = [i for i in range(int(max(l_real_z)), int(zmax), 10)]
    res_after = [0] * len(after_zmax)
    l_prior = [1 / 200] * len(l_z)
    l_final_z = before_zmin + l_real_z + after_zmax
    l_prior_final = res_before + l_prior + res_after
    l_res_final = res_before + l_posterior + res_after

    return l_res_final, l_final_z, l_prior_final


if __name__ == '__main__':
    sismo_ref = "/home/legentil/Documents/Article/first_group_simu/ls-730/FILM"
    sismo = "/home/legentil/Documents/Article/bdy_ok_group/ls-660/FILM"

    l_z = [-730, -720, -710, -700, -690, -680, -670, -660, -650, -640, -630, -620, -610, -600, -589, -580, -570, -560,
           -550, -540, -530]
    z_ref = -660

    print(proba_posteriori(l_z, z_ref))
