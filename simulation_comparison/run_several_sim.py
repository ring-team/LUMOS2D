"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
import os
import subprocess
from lib_cavity.write_local_mesh_file import build_2_local_meshes, merge
from lib_cavity.mesh_2d_file import write_2D_mesh_file_gmsh, write_2D_mesh_req_file
from simulation_comparison.condition_limites import write_boundary_condition
import shutil


def l_depth(z_min, z_max, gap):
    """
    Create the list of tested depths
    @param z_min: minimal depth
    @param z_max: maximal depth
    @param gap: step between two different tested depths
    @return: list of tested depths
    """
    l_depth = [z for z in range(z_min, z_max + gap, gap)]
    return l_depth


def mmg2d(l_args):
    """
    Run mmg2d executable
    @param l_args: list of option for mmg2d
    """
    subprocess.run(l_args)


def create_path_z(l_depth, root_path):
    """
    Create folder to save models at each depth
    @param l_depth: list of tested depths
    @param root_path: folder that will contain all the models
    """
    for z in l_depth:
        dir_name = root_path + "/ls" + str(z)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        os.chdir(dir_name)


def local_modif(dir_name_z, z, init_mesh_path, mmg_path, modified_regions, kept_regions):
    """
    Create the modified model with a contact depth z
    @param dir_name_z: folder name that will contain the model
    @param z: the contact depth
    @param init_mesh_path: entire file path of the initial meshed model
    @param mmg_path: mmg2d executable path
    @param modified_regions: list of indexes of modified regions
    @param kept_regions: list of indexes of regions that is not impacted by the modification
    """
    init_mesh_name = init_mesh_path.split("/")[-1]
    mesh_out_name = dir_name_z + "/" + init_mesh_name[:-4] + "o.mesh"
    print(mesh_out_name)

    # level_set discretization
    mmg2d([mmg_path, "-ls", str(z), "-nomove", "-noinsert", "-noswap",
           init_mesh_path, mesh_out_name])

    # regions isolation
    filepath_to_write_cut = dir_name_z + "/reg_cut.mesh"
    filepath_to_write_other = dir_name_z + "/reg_other.mesh"
    build_2_local_meshes(mesh_out_name, filepath_to_write_cut, filepath_to_write_other,
                         [modified_regions, kept_regions])

    # RequiredEdges
    filepath_to_write_cut_req = dir_name_z + "/reg_cut_req.mesh"
    write_2D_mesh_req_file(filepath_to_write_cut, filepath_to_write_cut_req)

    # Remeshing
    mmg2d([mmg_path, "-hmin", "90", "-hmax", "90", dir_name_z + "/reg_cut_req.mesh"])

    # merge
    filepath_mesh_merge = dir_name_z + "/mesh_merge.mesh"
    filepath_final_mesh = dir_name_z + "/mesh_sim.mesh"
    merge(dir_name_z + "/reg_cut_req.o.mesh", filepath_to_write_other, filepath_mesh_merge)
    write_2D_mesh_file_gmsh(filepath_mesh_merge, filepath_final_mesh)

    # tetgen export
    mmg2d([mmg_path, "-nomove", "-noinsert", "-noswap",
           filepath_final_mesh, dir_name_z + "/mesh_sim.1.node"])

    # rename
    os.rename(dir_name_z + "/mesh_sim.1.node", dir_name_z + "/mesh_sim.1_cl.node")

    # boudary_conditions
    write_boundary_condition(dir_name_z + "/mesh_sim.1_cl.node",
                             dir_name_z + "/mesh_sim.1.node",
                             "-5291.109375000000000", "10949.264648438000222", "-1837.562988281199978",
                             "3247.130371093699978")

    # cleaning
    os.remove(filepath_to_write_cut)
    os.remove(filepath_to_write_other)
    os.remove(filepath_to_write_cut_req)
    os.remove(filepath_mesh_merge)
    os.remove(dir_name_z + "/reg_cut_req.o.mesh")
    os.remove(dir_name_z + "/reg_cut_req.o.sol")


def local_modif_group(l_depth, root_path, init_mesh_path, mmg_path, modified_regions, kept_regions):
    """
    Create all the models, one for each depth
    @param l_depth: list of tested depths
    @param root_path: folder that will contain all the models
    @param init_mesh_path: entire file path of the initial meshed model
    @param mmg_path: mmg2d executable path
    @param modified_regions: list of indexes of modified regions
    @param kept_regions: list of indexes of regions that is not impacted by the modification
    """
    for z in l_depth:
        dir_name_z = root_path + "/ls" + str(z)
        local_modif(dir_name_z, z, init_mesh_path, mmg_path, modified_regions, kept_regions)
        shutil.copy(root_path + "/mesh_sim.model", dir_name_z + "/mesh_sim.model")
        name_param = dir_name_z + "/param_simu" + str(z) + ".txt"
        param_simu(dir_name_z, name_param)


def param_simu(dir_name_z, filename):
    """
    Create parameter file for seismic simulation (with Hou10ni2D)
    @param dir_name_z: folder name that will contain the model
    @param filename: name of the parameter file
    """
    param_file = open(filename, "w")
    param_file.write("restart=1 \nrestart_name=.restart \ndt_restart=0.1 \ndim=2 \n")
    param_file.write("mesh=" + dir_name_z + "/mesh_sim \n")
    param_file.write("kind_of_medium=0 \n\
    local_time_stepping=0 \n\
    type_data_input=.model \n\
    p_adaptative=1 \n\
    degree_max=6 \n\
    helmholtz=0 \n\
    simulation_time=4 Temps de simulation (en seconde) \n\
    dt_snap=0.02 time for the snapshots \n\
    type_src=3 Plane wave on diffracting obstacle (1) (puis mettre  theta incident) or point source (2) \
    or line of sources, one computation for all sources (3) or line of sources, one computation for each \
    source (4) or circle of source one computation for all sources (5) or circle of sources, \
    one computation for each source (6)  or ellipse of source one computation for all sources (7) \
    or ellipse of sources, one computation for each source (8) \n\
    nb_src=1 Nombre de sources \n\
    first_source_position=-3160 2470 Starting point \n\
    last_source_position=-3160 2470  end point (on place 4 sources équiréparties entre le premier et le dernier point) \n\
    type_src_time=1 Type de source : 1=Ricker, -1= fichier fourni par l'utilisateur \n\
    fpeak=15 Si paramètre précédent=1: fréquence centrale de la source \n\
    curved=0 Do you want curved elements : pas utile ici \n\
    symmetry=1 Do you want to use mumps with symmetry : pas utile ici \n\
    type_rcv=3 Pour les sismos : 3 correspond à des lignes de récepteurs \n\
    nb_lines=1 nombre de  lignes \n\
    nb_rcv_on_lines=150 nb de récepteurs sur la ligne \n\
    line_params=-5291 3245 10949 3245 \n\
    dt_sismo=0.001 pas de temps pour les sismo \n\
    type_cla=0 inutile \n\
    use_pml=.false. \n\
    #.FALSE. inutile \n\
    type_rtm=0 inutile \n\
    type_analytic=0 inutile \n")
    param_file.close()


def write_script_simu(filename, l_depth, root_path):
    """
    Create shell script to run all the simulations
    @param filename: name of the script file
    @param l_depth: list of tested depths
    @param root_path: folder that contains all the models
    """
    script_file = open(root_path + "/" + filename, "w")
    script_file.write("#! /bin/bash\n")
    script_file.write("cd /home/legentil/Programmation/hou10ni2d/build/\n")

    for z in l_depth:
        dir_name_z = root_path + "/ls" + str(z)
        script_file.write("mpirun -np 8 hou10ni.out " + dir_name_z + "/param_simu" + str(z) + ".txt\n")
        script_file.write("cp -r /home/legentil/Programmation/hou10ni2d/build/FILM/ " + dir_name_z + "/FILM/\n\n")

    script_file.close()


if __name__ == '__main__':
    mmg_path = "/home/legentil/Programmation/ringlab/third_parties/mmg/build/bin/mmg2d_O3"
    # l_depth_1 = l_depth(-730, -530, 10)
    # create_path_z(l_depth_1, "/home/legentil/Documents/Article/just_mesh")
    # local_modif("/home/legentil/Documents/Article/bdy_ok_group/ls-590b", -590,
    # "/home/legentil/Documents/Article/bdy_ok_group/rescale_test.mesh", mmg_path,
    # [6, 10], [1, 2, 3, 4, 5, 7, 8, 9])
    # local_modif_group(l_depth_1,"/home/legentil/Documents/Article/just_mesh", "/home/legentil/Documents/Article/just_mesh/article_model_f.mesh", mmg_path, [6, 10], [1, 2, 3, 4, 5, 7, 8, 9])
    # write_script_simu("simu.sh", l_depth_1, "/home/legentil/Documents/Article/bdy_ok_group")

    local_modif("/home/legentil/Bureau/raffinement/ls-700-nofem", 710,
                "/home/legentil/Bureau/raffinement/mesh_merge1234_2d.mesh", mmg_path,
                [2, 6], [1, 3, 4])
