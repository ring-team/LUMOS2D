"""
Association Scientifique pour la Geologie et ses Applications (ASGA)

Copyright (c) 2020 ASGA and Universite de Lorraine. All Rights Reserved.

 Author: Capucine Legentil - capucine.legentil@univ-lorraine.fr

 """
import subprocess


def mmg2d(l_args):
    """
    Run mmg2d executable
    @param l_args: list of inputs and options that is written in command line.
    """
    subprocess.run(l_args)
    print("MMG2D WRITING COMPLETED")


if __name__ == '__main__':
    mmg2d(["/home/legentil/Programmation/ringlab/third_parties/mmg/build/bin/mmg2d_O3", "-ls", "500", "-nomove",
           "-noinsert", "-noswap",
           "/home/legentil/Bureau/Test_jupyter/mesh_merge1234_2d.mesh"])
