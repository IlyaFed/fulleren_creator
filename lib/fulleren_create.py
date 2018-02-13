import numpy as np
import random
import requests
from subprocess import Popen
import os

def download_structure(N_C = 0):
    '''
    download_structure(N_C = 0)
        downloading structure from site, and adding it in fulleren structure
    return: fulleren
    '''
    if (N_C == 0):
        N_C = int(raw_input("set number of C in fulleren, see available http://www.nanotube.msu.edu/fullerene/fullerene-isomers.html"))
    url = 'http://www.nanotube.msu.edu/fullerene/C'+str(N_C)+'/C'+str(N_C)+'-0.xyz'
    xyz = requests.get(url, allow_redirects=True)

    fulleren = list()

    N_C = 0

    num_line = 0
    for line in xyz.content.split("\n"):#open(example_file, "r")):
        if num_line == 0:
            N_C = int(line)
        if num_line == N_C + 2:
            break
        if num_line > 1:
            fulleren.append(np.array([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]))
        num_line += 1

    r = np.linalg.norm(fulleren[0])

    return fulleren, r

def modificate_structure(fulleren, r = 0):
    '''
    modificate_structure(fulleren, r = 0)
        make fulleren in neccessary size
    return: fulleren, r
    '''
    if r == 0:
        r = int(raw_input("Put radius(existing = "+str(np.linalg.norm(fulleren[0]))+"): "))
    if r == 0:
        r = np.linalg.norm(fulleren[0])


    coeff = 1.0*r/np.linalg.norm(fulleren[0])

    for i in range(len(fulleren)):
        fulleren[i] *= coeff
    return fulleren, r

def get_structure(N_C = 0):
    '''
        getting structure from file, and adding it in fulleren structure
    return: fulleren
    '''


    fulleren_file = "external_source/C{:d}.xyz".format(N_C)

    fulleren = []

    N_C = 0

    num_line = 0
    for line in open(fulleren_file, "r"):
        if num_line == 0:
            N_C = int(line)
        if num_line == N_C + 2:
            break
        if num_line > 1:
            fulleren.append(np.array([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]))
        num_line += 1

    r = np.linalg.norm(fulleren[0])

    return fulleren, r

def add_hydrogen(r,  N_H = 0):
    '''
    add_hydrogen(r, N_H = 0)
        adding hydrogen in fulleren
    return: hydrogen
    '''
    if N_H == 0:
        N_H = int(raw_input("hydrogen quantity:"))
    r_H = 0.8*r
    grid = int(2*r_H/1) # step 1 angstrem

    interior_list = list()

    for x in np.linspace(-r_H, r_H, grid):
        for y in np.linspace(-r_H, r_H, grid):
            for z in np.linspace(-r_H, r_H, grid):
                if np.linalg.norm(np.array([x,y,z])) < r_H:
                    interior_list.append(np.array([x,y,z]))

    #print "interior_list size ", len(interior_list)
    hydrogen = list()

    rand_list = random.sample(range(len(interior_list)), N_H)

    for i in rand_list:
        r = interior_list[i]
        hydrogen.append(r)

    return  hydrogen


def write_lammps_input(fulleren, hydrogen, filename=0):
    '''
    write_lammps_input(fulleren, hydrogen, filename=0)
        write xyz files with fulleren and hydrogen
    return: 0
    '''
    if filename == 0:
        filename = "C"+str(len(fulleren))+"_H"+st(len(hydrogen))+".xyz"

    f = open(filename, "w")


    f.write(str(N_C+N_H)+"\nC"+str(N_C+N_H)+"\n")

    for i in range(len(fulleren)):
        f.write("C "+str(fulleren[i][0])+" "+str(fulleren[i][1])+" "+str(fulleren[i][2]) + "\n")

    for i in range(len(hydrogen)):
        f.write("H "+str(hydrogen[i][0])+" "+str(hydrogen[i][1])+" "+str(hydrogen[i][2]) + "\n")

    f.close()
    return 0


def write_qmmm(fulleren, hydrogen, firstname = 0, secondname = 0):
    '''
    write_qmmm(fulleren, hydrogen, firstname = 0, secondname = 0)
        Create QMMM file input
        write QM file to qm.xyz
        write MM file to mm.xyz
    return: 0
    '''
    f = open("qm.xyz", 'w')
    N_H = len(hydrogen)
    N_C = len(fulleren)

    f.write(str(N_H) +"\nHydrogen Molecule -- Xmol format\n")
    for i in range(len(hydrogen)):
        f.write("H "+str(hydrogen[i][0])+" "+str(hydrogen[i][1])+" "+str(hydrogen[i][2]) + "\n")
    f.close()

    f = open("mm.xyz", 'w')

    f.write(str(N_C) +"\Oxygen Molecule -- Xmol format\n")
    for i in range(len(fulleren)):
        f.write("O "+str(fulleren[i][0])+" "+str(fulleren[i][1])+" "+str(fulleren[i][2]) + "\n")

    f.close()
    return 0


def write_mol2(fulleren, type, filename = "pdb", start_n = 0, add_to_file=False, mol = "MOL", residue_n = 1):
    '''
    :param fulleren: structure with fulleren
    :param filename: ouput name of file
    :return: 0
    '''
    f = open("tmp_to_antechamber_mol2.pdb", "w")

    for i in range(len(fulleren)):
        f.write("{:>4s}  {:5d} {:>2s}{:<2d}{:1s}{:<3s} {:1s}{:<4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:>6s}{:6.2f}      {:<4s}{:>2s}\n".format(
            "ATOM", i+1+start_n, type, 1, "", mol, "", residue_n ,"", fulleren[i][0], fulleren[i][1], fulleren[i][2], "inf", 1, "",type))
    f.write("TER\n")
    f.close()

    p = Popen(["antechamber", "-i", "tmp_to_antechamber_mol2.pdb", "-fi", "pdb", "-o", filename, "-fo", "mol2"])
    p.wait()

    p = Popen(["rm", "tmp_to_antechamber_mol2.pdb"])
    p.wait()

    return len(fulleren)

def write_mol2_by_one(fulleren, type, filename = "pdb", add_to_file=False, mol = "MOL", residue_n = 1):
    '''
    :param fulleren: structure with fulleren
    :param filename: ouput name of file
    :return: 0
    '''

    for i in range(len(fulleren)):
        f = open("tmp_to_antechamber.pdb", "w")
        f.write("{:>4s}  {:5d} {:>2s}{:<2d}{:1s}{:<3s} {:1s}{:<4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:>6s}{:6.2f}      {:<4s}{:>2s}\n".format(
            "ATOM", 1, type, 1, "", mol, "",residue_n,"", fulleren[i][0], fulleren[i][1], fulleren[i][2], "inf", 1, "",type))
        f.close()
        p = Popen(["antechamber", "-i", "tmp_to_antechamber.pdb", "-fi", "pdb", "-o", filename.format(i), "-fo", "mol2"])
        p.wait()

        # Read in the file
        with open(filename.format(i), 'r') as file :
          filedata = file.read()

        # Replace the target string
        filedata = filedata.replace('DU', 'H ')

        # Write the file out again
        with open(filename.format(i), 'w') as file:
          file.write(filedata)

    p = Popen(["rm", "tmp_to_antechamber.pdb"])
    p.wait()

    return len(fulleren)

def test():
    fulleren, r = get_structure(180)
    #fulleren, r = modificate_structure(fulleren, r=6*r)
    hydrogen = add_hydrogen(r, 10)
    #write_amber_inpcrd(fulleren)
    filename = "/Users/scientist/work/amber/qmmm_just_amber/qmmm.pdb"
    #n = write_mol2(fulleren, "C", filename=filename, mol = "MOL", residue_n = 1)
    #write_pdb(hydrogen, "H", filename=filename, mol = "MOL", residue_n = 2, start_n = n, add_to_file=True)
    print "success"
    #write_amber_prmtop(fulleren)
    return 0

