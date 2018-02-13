from lib.fulleren_create import *
from subprocess import Popen
import sys

create_flag = ("new" in sys.argv) or ("all" in sys.argv)
tleap_flag = ("tleap" in sys.argv) or ("all" in sys.argv)
min_flag = ("min" in sys.argv) or ("all" in sys.argv)
just_run_flag = ("just_run" in sys.argv)

N_C = 500
N_H = 200

if create_flag:
    '''
    We read fulleren from file in external_source dirictory
    Then add in interior of fulerene hydrogen
    And in the end we write fullerene into mol2 file and every H atom in separated mol2 file (use antechamber)
    '''
    fulleren, r = get_structure(N_C)
    hydrogen = add_hydrogen(r, N_H)
    filepath = "./input/"
    try:
        os.stat(filepath)
    except:
        os.mkdir(filepath)

    write_mol2(fulleren, "C", filename=filepath+"C{:d}.mol2".format(N_C), mol = "MOL", residue_n = 1)
    write_mol2_by_one(hydrogen, "H", filename=filepath+"H{:d}.mol2" , mol = "MOL", residue_n = 1)

if tleap_flag:
    '''
    We create tleap script to create topology file to hydrogen
    '''
    tleap = open("tleap.in", "w")
    tleap.write("source leaprc.gaff\n"
            "loadamberparams external_source/parm10gaff.dat\n"
            "x = loadmol2 {:s}C{:d}.mol2\n".format(filepath, N_C))
    for i in range(len(hydrogen)):
        tleap.write("h{:d} = loadmol2 {:s}H{:d}.mol2\n".format(i, filepath, i) +
                "x = combine {{ x h{:d} }}\n".format(i) )
    tleap.write("savepdb x fulleren.pdb\n"
            "saveamberparm x input/prmtop input/inpcrd\n"
            "quit")
    tleap.close()

    #run creating topology file for amber
    p = Popen(["tleap", "-f", "tleap.in"])
    p.wait()
    print "tleap ---OK---"

if min_flag:
    #start minimization
    p = Popen(["sander", "-O", "-i", "minimization.in", "-p", "input/prmtop", "-c", "input/inpcrd", "-r", "01_min.rst", "-inf", "01_min.mdinfo","-o", "01_min.out"])
    p.wait()
    print "minimization ---OK---"

    #start mdrun
    p = Popen(["sander", "-O", "-i", "mdrun.in", "-p", "input/prmtop", "-c", "01_min.rst", "-x", "mdcrd", "-inf", "02_mdrun.mdinfo", "-o", "02_mdrun.out"])
    p.wait()
    print "molecular dynamic ---OK---"

if just_run_flag:
    p = Popen(["sander", "-O", "-i", "mdrun.in", "-p", "input/prmtop", "-c", "input/prmtop", "-x", "mdcrd", "-inf", "02_mdrun.mdinfo", "-o", "02_mdrun.out"])
    p.wait()
    print "molecular dynamic ---OK---"

