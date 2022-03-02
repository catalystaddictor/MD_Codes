import sys
import numpy as np 
import os

nanometer=10.0  # nenometer to Angstrom 
angstrom = 1.0 
class Molecule:
    def __init__(self, coord_file):
        self.coord_file = coord_file
        self._read()

    def _read(self):
        self._read_pdb()

    def _malloc(self, nAtoms):
        self.nAtoms = nAtoms
        self.coords = np.zeros((3, nAtoms), order="C")
        self.symbols = ["C"] * nAtoms
        self.resnames = ["MOL"] * nAtoms
        self.resids = [1] * nAtoms
        self.atomids = [i + 1 for i in range(nAtoms)]
        self.box = np.zeros(3)

    def set_box(self, box):
        self.box = np.array(box)

    def _read_pdb(self):

        with open(self.coord_file, "r") as file_handle:
            lines = file_handle.readlines()

        # find number of atoms
        nAtoms = 0
        for line in lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                nAtoms += 1
        assert nAtoms > 0
        self._malloc(nAtoms)

        # read coordinates
        atom_id = 0
        for line in lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                self.symbols[atom_id] = line[12:16].strip()
                self.resnames[atom_id] = "MOL"
                self.resids[atom_id] = 1
                self.x[atom_id] = float(line[30:38]) * angstrom
                self.y[atom_id] = float(line[38:46]) * angstrom
                self.z[atom_id] = float(line[46:54]) * angstrom

                atom_id += 1

    def write(self, out_file, write_mode="w"):
        self._write_gro(out_file, write_mode)

    def _write_gro(self, out_file, write_mode):
        file_handler = open(out_file, write_mode)

        file_handler.write("%s\n" % ("Created by MolAligner"))
        file_handler.write("%-10d\n" % (self.nAtoms))

        groFMT = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"
        for i in range(self.nAtoms):
            atomid = i + 1
            atomid = atomid
            resid = self.resids[i]
            file_handler.write(
                groFMT
                % (
                    resid,
                    self.resnames[i],
                    self.symbols[i],
                    atomid,
                    self.x[i] / nanometer,
                    self.y[i] / nanometer,
                    self.z[i] / nanometer,
                )
            )
        file_handler.write(
            #9*"%10.5f"+"\n"
            "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n"
            % (
                self.box[0],
                self.box[1],
                self.box[2],
                self.box[3],
                self.box[4],
                self.box[5],
                self.box[6],
                self.box[7],
                self.box[8],
            )
        )
        file_handler.close()


    @property
    def x(self):
        return self.coords[0, :]

    @property
    def y(self):
        return self.coords[1, :]

    @property
    def z(self):
        return self.coords[2, :]



def main(name):
    with open('CyclohexaneSolvation-pos-1.pdb','r') as F: 
        lines = F.readlines() 

    nlines = len(lines)
    nlines_per_frame= 428 
    nFrames= int(nlines/nlines_per_frame)
    
    out_traj = "traj.gro"
    os.system("rm -f "+ out_traj)

    print(f"Total frames found: {nFrames}")


    box = [
        2.21290,   
        1.19137,   
        4.00000,   
        0.00000,   
        0.00000,  
        -0.34816,   
        0.00000,   
        0.00000,   
        0.00000,
        ]

    for i in range(nFrames):
        print(f"Converting frame number: {i+1}/{nFrames}")
        start=i*nlines_per_frame 
        end = (i+1)*nlines_per_frame
        pdb_frame = f"frame-{i}.pdb"
        with open(pdb_frame,'w') as F:
            for j in range(start,end):
                F.write(lines[j])

        mol = Molecule(pdb_frame)
        mol.set_box(box)
        mol.write('traj.gro','a') 
        os.system(f"rm -f {pdb_frame}")

infile=
main(infile)
