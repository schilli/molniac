import sys, os
import numpy as np
import mdtraj as md
import mdtraj.core.trajectory
import mdtraj.core.topology
from mdtraj.core.trajectory import load_topology


def load(filename):
    mdtrajtrj = md.load(filename)
    trj       = Trajectory(other=mdtrajtrj)
    return trj


class Trajectory(mdtraj.core.trajectory.Trajectory):
    """
    Expand on mdtraj trajectory to make it mutatable
    """

    def __init__(self, xyz=None, topology=None, time=None, unitcell_lengths=None, unitcell_angles=None, other=None):

        # if created like the parent class
        if xyz is not None and topology is not None:
            top = Topology(topology)
            super(Trajectory, self).__init__(xyz, top, time=time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

        # if created from parent class instance
        elif other is not None:
            top = Topology(other._topology)
            super(Trajectory, self).__init__(other.xyz, top, time=time, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)



    def mutate_protein_residue(self, oldres, newres):
        """
        Mutate a protein residue into another protein residue.
        All but the first frame will be lost in the process.

        oldres: (string)
            A selection string that can be parsed by mdtraj topology.select()
        newres: (string)
            Either a three letter amino acid name or a filename that contains a residue
        """

        # get indices of old residue
        oldres_indices = self.topology.select(oldres)

        # make sure it is indeed a single residue
        nresidues = len(set([self.topology.atom(i).residue for i in oldres_indices]))
        if nresidues != 1:
            print("\"{}\" does not evaluate to a single residue, but rather to {}".format(oldres, nresidues))
            sys.exit(1)

        # load new residue as mdtrj trajectory
        if os.path.isfile(newres):
            newrestrj = md.load(newres)
        else:
            datapath = os.path.dirname(os.path.realpath(__file__)) + "/data"
            newres   = newres.split('/')[0].lower()
            newrestrj   = md.load(datapath + "/" + newres + ".pdb")


        ## overlay backbones of both residues

        oldresidue = self.top.atom(oldres_indices[0]).residue
        newresidue = newrestrj.top._residues[0] 

        resSeq = self.top.atom(oldres_indices[0]).residue.resSeq
        backbone_names = ['N', 'CA', 'C']
        if 'CB' in [atom.name for atom in newresidue.atoms] and 'CB' in [atom.name for atom in oldresidue.atoms]:
            backbone_names.append('CB')
        elif 'HA' in [atom.name for atom in newresidue.atoms] and 'HA' in [atom.name for atom in oldresidue.atoms]:
            backbone_names.append('HA')

        # ensure correct ordering of backbone atoms
        old_bb_indices = [self.top.select("resSeq {} and name {}".format(resSeq, name))[0] for name in backbone_names]
        new_bb_indices = [newrestrj.top.select("name {}".format(name))[0] for name in backbone_names]

        newrestrj.superpose(self, frame=0, atom_indices=new_bb_indices, ref_atom_indices=old_bb_indices)

        # set dihedrals by setting positions of Backbone O and H atoms
        old_O_ndx = self.top.select('resSeq {} and name O'.format(resSeq))[0]
        new_O_ndx = newrestrj.top.select('name O')[0]
        newrestrj._xyz[0,new_O_ndx,:] = self._xyz[0,old_O_ndx,:]
        try:
            old_H_ndx = self.top.select('resSeq {} and name H'.format(resSeq))[0]
            new_H_ndx = newrestrj.top.select('name H')[0] 
            newrestrj._xyz[0,new_H_ndx,:] = self._xyz[0,old_H_ndx,:]
        except IndexError:
            pass

        # If its a C-terminal residue, add OXT atom
        try:
            self.topology.residue(oldresidue.index+1)
        except IndexError:
            # add to topology
            atom_O = newresidue.atoms_by_name('O').next()
            newrestrj.top.add_atom('OXT', atom_O.element, newresidue, serial=newresidue._atoms[-1].serial+1)

            # compute position
            index_O  = newrestrj.top.select('name O')[0]
            index_C  = newrestrj.top.select('name C')[0]
            index_CA = newrestrj.top.select('name CA')[0]
            pos_O    = newrestrj.xyz[0,index_O,:]
            pos_C    = newrestrj.xyz[0,index_C,:]
            pos_CA   = newrestrj.xyz[0,index_CA,:]
            s        = pos_C - pos_CA
            v        = pos_O - pos_C
            cp       = np.dot(v, s) / np.dot(s, s)
            pos_OXT  = pos_O + 2*cp*s - 2*v

            # add coordinates
            xyz = np.zeros([1,newrestrj.xyz.shape[1]+1,3], dtype=self._xyz.dtype)
            xyz[0,:-1,:] = newrestrj._xyz
            xyz[0, -1,:] = pos_OXT
            newrestrj._xyz = xyz


        
        ## Update coordinates

        # set up some atom counters
        natoms_oldres = oldres_indices.shape[0]
        natoms_newres = newrestrj.n_atoms
        natoms_before = oldres_indices.min() # number of atoms before the residue starts
        natoms_final  = self.n_atoms - natoms_oldres + natoms_newres

        # copy over unaffected coordinates and insert new coordinates
        # if we have more than a single frame loaded, all but the first will be discarded
        xyz = np.zeros([1, natoms_final, 3], dtype=self._xyz.dtype)
        xyz[0,:natoms_before                           ,:] = self._xyz[0,:natoms_before              ,:]
        xyz[0,natoms_before+natoms_newres:             ,:] = self._xyz[0,natoms_before+natoms_oldres:,:]
        xyz[0,natoms_before:natoms_before+natoms_newres,:] = newrestrj._xyz
        self._xyz = xyz

        # set (chain, index, segment_id) of new residue
        newresidue.chain  = oldresidue.chain
        newresidue.index  = oldresidue.index
        newresidue.resSeq = oldresidue.resSeq
        try:
            oldresidue.segment_id = newresidue.segment_id
        except AttributeError:
            # before segment_id attribute was added
            pass

        # replace residue in chain and topology
        for i, r in enumerate(self.top._residues):
            if r is oldresidue:
                self.top._residues[i] = newresidue
                break
        for i, r in enumerate(newresidue.chain._residues):
            if r is oldresidue:
                newresidue.chain._residues[i] = newresidue
                break

        # replace atoms in topology
        self.top._atoms = self.top._atoms[:natoms_before] + newrestrj.top._atoms + self.top._atoms[natoms_before+natoms_oldres:]

        # renumber atoms
        for i, atom in enumerate(self.top._atoms):
            atom.index  = i
            atom.serial = i+1

        # regenerate bonds
        self.top._bonds = []
        self.top.create_standard_bonds()
        self.top.create_disulfide_bonds(self.xyz[0,:,:])






class Topology(mdtraj.core.topology.Topology):
    """
    Expand on mdtraj topology to make it mutatable
    """

    def __init__(self, other):
        super(Topology, self).__init__()

        for chain in other.chains:
            c = self.add_chain()
            for residue in chain.residues:
                try:
                    r = self.add_residue(residue.name, c, residue.resSeq, residue.segment_id)
                except AttributeError:
                    # before segment_id attribute was added
                    r = self.add_residue(residue.name, c, residue.resSeq)
                for atom in residue.atoms:
                    self.add_atom(atom.name, atom.element, r,
                                 serial=atom.serial)

        for a1, a2 in other.bonds:
            self.add_bond(a1, a2)


