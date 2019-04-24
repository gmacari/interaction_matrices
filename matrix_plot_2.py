#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 22:50:05 2019

@author: andrea
"""

#make some very import stuff

import os, sys
import numpy as np
import oddt
#import matplotlib
from matplotlib import pyplot as plt
from oddt.interactions import (close_contacts, hbonds, hydrophobic_contacts, pi_cation, pi_stacking, salt_bridges)


#handy function to automatize hbond identification
def hbonds_identifier(protein, ligand):

    '''
    input: protein obj and ligand obj from oddt
    return: a list containing [residue's number, residue's name, residue's atomtype, ligand's atomtype,
    ligand's atom nunmber]
    '''
    t = []
    protein_atoms, ligand_atoms, strict = hbonds(protein, ligand)
    for i in range(len(protein_atoms['resnum'])):
        t.append((protein_atoms['resnum'][i],protein_atoms['resname'][i],protein_atoms['atomtype'][i], ligand_atoms['atomtype'][i], ligand_atoms['id'][i]))
    return(t)


#handy function to automatize hydrophobic contact identification
def hydrophobics_identifier(protein, ligand):

    '''
    input: protein obj and ligand obj from oddt
    return: a list containing [residue's number, residue's name, residue's atomtype, ligand's atomtype,
    ligand's atom nunmber]
    '''
    t = []
    protein_atoms, ligand_atoms = hydrophobic_contacts(protein, ligand)
    for i in range(len(protein_atoms['resnum'])):
        t.append((protein_atoms['resnum'][i],protein_atoms['resname'][i],protein_atoms['atomtype'][i], ligand_atoms['atomtype'][i], ligand_atoms['id'][i]))
    return(t)


#return a binary list containing 1 if the interacting residue of  ligand is in the reference array
def feature_count(interaction_array, reference_array):

    '''
    Input: interaction_array = a list of residues interacting with a ligand
            reference_array = a non redundant list of all the residues interacting with all the ligands
    return: a binary array where 1 indicate that a residues is found in the interaction array, else 0
    '''

    feature_count = []
    for tupla in reference_array:
        print(tupla)
        if tupla in interaction_array:
            feature_count.append(1)
        else:
            feature_count.append(0)
    return feature_count




#START
sys.argv.append('/home/andrea/Dropbox/Herbicide_Bounds_xtal_Structures/Docking/1_2J8C_-cofctr_og/2j8c_clean.pdb')
sys.argv.append('/home/andrea/Dropbox/Herbicide_Bounds_xtal_Structures/Docking/5_2J8C_WT-I224S_ng/5a_2j8c_grid_wt/output/UQ0.pdb')
sys.argv.append('/home/andrea/Dropbox/Herbicide_Bounds_xtal_Structures/Docking/5_2J8C_WT-I224S_ng/5a_2j8c_grid_wt/output/atz.pdb')

#Take first argument as the target protein
protein = next(oddt.toolkit.readfile('pdb', sys.argv[1]))
protein.protein = True

#Take the next arguments as the ligands
ligand = []
for arg in sys.argv[2:]:
    ligand.append(next(oddt.toolkit.readfile('pdb', arg)))

#interaction calculation
hblist = []
hclist = []
hbs = []
hcs = []
hbtot = []
hctot = []
FC_hb = []
FC_hc = []
for i in range(len(ligand)):
    lig = ligand[i]
    #create a list with each ligands interactions as a sublist
    hblist.append(hbonds_identifier(protein, lig))
    hbs.append([(x[0],x[1]) for x in hblist[i]]) #retrieving only resid and resname for hbond

    hclist.append(hydrophobics_identifier(protein, lig))
    hcs.append([(x[0],x[1]) for x in hclist[i]])

    hbtot.extend(hbs[i])
    hctot.extend(hcs[i])

hbtot = set(hbtot)
hctot = set(hctot)
#binary arrays retrieving for interactions

for i in range(len(ligand)):
    FC_hb.append(feature_count(hbs[i], hbtot))
    FC_hc.append(feature_count(hcs[i], hctot))

#binary matrix for interactions
narray_hb = np.asarray(FC_hb)
narray_hc = np.asarray(FC_hc)

#taking in consideration the offset in residues numeration
res_hbtot = ["{1}_{0}".format(x[0], x[1]) for x in hbtot]
res_hctot = ["{1}_{0}".format(x[0], x[1]) for x in hctot]

#y axis label, the name of the ligand
lignames = []
for arg in sys.argv[2:]:
    lignames.append(os.path.basename(arg))


#plot the features presence matrix
#firstly for hbonds
fig, ax = plt.subplots()
fig.set_size_inches(12, 9)
ax.matshow(narray_hb, cmap=plt.cm.Blues)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(narray_hb.shape[1]), minor=False)
ax.set_yticks(np.arange(narray_hb.shape[0]), minor=False)
ax.invert_yaxis()

#labels
ax.set_xticklabels(res_hbtot, minor=False)
ax.set_yticklabels(lignames, minor=False)
plt.show()


##plot the features presence matrix
##firstly for hydrophobic
fig, ax = plt.subplots()
fig.set_size_inches(12, 9)
ax.matshow(narray_hc, cmap=plt.cm.Greens)
#
## put the major ticks at the middle of each cell
ax.set_xticks(np.arange(narray_hc.shape[1]), minor=False)
ax.set_yticks(np.arange(narray_hc.shape[0]), minor=False)
ax.invert_yaxis()
#
##labels
ax.set_xticklabels(res_hctot, minor=False)
ax.set_yticklabels(lignames, minor=False)
plt.show()



