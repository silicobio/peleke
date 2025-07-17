#!/usr/bin/env python

import sys
import argparse
import numpy as np 
import intercaat_functions as icaat

descr = '''
		Display all atomic interactions

		Required arguments:
			Program name: intercaat.py
			PDB file name (must be in same directory as intercaat.py): ex. -pdb 1cph.pdb
			Query chain (optional: specify query chain residues): ex1. -qc B   ex2. -qc B,20,21,25,27
			Interacting chain(s): ex1. -ic C ex2. -ic C,D

		DOI: 10.1093/bioinformatics/15.4.327 
		NUMBER  ATOM CLASS
		  1     Hydrophilic 
		  2     Acceptor 
		  3     Donor 
		  4     Hydrophobic 
		  5     Aromatic 
		  6     Neutral 
		  7     Neutral-donor 
		  8     Neutral-acceptor

		'''
parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description = descr)

parser.add_argument("-pdb", "--PDBFileName", help = "Name of saved PDB file.", required = True)
parser.add_argument("-qc", "--QueryChain", help = "Chain indicator of query chain. If you do not \
					want to look at the entire chain, you can specify particular residues", required = True)
parser.add_argument("-ic", "--InteractingChains", help = "Input example: -ic B,C. \
				   Chains interacting with the query chain.", required = True)
parser.add_argument( "-mi", "--MinInteractions", help = "Input example: -mi 3. Minimum number of interactions \
					required to be considered as part of the interface. default = 4", default = 4)
parser.add_argument("-di", "--DisplayInteractions", help = "Input: either -di yes or -di no. \
					If yes, displays interactions matrix. default = yes", default = 'yes')
parser.add_argument("-cc", "--ClassCompatibility", help = "Input: either -cc yes or -cc -no. If no, consider all \
					interactions. If yes, consider only class compatible interactions. default = yes", default = 'yes')
parser.add_argument("-sr", "--SolventRadius", help = "Input example: -sr 1.6. Solvent molecule radius. Decreasing this \
					value would require atoms to be closer together to interact. default = 1.4", default = '1.4')
parser.add_argument("-vi", "--VoronoiInclude", help = "Input example: -ec A,D. To increase the accuracy of this program, \
				   you can include other chains in the voronoi calculation. \
				   A consequence of this is a longer runtime. default = []", default = [])
parser.add_argument("-fp", "--FilePath", help = "Input example: /home/steven/. File path of PDB file. \
				   The default path is your current directory. default = ./", default = "./")

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()
arg2 = args.QueryChain.split(',')
arg3 = args.InteractingChains.split(',')
arg4 = int(args.MinInteractions)
arg5 = args.DisplayInteractions
arg6 = args.ClassCompatibility
arg7 = float(args.SolventRadius)
arg8 = args.VoronoiInclude
arg9 = args.FilePath
if arg8 != []:
	arg8 = arg8.split(',')

try:
    pdb = icaat.parse(args.PDBFileName, (arg2 + arg3 + arg8), arg9)
except FileNotFoundError:
    sys.exit(f'{pdb_filename} was not found')

coordinates = []
match = []

# Parse PDB file into list of PDB lines containing atoms and hetatms excluding hydrogen and water.
for line in pdb:
	coordinates.append([line[8], line[9], line[10]])

# Creates 3D voronoi diagram and returns indices of neighboring atoms
contacts = icaat.run_voro(coordinates)

# Creates a list (pdbAtomClass) that contains classes for all atoms analayzed 
pdbAtomClass = icaat.appendAtomClasses(pdb)

# Creates list of all atomic interactions
count1 = 0
count2 = 0
for buddies in contacts:
	while count2 < len(buddies):
		# XYZ = atom, X, Y, Z 
		XYZ1 = [pdb[count1][2][0:2], float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]
		XYZ2 = [pdb[buddies[count2]][2][0:2], \
		float(pdb[buddies[count2]][8]), float(pdb[buddies[count2]][9]), float(pdb[buddies[count2]][10])]
		# Returns the distance between two atoms and min distance reqeuired for solvent molecule
		Ad, Vd = icaat.inter(XYZ1,XYZ2,arg7)
		# Finds class of atom1 and atom2
		class1 = pdbAtomClass[count1]
		class2 = pdbAtomClass[buddies[count2]]
		# atom = residue, residue #, chain, atom
		atom1 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[count1][4], pdb[count1][6], pdb[count1][5], pdb[count1][2])
		atom2 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[buddies[count2]][4], pdb[buddies[count2]][6],\
				 pdb[buddies[count2]][5], pdb[buddies[count2]][2])
		Line = '{0} | {1} | {2:<4} |    {3}   {4}'.format(atom1, atom2, str(round(Ad,2)), str(class1), str(class2))
		# Only appends if classes are compatible or the user inputs that class compatibility does not matter
		# if class compatibility is unknown, the atomic interaction will be shown 
		if icaat.compatible(class1,class2) == True or arg6 == 'no':
			# line 1: appends if distance between atoms is within bound and accounts for occupancy < 1
			# line 2: appends only for specific chain
			if Ad < Vd and any((atom1 + atom2) in sub for sub in match) == False \
			and pdb[count1][5] == arg2[0]:
				# appends all residues from specified chain
				if len(arg2) == 1:
					# appends only against queried neighbor chains
					if pdb[buddies[count2]][5] in arg3:
						match.append(Line)
				#appends only specified residues from specified chain
				elif pdb[count1][6] in arg2:
					# appends only against queried neighbor chains
					if pdb[buddies[count2]][5] in arg3:
						match.append(Line)
		count2 += 1
	count2 = 0
	count1 += 1

# Filters the list generated above (match) based on arg2 and arg4
# Also displays interactions matrix based on arg5
newMatch = icaat.filterMatch(match, pdb, arg2, arg4, arg5)

print('  Query Chain    |Interacting Chains| Dist | AtomClasses')
for line in newMatch:
	print(line)
