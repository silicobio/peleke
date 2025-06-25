#!/usr/bin/env python
import numpy as np
from configparser import ConfigParser


def dist(x1, y1, z1, x2, y2, z2):
	'''
	Calculates distance between two points in 3D
	arg1:   (float) x coordinate of atom 1 
	arg2:   (float) y coordinate of atom 1 
	arg3:   (float) z coordinate of atom 1 
	arg4:   (float) x coordinate of atom 2 
	arg5:   (float) y coordinate of atom 2
	arg6:   (float) z coordinate of atom 2
	return: (float) Distance between two atoms
	'''
	D = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
	return D 


def compatible(c1,c2):
	'''
	Determines if two classes are compatible according to CSU
	Table 1 from DOI: 10.1093/bioinformatics/15.4.327
	arg1:   (int) class of atom 1
	arg2:   (int) class of atom 2
	return: (boolean) True/False
	'''
	if c1 == '?' or c2 == '?':
		return True

	compatibleMatrix = np.array([[1,1,1,0,1,1,1,1],\
								 [1,0,1,0,1,1,1,0],\
								 [1,1,0,0,1,1,0,1],\
								 [0,0,0,1,1,1,1,1],\
								 [1,1,1,1,1,1,1,1],\
								 [1,1,1,1,1,1,1,1],\
								 [1,1,0,1,1,1,0,1],\
								 [1,0,1,1,1,1,1,0]])

	if compatibleMatrix[c1-1][c2-1] == 0:
		return False
	if compatibleMatrix[c1-1][c2-1] == 1:
		return True


def run_voro(points):
	"""
	Decides whether to run the voronoi calculation with python or C
	arg1:   (list of lists of strings) each list contains the coordinate of one atom
	return: (list of list of ints) each list contains the indices of neighboring atoms to the atom at that particular index  
	"""
	config = ConfigParser()
	config.read('intercaat/intercaat_config.ini')
	python_version = config.get('qvoronoi_path', 'run_python_version')
	if python_version == 'no':
		return voroC(points)
	else:
		return voroPython(points)


def voroC(points):
	'''
	Creates 3D voronoi diagram and returns indices of neighboring atoms 
	For example: if neighbors = [[1,2],[0,2,3],[0,3],[1]], atom 0 would have atoms 1 and 2 as neighbors
	arg1:   (list of lists of strings) each list contains the coordinate of one atom
	return: (list of list of ints) each list contains the indices of neighboring atoms to the atom at that particular index  
	'''
	import random
	import subprocess

	# Puts arg1(points) in format qhull can read 
	newPoints = ['3', str(len(points))]
	for location in points:
		newPoints.append(location[0] + ' ' + location[1] + ' ' + location[2])

	# generate file with random tag to not interfere with others who may be using the program at the same time
	randomTag = str(random.randint(1,10000))
	textFile = 'save' + randomTag + '.txt'
	newFile = open(textFile, 'w')

	for line in newPoints:
		newFile.write(line)
		newFile.write('\n')

	newFile.close()
	config = ConfigParser()
	config.read('intercaat_config.ini')
	
	qhullPath = config.get('qvoronoi_path', 'qvoronoi_bin')
	qvoronoi  = config.get('qvoronoi_path', 'executable_name')
	debug  = config.get('qvoronoi_path', 'debug_qvoronoi')
	
	vorFi = subprocess.check_output('cat ' + textFile + ' | ' + qhullPath + qvoronoi, shell =True, text= True)
	if debug == 'yes':
		pass
	else:	
		subprocess.check_output('rm ' + textFile, shell =True, text= True)
	
	# Process qhull output to obtain only contacting atoms
	vorFi = vorFi.split('\n')
	vorFi.remove(vorFi[0])
	vorFi.remove(vorFi[-1])
	contacts = []
	for line in vorFi:
		LineList = line.split()
		contacts.append([int(LineList[1]),int(LineList[2])])

	# creates a list of lists of the indices of neighboring atoms to the atom at that particular index
	contacts = np.array(contacts)
	count1 = 0
	count2 = 0
	hold = []
	neighbors = []
	while count1 < len(points):
		while count2 < len(contacts):
			if contacts[count2,0] == count1:
				hold.append(contacts[count2,1])
			if contacts[count2,1] == count1:
				hold.append(contacts[count2,0])
			count2 += 1 
		hold = sorted(hold)
		neighbors.append(hold)
		hold = []
		count2 = 0
		count1 += 1
	return neighbors


def voroPython(points):
	from scipy.spatial import Voronoi
	'''
	Creates 3D voronoi diagram and returns indices of neighboring atoms 
	For example: if neighbors = [[1,2],[0,2,3],[0,3],[1]], atom 0 would have atoms 1 and 2 as neighbors
	arg1:   (numpy list of lists of strings) each list contains the coordinate of one atom
	return: (list of list of ints) each list contains the indices of neighboring atoms to the atom at that particular index  
	'''
	v = Voronoi(np.array(points))
	contacts = v.ridge_points

	count1 = 0
	count2 = 0
	hold = []
	neighbors = []
	while count1 < len(points):
		while count2 < len(contacts):
			if contacts[count2,0] == count1:
				hold.append(contacts[count2,1])
			if contacts[count2,1] == count1:
				hold.append(contacts[count2,0])
			count2 += 1 
		hold = sorted(hold)
		neighbors.append(hold)
		hold = []
		count2 = 0
		count1 += 1
	return neighbors


import os
import gzip

def parse(pdb_file, chains, path):
    """
    Parses a PDB or PDB.GZ file and extracts atom-level information 
    for the specified chains.

    Parameters:
        pdb_file (str): Filename of the PDB (e.g., "7vux.pdb" or "7vux.pdb.gz")
        chains (list): List of chain IDs to include (e.g., ['A', 'B'])
        path (str): Path to the directory containing the file (default: current dir)

    Returns:
        list: Parsed atom lines, where each line is a list of atom-level fields.
    """
    pdb_path = os.path.join(path, pdb_file)

    # Detect file type
    if pdb_file.endswith('.gz'):
        open_fn = gzip.open
        mode = 'rt'  # read text from gzip
    else:
        open_fn = open
        mode = 'r'

    pdb_lines = []
    with open_fn(pdb_path, mode) as f:
        for line in f:
            # Only process ATOM or HETATM lines (excluding hydrogens and water)
            if (line.startswith('ATOM') or line.startswith('HETATM')) and ' H ' not in line and 'HOH' not in line:
                chain_id = line[21].strip()
                if chain_id in chains:
                    try:
                        parsed_line = [
                            line[0:6].strip(),              # 0: Record type
                            int(line[6:11].strip()),         # 1: Atom serial number
                            line[12:16].strip(),             # 2: Atom name
                            line[17:20].strip(),             # 3: Residue name
                            chain_id,                        # 4: Chain ID
                            line[22:26].strip(),             # 5: Residue sequence number
                            float(line[30:38].strip()),      # 6: X
                            float(line[38:46].strip()),      # 7: Y
                            float(line[46:54].strip()),      # 8: Z
                            float(line[54:60].strip()),      # 9: Occupancy
                            float(line[60:66].strip()),      # 10: Temp factor (B-factor)
                            line[76:78].strip()              # 11: Element symbol
                        ]
                        pdb_lines.append(parsed_line)
                    except Exception as e:
                        print(f"Skipping malformed line: {line.strip()}\nReason: {e}")
                        continue

    return pdb_lines



def planar(z):
	'''
	Determines if a residue contains a ring and if that ring is planar
	If a residue has more than 1 ring, function assumes both rings are planar 
	arg1:   (np matrix of 3 floats) xyz coordinates of every atoms in a single residue
	return: (list of ints) list of the indices of atoms that are planar in a residue
	return: False if the residue has a ring but it is not planar (ex. proline)
	'''
	from math import ceil

	count1 = 0
	count2 = 0
	distMatrix = np.zeros((len(z[:,1]), len(z[:,1])))
	holdDistMatrix1 = np.zeros((len(z[:,1]), len(z[:,1])))
	holdDistMatrix2 = np.zeros((len(z[:,1]), len(z[:,1])))

	# create distance matrix (distance of residue atom to every other residue atom)
	# only keep distances less than 1.6 angstroms: covalent bond.
	for i in z:
		while count1 < len(z[:,1]):
			hold = dist(i[0], i[1], i[2], z[count1,0], z[count1,1], z[count1,2])
			if hold < 1.6:
				distMatrix[count1, count2] = hold
			count1 += 1
		count1 = 0
		count2 += 1

	holdDistMatrix1 = distMatrix

	moreThanOneRing = False
	# If a residue has more than 1 ring the program assumes both are planar
	if np.count_nonzero(np.tril(distMatrix)) > len(z):
		moreThanOneRing = True

	# If the below statement is True, the residue does not contain any rings
	if np.count_nonzero(np.tril(distMatrix)) < len(z):
		return False

	# iterativley remove atoms unitl all remaining atoms have at least two neighbots 
	# determines which atoms are in a ring
	count1 = 0
	count2 = 0
	while np.array_equal(holdDistMatrix1, holdDistMatrix2) == False:
		holdDistMatrix2 = np.copy(holdDistMatrix1)
		while count1 < len(holdDistMatrix1):
			if np.count_nonzero(holdDistMatrix1[count1,:]) == 1:
				while count2 < len(holdDistMatrix1):
					if holdDistMatrix1[count1,count2] != 0:
						holdDistMatrix1[count1,count2] = 0
						holdDistMatrix1[count2,count1] = 0
					count2 += 1
			count1 += 1
			count2 = 0
		count1 = 0

	# the indices of the matrix represent each atom in a residue
	# only need lower left half of matrix since top right half is redundant 
	bondMatrix = np.tril(holdDistMatrix1)

	# Obtain indices of atoms in the ring
	count1 = 0
	count2 = 0
	indices = []
	while count1 < len(bondMatrix):
		while count2 < len(bondMatrix):
			if bondMatrix[count1][count2] > 0:
				indices.append([count1,count2])
			count2 += 1
		count2 = 0
		count1 += 1

	# Creates lists of 4 atoms to check dihedral angle of 
	count1 = 0
	hold1 = 0
	hold2 = 0
	dihedralAtoms = []
	while count1 < ceil(len(indices)/2):
		for index in indices:
			if index == indices[count1]:
				continue
			elif indices[count1][0] in index or indices[count1][1] in index:
				if index[0] == indices[count1][0]:
					hold1 = index[1] 
				if index[1] == indices[count1][0]:
					hold1 = index[0]
				if index[0] == indices[count1][1]:
					hold2 = index[1] 
				if index[1] == indices[count1][1]:
					hold2 = index[0]
		dihedralAtoms.append([hold1,indices[count1][0],indices[count1][1],hold2])
		count1 += 1

	# Checks the dihedral angle and appends angle to dihedralAngle
	count1 = 0
	dihedralAngle = []
	while count1 < len(dihedralAtoms):
		coordinatesHold = np.array([[z[dihedralAtoms[count1][0]][0], z[dihedralAtoms[count1][0]][1], z[dihedralAtoms[count1][0]][2]], \
									[z[dihedralAtoms[count1][1]][0], z[dihedralAtoms[count1][1]][1], z[dihedralAtoms[count1][1]][2]], \
									[z[dihedralAtoms[count1][2]][0], z[dihedralAtoms[count1][2]][1], z[dihedralAtoms[count1][2]][2]], \
									[z[dihedralAtoms[count1][3]][0], z[dihedralAtoms[count1][3]][1], z[dihedralAtoms[count1][3]][2]]])
		dihedralAngle.append(dihe(coordinatesHold))
		count1 += 1

	# if the average of the dihedral angles is less that 3 returns list of planar atoms 
	planarAtoms = []
	if sum(dihedralAngle)/len(dihedralAngle) < 3 or moreThanOneRing == True:
		for atomIndex in indices:
			planarAtoms.append(atomIndex[0])
			planarAtoms.append(atomIndex[1])
		return list(set(planarAtoms))
	else:
		return False

def inter(atom1, atom2, solv):
	'''
	Takes two atoms are returns the distance between them and the minimum
	distance required to fit a solvent molecule in between the two atoms.
	Arbitrary default van der waal radius = 1.8
	Radii taken from doi: 10.1021/j100785a001
	arg1:   (list of string and 3 floats) atom1 [atom, x coordinate, y coordinate, z coordinate]
	arg2:   (list of string and 3 flaots) atom2 [atom, x coordinate, y coordinate, z coordinate]
	arg3:   (float) (optional) solvent van der waal radius
	arg1 or arg2 example: ['C' 20.32 19.45 17.53]
	return: (float) distance between two atoms
	return: (float) minimum distance required to fit a sovlent molecule between the two atoms
	'''
	rwS = 1.8
	rwO = 1.52
	rwN = 1.55
	rwC = 1.8
	rwI = 2.01
	rwB = 1.84
	rwF = 1.4
	rwCl = 1.73
	rwSolvent = float(solv)
	atom1rw = 1.8
	atom2rw = 1.8

	if atom1[0][0:1] == 'C':
		atom1rw = rwC
	if atom1[0][0:1] == 'N':
		atom1rw = rwN
	if atom1[0][0:1] == 'S':
		atom1rw = rwS
	if atom1[0][0:1] == 'O':
		atom1rw = rwO
	if atom1[0][0:1] == 'B':
		atom1rw = rwB
	if atom1[0][0:1] == 'F':
		atom1rw = rwF
	if atom1[0][0:1] == 'I':
		atom1rw = rwI
	if atom1[0][0:2] == 'CL' or atom2[0][0:2] == 'Cl':
		atom1rw = rwCl

	if atom2[0][0:1] == 'C':
		atom2rw = rwC
	if atom2[0][0:1] == 'N':
		atom2rw = rwN
	if atom2[0][0:1] == 'S':
		atom2rw = rwS
	if atom2[0][0:1] == 'O':
		atom2rw = rwO
	if atom2[0][0:1] == 'B':
		atom2rw = rwB
	if atom2[0][0:1] == 'F':
		atom2rw = rwF
	if atom2[0][0:1] == 'I':
		atom2rw = rwI
	if atom2[0][0:2] == 'CL' or atom2[0][0:2] == 'Cl':
		atom2rw = rwCl

	atom1atom2Dist = dist(atom1[1],atom1[2],atom1[3],atom2[1],atom2[2],atom2[3])
	vanDerWaalDist = atom1rw + atom2rw + 2*rwSolvent
	return(atom1atom2Dist, vanDerWaalDist)

def dihe(p):
	'''
	Obtains the dihedral angle from the coordinates of 4 atoms
	arg1:   (np matrix of floats) contains the X, Y, and Z coordinates of 4 atoms
	return: The dihedral angle between 4 points
	'''
	# obtain 3 vectors from 4 points
	vectors = p[:-1] - p[1:]
	vectors[0] *= -1
	# calculate normal to two vectors 
	normalVec = np.zeros((2,3))
	normalVec[0,:] = np.cross(vectors[0],vectors[1])
	normalVec[1,:] = np.cross(vectors[2],vectors[1])
	# normalize vectors
	normalVec /= np.sqrt(np.einsum('...i,...i', normalVec, normalVec)).reshape(-1,1)
	# return dihedrall angle
	return np.degrees(np.arccos( normalVec[0].dot(normalVec[1])))


def aClass(coordinates, atom, noPop):
	""" 
	Assigns each atom a class. Class definitions were copied from CSU
	Always assume nitrogen is class III. it can sometimes be class I but it is unclear when
	DOI: 10.1093/bioinformatics/15.4.327 
	arg1:   (list of list of floats) list of coordinates of each atom in a residue. Sometimes the coordinates of the atom 
	in the next residue. For example when assigning class to amino acids, you need to consider the N in the next residue.
	arg2:   (list of strings) list of atoms types for each atom in a residue
	arg3:	(boolean) Decides whether to return the atom in the next residue or not
	return: (dict) returns a dictionary with the atom index and the class it belongs to for
				   a single residue 
	return: (string) '?' if unsure of class type.    
	NUMBER  ATOM CLASS
	I       Hydrophilic 
	II      Acceptor 
	III     Donor 
	IV      Hydrophobic 
	V       Aromatic 
	VI      Neutral 
	VII     Neutral-donor 
	VIII    Neutral-acceptor 
	"""
	count1 = 0
	count2 = 0
	distMatrix = np.zeros((len(coordinates[:,1]), len(coordinates[:,1])))

	# create distance matrix (distance of residue atom to every other residue atom)
	# only keep distances less than 1.6 angstroms: covalent bond.
	for i in coordinates:
		while count1 < len(coordinates[:,1]):
			hold = dist(i[0], i[1], i[2], coordinates[count1,0], coordinates[count1,1], coordinates[count1,2])
			if hold < 1.7:
				distMatrix[count1, count2] = hold
			count1 += 1
		count1 = 0
		count2 += 1

	# Determines if residue is planar
	pla = False
	if np.count_nonzero(np.tril(distMatrix)) >= len(coordinates):
		pla = planar(coordinates)

	# Assigns atom class to oxygen atoms and as a default assings all carbons to class 4
	count1 = 0
	count2 = 0
	atomClasses = {}
	while count1 < len(atom):
		if atom[count1][0:1] == 'O':
			atomClasses[count1] = '1'
			if np.count_nonzero(distMatrix[count1,:]) == 1:
				while count2 < len(distMatrix):
					if distMatrix[count1,count2] != 0 and distMatrix[count1,count2] < 1.3:
						atomClasses[count1] = '2'
					count2 += 1
			if np.count_nonzero(distMatrix[count1,:]) > 1:
				atomClasses[count1] = '2'
			# The program incorrectly assigns atom O3' to class 1. This atom is specifically found in DNA and RNA
			# It assigns O3' to class 1 because the query residue does not contain the phosphate from the next residue
			# I could not think of a simple way to do this.
			# Can be fixed by upgrading function planar to be able to look at more than 1 ring at once on a single residue
			if atom[count1] == 'O3\'':
				atomClasses[count1] = '2'
		if atom[count1][0:1] == 'C':
			atomClasses[count1] = '4'
		count1 += 1
		count2 = 0

	# Assigns atom class to sulfur, potassium, fluorine atoms
	# Assigns atom class to chlorine, bromine, and iodine atoms
	# Assigns atom class to nitrogen. 
	count1 = 0
	while count1 < len(atom):
		if atom[count1][0:1] == 'S' or atom[count1][0:1] == 'P' or atom[count1][0:1] == 'F':
			atomClasses[count1] = '6'
		if atom[count1][0:2] == 'Cl' or atom[count1][0:2] == 'Br' or atom[count1][0:1] == 'I' or\
		   atom[count1][0:2] == 'CL' or atom[count1][0:2] == 'BR':
			atomClasses[count1] = '4'
		if atom[count1][0:1] == 'N':
			atomClasses[count1] = '3'
		count1 += 1

	# Assings '?' to all atoms not givin a class
	count1 = 0
	if len(atomClasses) < len(atom):
		while count1 < len(atom):
			try:
				atomClasses[count1]
			except:
				atomClasses[count1] = '?'
			count1 += 1

	# Assigns atom class to carbon
	count1 = 0
	count2 = 0
	holdClass= []
	while count1 < len(atom):
		if atom[count1][0:1] == 'C' and (atom[count1][0:2] != 'CL' or atom[count1][0:2] != 'Cl'):
			while count2 < len(distMatrix):
				if distMatrix[count1,count2] != 0:
					holdClass.append(atomClasses[count2])
				count2 += 1
			if holdClass.count('3') == 1:
				atomClasses[count1] = '7'
			if holdClass.count('2') == 1:
				atomClasses[count1] = '8'
			if holdClass.count('1') == 1:
				atomClasses[count1] = '6'
			if holdClass.count('2') + holdClass.count('3') > 1:
				atomClasses[count1] = '6'
		holdClass = []
		count1 += 1
		count2 = 0

	# Assings atom class to planar carbon atoms
	count1 = 0
	if pla != False:
		while count1 < len(pla):
			if atom[pla[count1]][0:1] == 'C':
				atomClasses[pla[count1]] = '5'
			count1 += 1

	# Decides whether to return the atom in the next residue or not
	if noPop == False:
		atomClasses.pop(len(atomClasses) - 1)

	return atomClasses	


def appendAtomClasses(pdb):
	"""
	Calls aClass to append all the atoms classes of each residue to a list (pdbAtomClass)
	arg1:   (list of lists of strings) parsed pdb file being analyzed
	return: (list) The atom classes of all atoms in a pdb file. Should be same length as pdb file
	"""
	pdb.append(['','','','','','','','',1,1,1])
	count1 = 1
	count2 = 0
	dictHold = {}
	pdbAtomClass = []
	classCoordinates = [[float(pdb[0][8]), float(pdb[0][9]), float(pdb[0][10])]]
	noPop = False
	atm = [pdb[0][2]]
	while count1 < len(pdb):
		# Second half of if statent is needed because the function planar can only handle a residue with one ring
		# Since DNA/RNA has at least two rings, the residue must be seperated into pieces, the backbone, and the nitrogenous base
		if pdb[count1 - 1][6] == pdb[count1][6] and \
		(('\'' in pdb[count1 - 1][2] and '\'' in pdb[count1][2]) or ('\'' not in pdb[count1 - 1][2] and '\'' not in pdb[count1][2]) \
		or ('\'' not in pdb[count1 - 1][2] and '\'' in pdb[count1][2])):
			classCoordinates.append([float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])])
			atm.append(pdb[count1][2])
		else:
			# For DNA, does not include the phosphate in the next residue with the nitrogenous base.
			if pdb[count1][2] != 'P':
				classCoordinates.append([float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])])
				atm.append(pdb[count1][2])
			else:
				noPop = True
			dictHold = aClass(np.array(classCoordinates),atm, noPop)
			while count2 < len(dictHold):
				try:
					pdbAtomClass.append(int(dictHold[count2]))
				except:
					pdbAtomClass.append(dictHold[count2])
				count2 += 1
			dictHold = {}
			classCoordinates = [[float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]]
			atm = [pdb[count1][2]]
			noPop = False
		count1 += 1
		count2 = 0
	pdb.pop(len(pdb)-1)
	return pdbAtomClass


def filterMatch(match, pdb, arg2, arg4, arg5):
	"""
	Filters the list generated above (match) based on arg2 and arg4
	Also displays interactions matrix based on arg5
	arg1:   (list of list of strings) All atomic interactions 
	arg2:   (list of list of strings) parsed pdb file being analyzed
	arg3:   (list of string(s)) Query chains and optional residues
	arg4:   (int) minimum number of atomic interactions
	arg5:   (string) option to display interaction matrix. 'yes' or 'no'
	return: (list of list of strings) filterd match
	"""
	# parse match (the list generated in the last for loop)
	l = []
	hold = []
	res = []
	for line in match:
		l = line.split(' ')
		for element in l:
			if len(element) > 0:
				hold.append(element)
		res.append(hold)
		hold = []

	# need to add empty list to the end of res for next loop to work
	res.append(['','','','','','','','','','','','','',''])
	res = np.array(res)

	# Sorts arg2. Need to turn into int to sort. Also remove chain from arg2
	count4 = 0
	# if arg2 only contains chain identifier
	if len(arg2) == 1:
		for AA in pdb:
			if AA[5] == arg2[0] and int(AA[6]) not in arg2:
				arg2.append(int(AA[6]))
		arg2.remove(arg2[0])
	# if arg2 contains the chain identifier and specific residues
	else:
		arg2Hold = []
		for AA in arg2[1:]:
			arg2Hold.append(int(AA)) 
		arg2 = arg2Hold
	arg2 = sorted(arg2)
	# turns values in arg2 back to strings
	while count4 < len(arg2):
		arg2[count4] = str(arg2[count4])
		count4 += 1

	# obtain matrix of # of interactions/residue
	# first column: residue number, second column: # of interactions
	# count1 loops through while loop, count2 is an index for arg2, count3 counts # of interactions
	count1 = 1
	count2 = 0
	count3 = 1
	interactions = []
	while count1 < len(res):
		if res[count1 -1][1] != res[count1][1]:
			interactions.append([int(res[count1 - 1][1]), count3])
			count2 += 1
			count3 = 1
		else:
			count3 +=1
		count1 += 1
	interactions = np.array(interactions)

	# Creates new list of matches that have min # of interactions
	count1 = 0
	count2 = 0
	newMatch = []
	while count1 < len(res):
		while count2 < len(interactions):
			if res[count1][1] == str(int(interactions[count2][0])):
				if interactions[count2][1] >= arg4:
					newMatch.append(match[count1])
			count2 += 1
		count2 = 0
		count1 += 1
	
	# Prints a matrix of interactions/residue
	newInteractions = []
	count1 = 0
	if arg5 == 'yes':
		for i in interactions:
			if i[1] >= arg4:
				newInteractions.append(list(i))
			count1 += 1
		newMatch.append('--------------------------------')
		newInteractionRes = []
		count1 = 0
		while count1 < (len(newMatch)-1):
			if newMatch[count1][5:9] != newMatch[count1+1][5:9]:
				newInteractionRes.append(newMatch[count1][0:3])
			count1 += 1
		count1 = 0
		print('Res #   Interactions')
		while count1 < len(newInteractions):
			print('{0} {1:<5}  {2:<4}'.format(newInteractionRes[count1], newInteractions[count1][0], newInteractions[count1][1]))
			count1 += 1
		newMatch.pop(len(newMatch)-1)

	return newMatch
