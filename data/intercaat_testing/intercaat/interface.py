from . import intercaat_functions as icaat


class InterfaceAnalyzer:
    def __init__(self, pdb_file, query_chain, interact_chains, path="./", min_contacts=4, solvent_radius=1.4):
        self.pdb_file = pdb_file
        self.query_chain = query_chain
        self.interact_chains = interact_chains
        self.path = path
        self.min_contacts = min_contacts
        self.solvent_radius = solvent_radius

    def get_interface_residues(self):
        chains = [self.query_chain] + self.interact_chains
        pdb = icaat.parse(self.pdb_file, chains, self.path)

        coordinates = [[line[8], line[9], line[10]] for line in pdb]
        contacts = icaat.run_voro(coordinates)
        atom_classes = icaat.appendAtomClasses(pdb)

        matches = []
        for i, buddies in enumerate(contacts):
            for b in buddies:
                Ad, Vd = icaat.inter(
                    [pdb[i][2][:2], float(pdb[i][8]), float(pdb[i][9]), float(pdb[i][10])],
                    [pdb[b][2][:2], float(pdb[b][8]), float(pdb[b][9]), float(pdb[b][10])],
                    self.solvent_radius  # or whatever value is appropriate for 'solv'
                )
                class1 = atom_classes[i]
                class2 = atom_classes[b]
                if icaat.compatible(class1, class2):
                    if Ad < Vd and pdb[i][5] == self.query_chain:
                        if pdb[b][5] in self.interact_chains:
                            matches.append('{0:<3} {1:>5} {2} {3:<4} | {4:<3} {5:>5} {6} {7:<4} | {8:<4} | {9} {10}'.format(
                                pdb[i][4], pdb[i][6], pdb[i][5], pdb[i][2],
                                pdb[b][4], pdb[b][6], pdb[b][5], pdb[b][2],
                                str(round(Ad, 2)), str(class1), str(class2)
                            ))

        return icaat.filterMatch(matches, pdb, [self.query_chain], self.min_contacts, 'no')
