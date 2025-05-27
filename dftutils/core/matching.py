from pymatgen.core.structure import Structure

class StructureMatcher:
    # Config
    # 1. type: 'all' or 'select' (all sites or some sites)
    # 2. selection: 'None' if match_type is 'all', otherwise a list of sites or a list of elements
    # 3. sort_first: sort the structures before matching

    def __init__(self, structure1, structure2, config):
        self.s1 = structure1.copy()
        self.s2 = structure2.copy()
        self.s1_backup = structure1.copy()
        self.s2_backup = structure2.copy()
        self.config = config

    def from_paths_and_config(structure1_path, structure2_path, config):
        structure1 = Structure.from_file(structure1_path)
        structure2 = Structure.from_file(structure2_path)
        return StructureMatcher(structure1, structure2, config)

    def match(self):
        if self.config['type'] == 'all':
            return self.match_all()
        elif self.config['type'] == 'select':
            if isinstance(self.config['selection'][0], int):
                return self.match_select(self.config['selection'])
            elif isinstance(self.config['selection'][0], str):
                sites = []
                for site in self.s1.sites:
                    if site.species_string in self.config['selection']:
                        sites.append(site)
                return self.match_select(sites)
            
    def match_all(self):
        if self.config['sort_first']:
            self.s1 = self.s1.sort()
            self.s2 = self.s2.sort()

        self.s1_backup = self.s1.copy()
        self.s2_backup = self.s2.copy()

        pairings = []
        for si, sa in enumerate(self.s1.sites):
            clj = -1
            mind = 10000000
            for sj, sb in enumerate(self.s2.sites):
                d = sa.distance(sb)
                if d <= mind:
                    mind = d
                    clj = sj
            pairings.append((si, clj))

        s3 = self.s2.copy()

        for pair in pairings:
            self.s2.sites[pair[0]] = s3.sites[pair[1]]

        return self.s1, self.s2
    
    def match_select(self, sites):
        self.s1_backup = self.s1.copy()
        self.s2_backup = self.s2.copy()

        pairings = []
        for si, sa in enumerate(self.s1.sites[sites]):
            clj = -1
            mind = 10000000
            for sj, sb in enumerate(self.s2.sites[sites]):
                d = sa.distance(sb)
                if d <= mind:
                    mind = d
                    clj = sj
            pairings.append((si, clj))

        s3 = self.s2.copy()

        for pair in pairings:
            self.s2.sites[pair[0]] = s3.sites[pair[1]]

        return self.s1, self.s2

