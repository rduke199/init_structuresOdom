import os, re, json
from pymatgen.io.gaussian import GaussianOutput


class CollectGausData:
    """
    For a molecule directory, mol_dir, this collects the data from the geometry optimization, frequency calculation,
    and tddft calculation. It can also aggregate all this data
    """

    def __init__(self, mol_dir):
        self.mol_dir = mol_dir
        self.opt_file = os.path.join(self.mol_dir, 'opt.log')
        self.freq_file = os.path.join(self.mol_dir, 'freq.log')
        self.tddft_file = os.path.join(self.mol_dir, 'tddft/', 'tddft.log')
        self.wtuning_dir = os.path.join(self.mol_dir, 'wtuning/')

    @property
    def omega(self):
        output_path = os.path.join(self.wtuning_dir, 'output.log')
        try:
            with open(output_path, 'r') as fn:
                w_data = fn.readlines()[-2].split()[1]
            return "0{}".format(w_data.split('.')[1])
        except FileNotFoundError:
            return None

    @property
    def mol_name(self):
        """
        Get molecule name from file name
        """
        return '_'.join(self.mol_dir.split('/')[-2:-1])

    @staticmethod
    def mol_object(log_fn):
        """
        Get pymatgen object in preparation for further processing
        Returns:
            A pymatgen object
        """
        try:
            return GaussianOutput(log_fn)
        except FileNotFoundError:
            print("File does not exist: {}".format(log_fn))

    @staticmethod
    def runtime_from_log(logfile):
        """
        Collects runtime in core hours from a logfile
        """
        time_patt = re.compile(r"\d+\.d+|\d+")
        time_data = []
        with open(logfile, "r") as f:
            line = f.readline()
            while line != "":
                if re.match("Job cpu time", line.strip()):
                    time_data.extend(time_patt.findall(line))
                line = f.readline()
        if time_data:
            time_data = [float(time) for time in time_data]
            runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
        else:
            runtime = 0
        return round(runtime, 3)

    @staticmethod
    def get_tddft_excitations(logfile):
        """
        Read a excitation energies after a TD-DFT calculation.
        :param logfile, path to logfile
        Returns:
            A list: A list of tuple for each transition such as
                    [(energie (eV), lambda (nm), oscillatory strength), ... ]
        """

        float_patt = re.compile(r"\s*([+-]?\d+\.\d+)")
        state_patt = re.compile(r"[a-zA-Z]*let")
        transitions = {"Singlet": [],
                       "Triplet": []
                       }

        # read in file
        with open(logfile, "r") as f:
            line = f.readline()
            td = False
            while line != "":
                if re.search(r"^\sExcitation energies and oscillator strengths:", line):
                    td = True

                if td:
                    if re.search(r"^\sExcited State\s*\d", line):
                        val = [float(v) for v in float_patt.findall(line)]
                        try:
                            state = state_patt.findall(line)[0]
                        except:
                            state_val = val.pop(0)
                            if state_val == 1:
                                state = "Singlet"
                            else:
                                state = "Triplet"
                        # transitions.append(tuple(val[0:3]))
                        transitions[state].append(tuple(val))
                line = f.readline()
        return transitions

    @staticmethod
    def write_json(data, filename):
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)

    @property
    def opt_data(self):
        """
        Collect data from geometry optimization
        Returns:
            A dict: a dictionary containing optimization data
        """
        mol = self.mol_object(self.opt_file)
        if mol and mol.properly_terminated:
            structure = mol.final_structure
            num_electrons = mol.electrons[0]
            eigens = list(mol.eigenvalues.values())[0]
            homo = eigens[num_electrons - 1] * 27.2114
            lumo = eigens[num_electrons] * 27.2114
            homo_lumo_gap = lumo - homo

            json_data = {
                "charge": mol.charge,
                "structure": structure.to(fmt="xyz"),
                "energy": mol.final_energy * 27.2114,
                "homo": homo,
                "lumo": lumo,
                "homo_lumo_gap": homo_lumo_gap,
                "opt_runtime": self.runtime_from_log(self.opt_file)
            }
            return json_data
        else:
            print("Error. Calculation did not properly terminate for {}".format(self.mol_name))
            return {}

    @property
    def freq_data(self):
        """
        Collect data from frequency calculation
        Returns:
            A dict: a dictionary containing frequency calculation data
        """
        mol = self.mol_object(self.freq_file)
        if mol and mol.properly_terminated:
            try:
                frequencies = [f["frequency"] for f in mol.frequencies[0]]
                frequencies_dict = mol.frequencies
                # frequencies = len([f["frequency"] for f in mol.frequencies[0] if f["frequency"] < 0])
            except IndexError:
                frequencies = None
                frequencies_dict = None
            return {"frequencies": frequencies,
                    "frequencies_dict": frequencies_dict,
                    "freq_runtime": self.runtime_from_log(self.freq_file)}
        else:
            print("Error. Frequency calculation did not properly terminate for {}".format(self.mol_name))
            return {}

    @property
    def tddft_data(self):
        """
        Collect data from tddft calculation
        Returns:
            A dict: a dictionary containing tddft calculation data
            :return:
        """
        mol = self.mol_object(self.tddft_file)
        if mol and mol.properly_terminated:
            json_data = {
                "tddft_energy": mol.final_energy * 27.2114,
                "excitations": self.get_tddft_excitations(self.tddft_file),
                "tddft_runtime": self.runtime_from_log(self.tddft_file)
            }
            return json_data
        else:
            print("Error. Frequency calculation did not properly terminate for {}".format(self.mol_name))
            return {}

    @property
    def all_data(self):
        """
        Aggregates all calculation data
        Returns:
            A dict: a dictionary containing all data
        """
        return {**self.opt_data, **self.freq_data, **self.tddft_data}


def main():
    # Establish paths that will be used
    home = os.getcwd()
    opt_runs_path = os.path.join(home, 'opt_omega/')
    master_json_file = os.path.join(home, "master_data.json")

    # Iterate through mol directories and gather information
    master_data = {}
    mols = [m for m in os.listdir(opt_runs_path) if m.startswith('mol')]
    for mol in mols:
        mol_name = mol.split('.')[0]
        master_data[mol_name] = {}
        for charge in os.listdir(os.path.join(opt_runs_path, mol)):
            mol_dir_path = os.path.join(opt_runs_path, mol, charge)
            molecule = CollectGausData(mol_dir_path)
            # master_data[mol_name] = molecule.omega
            master_data[mol_name][charge] = molecule.all_data
            print("Data collected for {} {}".format(mol_name, charge))

    CollectGausData.write_json(master_data, master_json_file)


if __name__ == "__main__":
    main()
