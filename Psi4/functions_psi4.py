# class Psi4Output:
#     def __init__(self, filename):
#         """
#         Args:
#             filename: Filename of Psi4 output file.
#         """
#         self.filename = filename
#         self._parse(filename)
#
#     @property
#     def final_energy(self):
#         """
#         :return: Final energy in Gaussian output.
#         """
#         return self.energies[-1]
#
#     @property
#     def final_structure(self):
#         """
#         :return: Final structure in Gaussian output.
#         """
#         return self.structures[-1]
#
#     def parse_dat(out_file):

def get_final_geometry(filename):
    structures = []
    with open(filename,'r'):
        data = filename.readlines()
    for line in data:
        if re.search('Geometry (in Angstrom)'):
            geom_start_line = line
            count = 4

            components = row.split()