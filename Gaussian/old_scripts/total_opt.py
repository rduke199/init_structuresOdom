def generate_gaussian_input(paramset=None,mol=None):
    route_parameters = paramset.route_parameters
    input_parameters = paramset.input_parameters
    charge = paramset.charge
    multiplicity = paramset.multiplicity
    functional = paramset.functional
    basis_set = paramset.basis_set

    ginput = GaussianInput(mol=mol,charge=charge, spin_multiplicity=multiplicity,
                           title=None, functional=functional, basis_set=basis_set,
                           route_parameters=route_parameters, input_parameters=input_parameters, dieze_tag="#")
    return ginput

def run_opt(calc_dir, functional, basis_set, calculation, omega=None, charge=0):
    identifier = mol_dir.split('/')[-1][:14]
    mol = Molecule.from_file(in_fn)
    mol = self.get("mol", )
    paramset = self["paramset"]
    gaussian_cmd = env_chk(self["g16_cmd"], fw_spec)


    # get the type of gaussian optimization calculation
    if self.get("type", ):
        if self["type"] == "coupling":
            subwdir = '{}/{}/opt/{}/{}'.format(calc_dir, self["type"], self["label"], self["subtype"])
        elif self["type"] == "reorganization":
            subwdir = '{}/{}/opt/{}'.format(calc_dir, self["type"], self["subtype"])
        else:
            subwdir = '{}/{}/opt'.format(calc_dir, self["type"], )
    else:
        subwdir = "{}/mol_opt".format(calc_dir)

    # create folders
    os.system('mkdir -p {}'.format(subwdir))
    os.chdir(subwdir)
    dir = os.getcwd()

    # set prefix for the file names
    if self.get("prefix"):
        prefix = self["prefix"]
    elif self.get("subtype" ,):
        prefix = "mol_" + self["subtype"]
    else:
        prefix = "mol_opt"

    file_com = prefix + ".com"
    file_chk = prefix + ".chk"
    file_log = prefix + ".log"

    # run gaussian optimize until normal termination is achieved. max runs is 5
    converged = False
    runs = 0

    while not converged and runs < 5:
        runs += 1

        # write input files for gaussian optimization calculation
        gauss_inp = generate_gaussian_input(paramset=paramset, mol=mol)
        gauss_inp.link0_parameters ={"%chk": file_chk, "%mem": "48GB" ,"%nprocshared": nprocs}
        gauss_inp.write_file(file_com, cart_coords=True)

        # run gaussian optimization
        return_code = subprocess.call(gaussian_cmd + " " + file_com, shell=True)
        gout = GaussianOutput(file_log)

        # check for normal termination
        if not gout.properly_terminated:
            mol = gout.final_structure
        else:
            # if netural molecule fails to converge rerun with last geometry
            if gauss_inp.charge == 0 and gauss_inp.spin_multiplicity == 1:
                if gout.frequencies[0][0]["frequency"] < 0:
                    pass
                else:
                    converged = True

                mol = gout.final_structure

            else:
                os.system("formchk " + file_chk)
                mol = gout.final_structure
                converged = True

    if not converged:
        raise RuntimeError("Structure not converged. calc_dir " + dir)

    return FWAction(update_spec={"gaussrun_dir": calc_dir, "identifier": identifier ,"mol_opt": mol})
