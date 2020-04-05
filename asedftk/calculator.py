from ase.units import Bohr, Hartree
from ase.calculators.calculator import Calculator, all_changes


class DFTK(Calculator):
    """Class for doing DFTK calculations."""

    # Parameters: model (some dummy way to set it first, ideally all of DFTK supported somehow)


    implemented_properties = ["energy", "free_energy", "forces"]
    ignored_changes = {'pbc'}  # DFTK always uses periodic BCs
    default_parameters = dict(pps="gth")
    discard_results_on_any_change = True

    def __init__(self, restart=None, ignore_bad_restart=False,
                 label="dftk", atoms=None, **kwargs):
        Calculator.__init__(self, restart=restart,
                            ignore_bad_restart=ignore_bad_restart, label=label,
                            atoms=atoms, command=None, **kwargs)
        self.setup_dftk(atoms=atoms)

    def setup_dftk(self, atoms=None):
        from julia import DFTK, Main

        # Set atoms consistently everywhere
        if atoms is None:
            if self.atoms is None:
                return None
            else:
                atoms = self.atoms
        if self.atoms is None:
            self.atoms = atoms

        # build Model:
        # TODO Temperature, spin polarisation, smearing
        Main.eval("""
            using DFTK

            function build_model(atoms, functional, family, core)
                lattice = load_lattice(atoms)

                atoms = [(ElementPsp(element.symbol, psp=load_psp(element.symbol,
                                                          functional=functional,
                                                          family=family,
                                                          core=Symbol(core)))
                  => positions) for (element, positions) in load_atoms(atoms)]

                if functional == "lda"
                    model_LDA(lattice, atoms)
                else
                    model_DFT(lattice, atoms, [:gga_x_pbe, :gga_c_pbe])
                end
            end
        """)
        self.model = Main.build_model(atoms, "lda", "hgh", "fullcore")

        Ecut = 15
        kgrid = [1, 1, 1]
        temperature = 0.0
        smearing = DFTK.Smearing.FermiDirac()
        spin_polarisation = DFTK.eval(":none")
        self.basis = DFTK.PlaneWaveBasis(self.model, Ecut, kgrid=kgrid)

    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def read(self, label):
        print("read called with ", label)
        pass

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        from julia import DFTK, Main

        Calculator.calculate(self, atoms, properties, system_changes)

        if self.atoms is None:
            raise CalculatorSetupError('An Atoms object must be provided to '
                                       'perform a calculation')
        self.setup_dftk(atoms)

        scfres = DFTK.self_consistent_field(self.basis, tol=1e-6)

        # TODO Unit conversion!
        self.results["energy"] = sum(Main.values(scfres.energies))

        # TODO Unit conversion!
        if 'forces' in properties:
            self.results["forces"] = DFTK.forces(scfres)
