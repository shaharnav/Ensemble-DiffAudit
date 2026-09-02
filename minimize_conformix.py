"""
Energy-minimize ConforMix receptor conformers to relieve local side-chain
steric distortion that breaks Meeko's bond-perception during receptor prep
(mk_prepare_receptor sees residues with side-chain atoms close enough in
space to look covalently bonded, even though backbone connectivity is clean).

Uses PDBFixer (adds missing atoms/hydrogens) + a short OpenMM restrained
minimization (restrain heavy atoms lightly so the fold isn't erased, just
local clashes relieved).
"""
import sys

from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation, HBonds, NoCutoff
from openmm import unit, CustomExternalForce, LangevinMiddleIntegrator


def minimize(in_pdb, out_pdb, restraint_k=100.0):
    fixer = PDBFixer(filename=in_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    forcefield = ForceField("amber14-all.xml")
    system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

    # Light positional restraint on heavy atoms so minimization relieves local
    # clashes without unfolding the ConforMix-predicted conformation.
    restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", restraint_k * unit.kilojoules_per_mole / unit.nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for atom in fixer.topology.atoms():
        if atom.element is not None and atom.element.symbol != "H":
            pos = fixer.positions[atom.index]
            restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
    system.addForce(restraint)

    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    sim = Simulation(fixer.topology, system, integrator)
    sim.context.setPositions(fixer.positions)

    state0 = sim.context.getState(getEnergy=True)
    sim.minimizeEnergy(maxIterations=2000)
    state1 = sim.context.getState(getEnergy=True, getPositions=True)

    print(f"  {in_pdb}: E {state0.getPotentialEnergy()} -> {state1.getPotentialEnergy()}")

    with open(out_pdb, "w") as f:
        PDBFile.writeFile(sim.topology, state1.getPositions(), f)


if __name__ == "__main__":
    in_pdb, out_pdb = sys.argv[1], sys.argv[2]
    minimize(in_pdb, out_pdb)
