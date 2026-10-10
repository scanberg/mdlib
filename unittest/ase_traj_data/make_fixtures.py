"""Regenerate the small ASE .traj fixtures (requires ASE and NumPy)."""

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import Trajectory, ulm


here = Path(__file__).parent

atoms = Atoms(
    "H2O",
    positions=[(0, 0, 0), (0.9, 0, 0), (0, 0.8, 0)],
    cell=[10, 10, 10],
    pbc=True,
)
with Trajectory(here / "ase_fixed.traj", "w") as trajectory:
    for step in range(3):
        atoms.positions[:, 0] += 0.1 * step
        atoms.set_cell([10 + step, 10, 10])
        if step == 2:
            atoms.set_cell([[12, 0, 0], [1.5, 10, 0], [-0.5, 0.75, 10]])
        atoms.info["time_ps"] = step * 0.25
        trajectory.write(atoms)

with ulm.open(here / "ase_float32.traj", "w", tag="ASE-Trajectory") as writer:
    writer.write(
        version=1,
        ase_version="fixture",
        pbc=[False] * 3,
        numbers=np.array([0, 8], dtype=np.int32),
        positions=np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32),
        cell=[[0.0] * 3] * 3,
    )
    writer.sync()
    writer.write(
        positions=np.array([[1, 2, 3], [7, 8, 9]], dtype=np.float32),
        cell=[[0.0] * 3] * 3,
    )

changed = Atoms("H", positions=[(0, 0, 0)], cell=[10] * 3, pbc=True)
with Trajectory(here / "ase_changing_atoms.traj", "w") as trajectory:
    trajectory.write(changed)
    changed += Atoms("H", positions=[(1, 0, 0)])
    trajectory.write(changed)

rotated = Atoms("H", positions=[(1, 2, 3)], cell=[[9, 1, 0], [-1, 9, 0], [0, 0, 10]], pbc=True)
with Trajectory(here / "ase_rotated_cell.traj", "w") as trajectory:
    trajectory.write(rotated)

partial_pbc = Atoms("H", positions=[(1, 2, 3)], cell=[10, 10, 10], pbc=[True, False, True])
with Trajectory(here / "ase_partial_pbc.traj", "w") as trajectory:
    trajectory.write(partial_pbc)

# One frame is a structure, not a trajectory
with Trajectory(here / "ase_single_frame.traj", "w") as trajectory:
    trajectory.write(Atoms("H2O", positions=[(0, 0, 0), (0.9, 0, 0), (0, 0.8, 0)], cell=[10, 10, 10], pbc=True))

# Not periodic: the box ASE keeps around the molecule is not a cell to wrap in
molecule = Atoms("H2O", positions=[(0, 0, 0), (0.9, 0, 0), (0, 0.8, 0)], pbc=False)
molecule.center(vacuum=4.0)
with Trajectory(here / "ase_vacuum_box.traj", "w") as trajectory:
    for step in range(2):
        molecule.positions[0, 2] += 0.5 * step
        trajectory.write(molecule)

# The same atoms throughout, with the header written again: for a constraint, then for pbc
rewritten = Atoms("H2O", positions=[(0, 0, 0), (0.9, 0, 0), (0, 0.8, 0)], cell=[10, 10, 10], pbc=True)
with Trajectory(here / "ase_header_rewrite.traj", "w") as trajectory:
    trajectory.write(rewritten)
    rewritten.set_constraint(FixAtoms([0]))
    trajectory.write(rewritten)
    rewritten.pbc = False
    trajectory.write(rewritten)
