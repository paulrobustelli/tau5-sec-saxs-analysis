#!/usr/bin/env python
"""Harmonic Gaussian-distance variant; original Jeffreys runner preserved."""

import argparse
import re
import sys
import numpy as np
from openmm import (CMAPTorsionForce, CustomBondForce, LangevinMiddleIntegrator,
                    Platform)
from openmm.app import (PDBFile, ForceField, Simulation, CutoffNonPeriodic,
                        HBonds, XTCReporter, StateDataReporter,
                        CheckpointReporter)
from openmm.unit import kelvin, picosecond, picoseconds


def read_cmap_file(path):
    out, title, values = {}, None, []
    for raw in open(path):
        line = raw.rstrip("\n")
        if re.search(r"-\(", line) and re.search(r", type \d+", line):
            if title is not None:
                out[title] = np.reshape(values, (24, 24))
            title, values = line.split()[1].rstrip(","), []
        elif title is not None and line and not line.startswith("# "):
            values += [float(x) for x in line.split()]
    if title is not None:
        out[title] = np.reshape(values, (24, 24))
    return out


def shift_matrix_values(matrix, step_0=1, step_1=0):
    return np.roll(np.roll(matrix, step_0, axis=0), step_1, axis=1)


def choose_platform(name):
    available = [Platform.getPlatform(i).getName()
                 for i in range(Platform.getNumPlatforms())]
    if name == "auto":
        name = next((x for x in ("CUDA", "OpenCL", "CPU", "Reference")
                     if x in available), None)
    if name not in available:
        raise ValueError(f"Platform {name!r} unavailable; found {available}")
    platform = Platform.getPlatformByName(name)
    properties = {}
    if name == "CUDA":
        properties["Precision"] = "double"
    elif name == "OpenCL":
        properties["Precision"] = "single"
    return platform, properties, available


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("-ff", required=True)
    p.add_argument("-pdb", required=True)
    p.add_argument("-xtc", required=True)
    p.add_argument("-cmap")
    p.add_argument("-baies")
    p.add_argument("-temp", type=float, default=298.15)
    p.add_argument("-friction", type=float, default=0.01)
    p.add_argument("-tstep", type=float, default=0.002)
    p.add_argument("-Nsteps", type=int, default=20000)
    p.add_argument("-save_xtc", type=int, default=500)
    p.add_argument("-save_stdout", type=int, default=500)
    p.add_argument("-checkpoint")
    p.add_argument("-restart", help="OpenMM binary checkpoint to resume")
    p.add_argument("-restart_state", help="portable OpenMM State XML to resume")
    p.add_argument("-state_out", help="write a portable OpenMM State XML at exit")
    p.add_argument("-platform", default="auto")
    p.add_argument("-seed", type=int, default=20261008)
    return p


def main():
    args = parser().parse_args()
    pdb = PDBFile(args.pdb)
    ff = ForceField(args.ff)
    system = ff.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                             constraints=HBonds)

    if args.cmap:
        cmap_force = CMAPTorsionForce()
        cmap_list = {}
        for key, matrix in read_cmap_file(args.cmap).items():
            formatted = shift_matrix_values(matrix*4.184, 12, 12).transpose()
            cmap_list[key] = cmap_force.addMap(24, formatted.flatten())
        residues = {}
        for index, atom in enumerate(pdb.topology.atoms()):
            key = (int(atom.residue.id), atom.residue.chain.id)
            residues.setdefault(key, {"resname": atom.residue.name})[atom.name] = index
        for resid, chain in residues:
            if ((resid-1, chain) in residues and (resid+1, chain) in residues
                    and "C" in residues[(resid-1, chain)]
                    and "N" in residues[(resid+1, chain)]):
                cur, prev, nxt = residues[(resid, chain)], residues[(resid-1, chain)], residues[(resid+1, chain)]
                a1, a2, a3, a4 = prev["C"], cur["N"], cur["CA"], cur["C"]
                suffix = "PRO" if nxt["resname"] == "PRO" else "XXX"
                cmap_force.addTorsion(cmap_list[f'{cur["resname"]}-({suffix})'],
                                      a1, a2, a3, a4, a2, a3, a4, nxt["N"])
        system.addForce(cmap_force)
        print(f"CMAP torsions: {cmap_force.getNumTorsions()}")

    if args.baies:
        kt = args.temp*8.314462618e-3
        force = CustomBondForce("0.5*kT*((r-mu)/sigma)^2")
        force.addGlobalParameter("kT", kt)
        force.addPerBondParameter("mu")
        force.addPerBondParameter("sigma")
        data = np.atleast_2d(np.loadtxt(args.baies))
        if not np.isfinite(data).all() or np.any(data[:,4] <= 0):
            raise ValueError("Restraint parameters must be finite with positive sigma")
        for row in data:
            force.addBond(int(row[1]-1), int(row[2]-1), [row[3], row[4]])
        system.addForce(force)
        print(f"Harmonic restraints: {force.getNumBonds()}; U=0.5*kT*((r-mu)/sigma)^2; sigma is SD")

    integ = LangevinMiddleIntegrator(args.temp*kelvin,
                                     args.friction/picosecond,
                                     args.tstep*picoseconds)
    integ.setRandomNumberSeed(args.seed)
    platform, properties, available = choose_platform(args.platform)
    print(f"OpenMM platforms: {available}; using {platform.getName()} {properties}")
    sim = Simulation(pdb.topology, system, integ, platform, properties)
    sim.context.setPositions(pdb.positions)
    if args.restart and args.restart_state:
        raise ValueError("Use only one of -restart or -restart_state")
    if args.restart_state:
        sim.loadState(args.restart_state)
        print(f"Restarted from {args.restart_state} at step {sim.currentStep}")
    elif args.restart:
        sim.loadCheckpoint(args.restart)
        print(f"Restarted from {args.restart} at step {sim.currentStep}")
    else:
        print("Minimizing...")
        sim.minimizeEnergy(maxIterations=10000)
        sim.context.setVelocitiesToTemperature(args.temp*kelvin, args.seed)
    sim.reporters.append(XTCReporter(args.xtc, args.save_xtc))
    sim.reporters.append(StateDataReporter(sys.stdout, args.save_stdout,
                           step=True, potentialEnergy=True, temperature=True,
                           speed=True, remainingTime=True,
                           totalSteps=sim.currentStep + args.Nsteps, separator="\t"))
    if args.checkpoint:
        sim.reporters.append(CheckpointReporter(args.checkpoint,
                                                 args.save_stdout))
    sim.step(args.Nsteps)
    if args.state_out:
        sim.saveState(args.state_out)
        print(f"Saved portable state to {args.state_out}")
    print("Completed", args.Nsteps, "steps")


if __name__ == "__main__":
    main()
