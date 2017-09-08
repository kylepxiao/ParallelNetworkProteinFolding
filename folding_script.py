from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
from MMTK.Dynamics import VelocityVerletIntegrator, Heater, TranslationRemover, RotationRemover
from MMTK.Visualization import view
from MMTK.Trajectory import Trajectory, TrajectoryOutput, RestartTrajectoryOutput, StandardLogOutput, trajectoryInfo

from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain, Protein

from MMTK.Trajectory import Trajectory

from MMTK.DCD import writeDCDPDB

import time

file = open("output.txt", 'w')

start = time.time()

#
# First problem: construct an all-atom model from a structure without
# hydrogens. This is the standard problem when using an all-atom force

# field with crystallographic structures.
#
#
# Load the PDB file.
configuration = PDBConfiguration('insulin.pdb')


# Construct the peptide chain objects. This also constructs positions
# for any missing hydrogens, using geometrical criteria.
chains = configuration.createPeptideChains()

# Make the protein object.
#insulin = Protein(chains)

# Define system
universe = InfiniteUniverse(Amber99ForceField(mod_files=['frcmod.ff99SB']))
universe.protein = Protein(chains)

# Initialize velocities
universe.initializeVelocitiesToTemperature(50.*Units.K)
print 'Temperature: ', universe.temperature()
print 'Momentum: ', universe.momentum()
print 'Angular momentum: ', universe.angularMomentum()
file.write('Temperature: ' + str(universe.temperature()) + "\n")
file.write('Momentum: ' + str(universe.momentum()) + "\n")
file.write('Angular momentum: ' + str(universe.angularMomentum()) + "\n")

# Create integrator
integrator = VelocityVerletIntegrator(universe, delta_t=1.*Units.fs)

# Heating and equilibration
integrator(steps=1000,
                    # Heat from 50 K to 300 K applying a temperature
                    # change of 0.5 K/fs; scale velocities at every step.
	   actions=[Heater(50.*Units.K, 300.*Units.K, 0.5*Units.K/Units.fs,
                           0, None, 1),
                    # Remove global translation every 50 steps.
		    TranslationRemover(0, None, 50),
                    # Remove global rotation every 50 steps.
		    RotationRemover(0, None, 50),
                    # Log output to screen every 100 steps.
                    StandardLogOutput(100)])

					
print("Time: " + str(time.time() - start)),
file.write("Time: " + str(time.time() - start))
					
# "Production" run
trajectory = Trajectory(universe, "insulin.nc", "w", "A simple test case")

universe.protein.writeToFile('folded_protein.pdb')

'''integrator(steps=1000,
                      # Remove global translation every 50 steps.
           actions = [TranslationRemover(0, None, 50),
                      # Remove global rotation every 50 steps.
                      RotationRemover(0, None, 50),
                      # Write every second step to the trajectory file.
                      TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "configuration"),
                                       0, None, 2),
                      # Write restart data every fifth step.
                      RestartTrajectoryOutput("restart.nc", 5),
                      # Log output to screen every 10 steps.
                      StandardLogOutput(10)])'''

writeDCDPDB(trajectory.configuration, 'insulin_new.dcd', 'insulin_new.pdb')					  

trajectory.close()

# Print information about the trajectory file
print "Information about the trajectory file 'insulin.nc':"
print trajectoryInfo('insulin.nc')

file.write("Information about the trajectory file 'insulin.nc':")
file.write(trajectoryInfo('insulin.nc'))

# Reopen trajectory file
trajectory = Trajectory(universe, "insulin.nc", "r")

writeDCDPDB(trajectory.configuration, 'insulin_new.dcd', 'insulin_new.pdb')

# Read step 10 and display configuration
#step10 = trajectory[10]
#view(universe, step10['configuration'])

# Print the kinetic energy along the trajectory
'''print "Kinetic energy along trajectory:"
file.write("Kinetic energy along trajectory:")
file.write(trajectory.kinetic_energy)'''

trajectory.close()
file.close()