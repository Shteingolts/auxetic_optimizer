# NOTE: This script can be modified for different atomic structures, 
# units, etc. See in.elastic for more info.
#

# Define the finite deformation size. Try several values of this
# variable to verify that results do not depend on it.
# variable up equal 1.0e-1
variable up equal 1.0e-4
 
# Define the amount of random jiggle for atoms
# This prevents atoms from staying on saddle points
variable atomjiggle equal 1.0e-5

# Uncomment one of these blocks, depending on what units
# you are using in LAMMPS and for output

# metal units, elastic constants in eV/A^3
#units		metal
#variable cfac equal 6.2414e-7
#variable cunits string eV/A^3

# metal units, elastic constants in GPa
# units		metal
# variable cfac equal 1.0e-4
# variable cunits string GPa

# real units, elastic constants in GPa
#units		real
#variable cfac equal 1.01325e-4
#variable cunits string GPa

# COPIED FROM MODIFIED SCRIPT
#---------------------------------
# Choose potential
# metal units, elastic constants in GPa
units		metal
atom_style	full
bond_style      harmonic
angle_style     harmonic
dihedral_style  zero
special_bonds   amber
variable cfac equal 1.0e-4
variable cunits string GPa
boundary	p p p
#---------------------------------
read_data network.lmp


# Define minimization parameters
variable etol equal 0.0 
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2

# Need to set mass to something, just to satisfy LAMMPS
mass 1 1.0e-20

