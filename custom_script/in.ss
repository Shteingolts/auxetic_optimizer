# FENE beadspring 

############### deformation ################
variable T equal 1.00
variable dt equal 0.01
## can be axial or shear, in any direction
## set strain, strain rate, dimension
## set replica number
variable strain equal 0.2
variable srate equal 1e-4
variable dim equal 3
variable rep equal 0
variable M equal 1000 # averaging time for properties
############################################

units		lj
atom_style	bond
special_bonds lj/coul 0.0 1.0 1.0
#for optional breakable bonds
#special_bonds lj/coul 1 1 1 # quartic

read_data	data.eq.f3.n16

neighbor	1.0 bin
neigh_modify	every 1 delay 0

bond_style      fene 
bond_coeff	* 30.0 1.5 1.0 1.0

#for optional breakable bonds
#bond_style quartic
#bond_coeff * 1200 -0.55 0.25 1.3 34.6878

pair_style	lj/cut 2.5
pair_coeff	* * 1.0 1.0 2.5
pair_modify	shift yes

fix		1 all nve
fix		2 all langevin $T $T 100 904297

timestep ${dt}

############### computes and analysis  ###################

compute 1 all displace/atom
thermo         $M 

############### deformation ################

variable pxx equal pxx
variable pyy equal pyy
variable pzz equal pzz
variable sxx equal -pxx
variable syy equal -pyy
variable szz equal -pzz
variable pxy equal pxy
variable pxz equal pxz
variable pyz equal pyz
variable von equal sqrt(0.5*((v_pxx-v_pyy)^2+(v_pyy-v_pzz)^2+(v_pzz-v_pxx)^2)+3*(v_pxy^2+v_pyz^2+v_pxz^2))
variable hyd equal (v_pxx+v_pyy+v_pzz)/3
variable lx  equal lx
variable ly  equal ly
variable lz  equal lz
variable vol equal vol
variable lx0 equal ${lx}
variable ly0 equal ${ly}
variable lz0 equal ${lz}
variable xy equal xy
variable yz equal yz
variable xz equal xz
variable vol equal vol
variable xStrain equal (v_lx-v_lx0)/v_lx0
variable yStrain equal (v_ly-v_ly0)/v_ly0
variable zStrain equal (v_lz-v_lz0)/v_lz0
variable xyStrain equal v_xy/v_ly0
variable yzStrain equal v_yz/v_lz0
variable zxStrain equal v_xz/v_lx0
variable xelong equal ((v_pyy+v_pzz)/2.0-v_pxx) #elongation stress in GPa
variable yelong equal ((v_pxx+v_pzz)/2.0-v_pyy)
variable zelong equal ((v_pxx+v_pyy)/2.0-v_pzz)
variable strainsteps equal ${strain}/${dt}/${srate}
print "strainsteps ${strainsteps}"

if "${dim} > 2" then "change_box all triclinic"

if "${dim} == 0" then "fix     3 all deform 1 x erate ${srate} units box remap x"
if "${dim} == 1" then "fix     3 all deform 1 y erate ${srate} units box remap x"
if "${dim} == 2" then "fix     3 all deform 1 z erate ${srate} units box remap x"
if "${dim} == 3" then "fix     3 all deform 1 xy erate ${srate} units box remap x"
if "${dim} == 4" then "fix     3 all deform 1 yz erate ${srate} units box remap x"
if "${dim} == 5" then "fix     3 all deform 1 xz erate ${srate} units box remap x"

fix    sysstr all ave/time 1 $M $M v_pxx v_pyy v_pzz v_pxy v_pxz v_pyz v_xelong v_yelong v_zelong v_von v_hyd v_lx v_ly v_lz v_vol v_xStrain v_yStrain v_zStrain v_xyStrain v_yzStrain v_zxStrain file press.${rep}.d${dim}.${srate}.$T

variable nsteps equal round(${strainsteps})
print "nsteps ${nsteps}"

thermo_style custom step pe ke temp press v_vol xy yz xz v_xStrain v_yStrain v_zStrain v_xyStrain v_yzStrain v_zxStrain
thermo_modify flush yes

dump 1 all atom $M dump.$T.${dim}.${rep}.lammpstrj
dump 2 all custom $M dump.$T.${dim}.${rep}.displacements.lammpstrj x y z xu yu zu c_1[1] c_1[2] c_1[3] c_1[4]

run  ${nsteps}
