# GROMACS (Fluctuating Hydrodynamics / Molecular Dynamics two-way coupling)

Modified by Ivan Korotkin.

# Energy Minimisation

	gmx grompp -f minim.mdp -c conf.gro -p topol.top -o em.tpr -maxwarn 1
	gmx mdrun -v -deffnm em

	gmx grompp -f grompp.mdp -c em.gro -o afterminim.tpr
	gmx mdrun -nt 1 -s afterminim.tpr -o afterminim_traj.trr -c afterminim_confout.gro -v


