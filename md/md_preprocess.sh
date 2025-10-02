#!/usr/bin/bash

# Need to prepare: protein.pdb (added hydrogen, saved by Swiss-PDBViewer); ligand.pdbqt; *.mdp
# Usage: bash md_preprocess.sh <protein_pdb> <docked_pdbqt> <gpu_id>
# Example: bash md_preprocess.sh 2rku.pdb top1M_split1_60280_ZINC001783543529_1_out.pdbqt 0 &

OBABEL=/path/bin/obabel
SOBTOP_PATH=/path/sobtop_1.0/
CODE_PATH=/path/sobtop_1.0/
WD=$(pwd)

printf '6\n' | /path/install/bin/gmx pdb2gmx -f $1 -o protein.gro -water tip3p -ignh

$OBABEL -i pdbqt $2 -o mol2 -O ligand.mol2 -h -l 1

cd $SOBTOP_PATH
printf "2\n$WD/ligand.gro\n1\n2\n4\n\n$WD/ligand.itp\n0\n" | ./sobtop $WD/ligand.mol2

cd $WD
sed -i 's/1  UNL1/MOL/' ligand.mol2

python $CODE_PATH/md_preprocess.py 1 $WD protein.gro

/path/install/bin/gmx editconf -f complex.gro -o newbox.gro -d 1.0 -c -bt dodecahedron

/path/install/bin/gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

/path/install/bin/gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2

printf '15\n' | /path/install/bin/gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

/path/install/bin/gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2

/path/install/bin/gmx mdrun -ntomp 8 -ntmpi 1 -gpu_id $3 -v -deffnm em

printf '0 & ! a H*\nq\n' | /path/install/bin/gmx make_ndx -f ligand.gro -o index_ligand.ndx

printf '3\n' | /path/install/bin/gmx genrestr -f ligand.gro -n index_ligand.ndx -o posre_ligand.itp -fc 1000 1000 1000

python $CODE_PATH/md_preprocess.py 2 $WD

printf '1|13\nq\n' | /path/install/bin/gmx make_ndx -f em.gro -o index.ndx

/path/install/bin/gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
/path/install/bin/gmx mdrun -ntomp 8 -deffnm nvt -gpu_id $3 -v

/path/install/bin/gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
/path/install/bin/gmx mdrun -ntomp 8 -deffnm npt -gpu_id $3 -v

# Change `nsteps` in `md.mdp`. 50000000 for 100ns; 25000000 for 50ns; 12500000 for 25ns  

/path/install/bin/gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_100.tpr -maxwarn 2
/path/install/bin/gmx mdrun -ntomp 8 -gpu_id $3 -s md_0_100.tpr -v


printf '1 0\n' | /path/install/bin/gmx trjconv -s md_0_100.tpr -f traj_comp.xtc -o md_0_100_center.xtc -center -pbc mol -ur compact

printf '4 0\n' | /path/install/bin/gmx trjconv -s md_0_100.tpr -f md_0_100_center.xtc -o md_0_100_fit.xtc -fit rot+trans

printf '22 & ! a H*\nname 23 complex_Heavy\n13 & ! a H*\nname 24 ligand_Heavy\nq\n' | /path/install/bin/gmx make_ndx -f em.gro -n index.ndx

printf '23 24\n' | /path/install/bin/gmx rms -s em.tpr -f md_0_100_center.xtc -n index.ndx -tu ns -o rmsd_lig.xvg

# printf '23 23\n' | /path/install/bin/gmx rms -s em.tpr -f md_0_100_center.xtc -n index.ndx -tu ns -o rmsd_complex.xvg

printf '1\n' | /path/install/bin/gmx rmsf -s md_0_100.tpr -f md_0_100_center.xtc -o rmsf_protein.xvg -res

/path/install/bin/gmx trjconv -f md_0_100_fit.xtc -o trj_80-100ns -b 80000 -e 100000 -dt 10

printf '1 13\n' | /path/install/bin/gmx hbond -s md_0_100.tpr -f trj_80-100ns.xtc -num hbond.xvg

# conda activate python39
# $GMX_MMPBSA -O -i mmpbsa.in -cs md_0_100.tpr -ci index.ndx -cg 1 13 -ct trj_80-100ns.xtc -cp topol.top
