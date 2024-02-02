#!/bin/bash

cd /home/sc_enigma/Projects/MD_UIO66/analysis/rotacf/

for JOB_NAME in uio66_tempo uio66_tempo_water128 uio66_tempo_water256 uio66_tempo_water1038
do
	echo 4 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME.xvg
done

for JOB_NAME in uio66_tempo_removed_linkers uio66_tempo_replaced_linkers uio66_tempo_water_removed
do
	echo 5 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed1/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed1/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME/removed1.xvg

	echo 6 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed2/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed2/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME/removed2.xvg

	echo 7 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed3/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed3/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME/removed3.xvg
done

for JOB_NAME in uio66_tempo_removed_linker_and_cluster
do
	echo 8 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed1/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed1/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME/removed1.xvg

	echo 9 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed2/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed2/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME/removed2.xvg

	echo 10 | gmx_mpi rotacf \
	-f /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed3/prod/traj_comp.xtc \
	-s /media/sc_enigma/_data/Projects/MD_UIO66/production/$JOB_NAME/removed3/prod/prod.tpr \
	-n index.ndx -o N_O/$JOB_NAME/removed3.xvg
done