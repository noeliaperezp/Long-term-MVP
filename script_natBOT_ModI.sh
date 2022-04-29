#!/bin/bash
#$ -cwd

rm script_natBOT_f_ModI.sh.*

#script_natBOT.sh <d> <LAMB> <NCRO> <BETA> <AVEs> <AVEh>
#To be used only in a machine with /state/partition1 directory

#Check number of arguments
if [ $# -ne 2 ]  
then
	echo "Usage: $0 <d> <REPS>" 
	exit 1
fi

#Set arguments
d=$1
REPS=$2


#Variables
LAMB=0.34
NIND=10000
NCRO=1000
BETA=0.2
AVEs=0.02
AVEh=0.2
LAMBL=0.015
Vs=0


#Working directory
WDIR=$PWD 
mkdir -p $WDIR/NATBOT/ModI_L$LAMB.k$NCRO.s$AVEs.h$AVEh.N$NIND.Vs$Vs
DIR="NATBOT/ModI_L$LAMB.k$NCRO.s$AVEs.h$AVEh.N$NIND.Vs$Vs"

#Scratch directory
mkdir -p /state/partition1/noeliaNATf$d/$SLURM_JOBID/

#Copy all files in scratch directory
cp seedfile /state/partition1/noeliaNATf$d/$SLURM_JOBID/
cp naturalfBOT /state/partition1/noeliaNATf$d/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd /state/partition1/noeliaNATf$d/$SLURM_JOBID

###########################################################################

for ((i=$REPS; i<= $REPS ; i++)); do
START=$(date +%s)
SEED=$(( 1234 + $i ))

time ./naturalfBOT>>out<<@
0
$SEED
$NIND	N
0	PS(99=random)
20	Lenght genome (99=free)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_a
$LAMBL	Lambda_L
0.0	absolute effect of lethal (QT): normal (aL,aL)
0	random proportion(0) or large mutants (1)
1.0	Psi
$BETA	beta_s
$BETA	beta_a
$AVEs	ave |s|
$AVEs	ave |a|
0.0	PP_s
0.5	PP_a
2	dom model (0=cnt; 1:Deng, 2:CK94 gamma)
$AVEh	h_s (mod 0), k_s (mod 1)
$AVEh	ave h_s (mod 2)
$AVEh	h_a (mod 0), k_a (mod 1)
$AVEh	ave h_a (mod 2)
99	rho (99:a=s)
$Vs	Vs
1	multi(1), add(2)
10000	generations
2000	gen/block
0	GENBOT
@

cat popfile >> rep${i}.POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh
cat datafile >> rep${i}.DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh
cat genfile >> rep${i}.GENFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "naturalf took 		$DIFF seconds" >> rep${i}.timefile

#Copy output files to main directory
cp -r /state/partition1/noeliaNATf$d/$SLURM_JOBID/rep${i}.POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $WDIR/$DIR/
cp -r /state/partition1/noeliaNATf$d/$SLURM_JOBID/rep${i}.DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $WDIR/$DIR/
cp -r /state/partition1/noeliaNATf$d/$SLURM_JOBID/rep${i}.GENFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $WDIR/$DIR/
cp -r /state/partition1/noeliaNATf$d/$SLURM_JOBID/rep${i}.timefile $WDIR/$DIR/

done

###########################################################################

#Cleaning of scratch
rm -r /state/partition1/noeliaNATf$d/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
