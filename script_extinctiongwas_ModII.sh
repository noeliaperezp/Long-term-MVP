#!/bin/bash
#$ -cwd

#...............................................................................#

#  System      | Mating             | Contributions    | Couples  | Limitation of K  | Selfing  |
# -------------|--------------------|------------------|----------|------------------|----------|
#  (0) RC      | Random             | Random           | Monogamy | Yes              | No       |
#  (1) EC      | Random or avoid FS | Equalization (2) | Monogamy | Yes              | No       |
#  (2) CM      | Circular           | Equalization (2) | Monogamy | Yes              | No       |
#  (3) PFS     | Random (% FS)      | Random           | Monogamy | Yes              | No       |
#  (4) RCpol   | Random             | Random           | Polygamy | No               | Optional |
#  (5) RCW     | Random             | Random           | Polygamy | Yes              | Optional |

# FS = full-sib mating
#...............................................................................#

rm script_extinctiongwas_ModII.sh.*

################################# ARGUMENTS ###############################

#script_extinctiongwas.sh <d> <NP> <Vs> <NMAX> <MFEC> <NEU>
#To be used only in a machine with /state/partition1 directory

#Check number of arguments
if [ $# -ne 6 ]  
then
	echo "Usage: $0 <d> <NP> <Vs> <NMAX> <MFEC> <REPS_BP>" 
	exit 1
fi

#Set arguments
d=$1
NP=$2
Vs=$3
NMAX=$4		
MFEC=$5
REPS_BP=$6

################################ VARIABLES ################################

TYPE=5		#(R=0; EC=1; CM=2; PFS=3; RCpol=4; RCW=5)
FV=1		#fec_via (0: equal, 1: 1/3-2/3, 2: viability, 3: fecundity)

#Parameters
LAMB=0.22
NCRO=1000
AVEs=0.05
AVEh=0.2
BETA=0.45
LAMBL=0.015

############################## DIRECTORIES ################################

#Working directory
WDIR=$PWD 

##################
#Natural population
NATDIR="NATBOT/ModII_L$LAMB.k$NCRO.s$AVEs.h$AVEh.N$NP.Vs$Vs"

##################
#Output directory
if [[ $TYPE == 0 ]]; then case="RC" 
elif [[ $TYPE == 1 ]]; then case="EC"
elif [[ $TYPE == 2 ]]; then case="CM"
elif [[ $TYPE == 3 ]]; then case="PFS.$PFS"
elif [[ $TYPE == 4 ]]; then case="RCpol"
else case="RCW"
fi

if [[ $FV == 0 ]]; then fit="fec&via"
elif [[ $FV == 1 ]]; then fit="fec&via_13_23"
elif [[ $FV == 2 ]]; then fit="via"
else fit="fec"
fi

mkdir -p $WDIR/extinctiongwas_minNe_Results/N$NP.Vs$Vs.L$LAMB.K$NCRO.s$AVEs.h$AVEh.N$NMAX/$case.$fit.mfec$MFEC
DIR="extinctiongwas_minNe_Results/N$NP.Vs$Vs.L$LAMB.K$NCRO.s$AVEs.h$AVEh.N$NMAX/$case.$fit.mfec$MFEC"

##################
#Scratch directory
mkdir -p /state/partition1/noeliaEXT$d/$SLURM_JOBID/
SCDIR="/state/partition1/noeliaEXT$d" 

######################## TRANSFER OF FILES TO SCRATCH #####################
 
#Copy all files in scratch directory
cp seedfile $SCDIR/$SLURM_JOBID/
cp extinctiongwas_minNe $SCDIR/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd $SCDIR/$SLURM_JOBID

###########################################################################

for ((i=$REPS_BP; i<= $REPS_BP ; i++)); do
START=$(date +%s)

cp $WDIR/$NATDIR/rep${i}.POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/popfile
cp $WDIR/$NATDIR/rep${i}.DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/datafile

############# EXTINCTIONGWAS #############
 
time ./extinctiongwas_minNe>>out<<@
0
-99
$NMAX	NMAX
$MFEC	MAXFEC
$TYPE	type (R=0; EC=1; CM=2; PFS=3; RCpol=4; RCW=5)
1	Avoidance of full-sib matings (EC)(0 no, 1 yes)
0	Selfing in case of polygamy (0 no, 1 yes)
0	PFS (For PFS, proportion of partial full-sib mating; 0=random)
20	L (morgans; 99=FreeRec)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_s
$LAMBL	Lambda_L
$BETA	beta_s
$AVEs	ave |s|
1	dom (0 constant, 1 CK94)
$AVEh	ave h_s
0.0	ave_aL (rho != 99)
$BETA	beta_a (rho != 99)
$AVEs	ave_a (rho != 99)
$AVEh	ave_ha (rho != 99)
99.0	rho
0.5	Positives for the trait
$Vs	Stabilizing selection (Vs)
0.0	optimal
1.0	VE
0	relaxation factor (0: no, 1:yes)
0	neutral (0: no, 1:yes)
1	scaling (0: no, 1:yes)
$FV	fec_via (0: equal, 1: 1/3-2/3, 2: viability, 3: fecundity)
0	Fecundity as a maternal component (0: no; 1: yes)
0	Pedfile information (0: survivors; 1: zygotes)
1	Environmental source of mortality (A=0.3) (0: no; 1: yes)
0	Natural catastrophes (P=0.058) (0: no; 1: yes)
1000	generations
99	Change type to: RC(0), EC(1), CM(2), PFS(3), without change(99)
99	Generation of change
100	replicates
@

cp -r $SCDIR/$SLURM_JOBID/genfile.dat $WDIR/$DIR/GENFILE${i}
rm genfile.dat
rm datafile
rm popfile

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "extinctiongwas BP ${i} took $DIFF seconds\n" >> timefile${i}
cp -r $SCDIR/$SLURM_JOBID/timefile${i} $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/out $WDIR/$DIR/out${i}

done

################# TRANSFER OF OTHER FILES TO DIRECTORY ####################

cp -r $SCDIR/$SLURM_JOBID/seedfile $WDIR
cp -r $SCDIR/$SLURM_JOBID/dfilename*.dat $WDIR/$DIR/

###########################################################################
######################### CLEANING OF SCRATCH #############################

rm -r $SCDIR/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
