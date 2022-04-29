
/* naturalfBOT.c (modelo de fecundidad) BOTTLENECK N=10,100 generations t=10-20-100*/

/* *********************************************** */

#include "libhdr"
#include "ranlib.h"
#define NN 60001  /* max number of (gen/blocks) */
#define MM 10001  /* max number of NIND */
#define MMM 2001  /* max NCRO */
#define maxmpt 5001
#define normalthreshold 30

int muts, mutsL, marker, whatproportion, GENBOT;
int lastinmutantspoissontable, lastinmutantspoissontableL;
int NIND, NCRO, NLOCI, TOTLOCI, block, model_selection;
int numberofrecom[NN], gen, generations, i, j, k, l, m;
int totalmutant[NN], doublemutant[NN], recurrentmutant[NN];
int mutants_ocurred, RM[NN], initialgen[MMM][31], dom;
int ran_i, ran_k, ran_l, ran_h, countNoSS, NoSS_k[6001], NoSS_l[6001];

int gm[MM][MMM][2], sm[MM][MMM][2];
double genvalm_s[MM], genvalm_a[MM], pm_a[MM], pm_s[MM];

double Lambda_a, Lambda_L, Psi, threshold, h2a;
double PS, L, AA, Aa, aa, q[MMM][31];
double meanmuts_a[NN] , meanmuts_s[NN] , meanmutsL_s[NN];
double ave_s, addedfixed_s, alpha_s, beta_s, k_s, ave_hs, epsilon_s, P_s, zs, h1_s, h2_s, h3_s, h4_s, h_sa;
double ave_a, ave_aL, addedfixed_a, alpha_a, beta_a, k_a, ave_ha, epsilon_a, P_a, za, h1_a, h2_a, h3_a, h4_a;
double rho, Vs, bSS, VE;
double mean_s[NN], mean_s2[NN], mean_s3[NN], mean_s4[NN], mean_hs[NN], mean_sa[NN];
double mean_real_a[NN], mean_a[NN], mean_a2[NN], mean_a3[NN], mean_a4[NN], mean_ha[NN];
double mean_a2_s[NN], h1_a_s, h2_a_s, mean_a_s[NN];
double mutantspoissontable[maxmpt], mutantspoissontableL[maxmpt];
double freqn[NN], CHn[NN], fixn[NN], lostn[NN], othern[NN];
double freq[NN], CH[NN], fix[NN], lost[NN], other[NN];
double s[MMM][31], hs[MMM][31], a[MMM][31], ha[MMM][31];
double cumH[MMM][31];
double timetolossn[NN], timetofixn[NN], timetoloss[NN], timetofix[NN];
double Ave_CHn, Ave_NSSn, Ave_timetolossn, Ave_timetofixn;
double Ave_freq, Ave_CH, Ave_NSS, Ave_timetoloss, Ave_timetofix;
double countgenes, d_s, d_a, alfa_a, alfa_s;

double bvgam();

/* ********************** Distribution on q ************************ */

struct acc ng_all, a_all, s_all, ha_all, hs_all, VA_all, VD_all;
struct acc Ave_ng_all, Ave_a_all, Ave_s_all, Ave_ha_all, Ave_hs_all, Ave_VA_all, Ave_VD_all;

struct acc q_ng_00_01, q_ng_01_02, q_ng_02_03;
struct acc q_ng_03_04, q_ng_04_05, q_ng_05_06, q_ng_06_07;
struct acc q_ng_07_08, q_ng_08_09, q_ng_09_10;
struct acc q_s_00_01, q_s_01_02, q_s_02_03, q_s_03_04;
struct acc q_s_04_05, q_s_05_06, q_s_06_07, q_s_07_08;
struct acc q_s_08_09, q_s_09_10;
struct acc q_a_00_01, q_a_01_02, q_a_02_03, q_a_03_04;
struct acc q_a_04_05, q_a_05_06, q_a_06_07, q_a_07_08;
struct acc q_a_08_09, q_a_09_10;
struct acc q_hs_00_01, q_hs_01_02, q_hs_02_03, q_hs_03_04;
struct acc q_hs_04_05, q_hs_05_06, q_hs_06_07, q_hs_07_08;
struct acc q_hs_08_09, q_hs_09_10;
struct acc q_ha_00_01, q_ha_01_02, q_ha_02_03, q_ha_03_04;
struct acc q_ha_04_05, q_ha_05_06, q_ha_06_07, q_ha_07_08;
struct acc q_ha_08_09, q_ha_09_10;
struct acc q_VA_00_01, q_VA_01_02, q_VA_02_03, q_VA_03_04;
struct acc q_VA_04_05, q_VA_05_06, q_VA_06_07, q_VA_07_08;
struct acc q_VA_08_09, q_VA_09_10;
struct acc q_VD_00_01, q_VD_01_02, q_VD_02_03, q_VD_03_04;
struct acc q_VD_04_05, q_VD_05_06, q_VD_06_07, q_VD_07_08;
struct acc q_VD_08_09, q_VD_09_10;

struct acc Ave_q_ng_00_01, Ave_q_ng_01_02, Ave_q_ng_02_03;
struct acc Ave_q_ng_03_04, Ave_q_ng_04_05, Ave_q_ng_05_06, Ave_q_ng_06_07;
struct acc Ave_q_ng_07_08, Ave_q_ng_08_09, Ave_q_ng_09_10;
struct acc Ave_q_s_00_01, Ave_q_s_01_02, Ave_q_s_02_03, Ave_q_s_03_04;
struct acc Ave_q_s_04_05, Ave_q_s_05_06, Ave_q_s_06_07, Ave_q_s_07_08;
struct acc Ave_q_s_08_09, Ave_q_s_09_10;
struct acc Ave_q_a_00_01, Ave_q_a_01_02, Ave_q_a_02_03, Ave_q_a_03_04;
struct acc Ave_q_a_04_05, Ave_q_a_05_06, Ave_q_a_06_07, Ave_q_a_07_08;
struct acc Ave_q_a_08_09, Ave_q_a_09_10;
struct acc Ave_q_hs_00_01, Ave_q_hs_01_02, Ave_q_hs_02_03, Ave_q_hs_03_04;
struct acc Ave_q_hs_04_05, Ave_q_hs_05_06, Ave_q_hs_06_07, Ave_q_hs_07_08;
struct acc Ave_q_hs_08_09, Ave_q_hs_09_10;
struct acc Ave_q_ha_00_01, Ave_q_ha_01_02, Ave_q_ha_02_03, Ave_q_ha_03_04;
struct acc Ave_q_ha_04_05, Ave_q_ha_05_06, Ave_q_ha_06_07, Ave_q_ha_07_08;
struct acc Ave_q_ha_08_09, Ave_q_ha_09_10;
struct acc Ave_q_VA_00_01, Ave_q_VA_01_02, Ave_q_VA_02_03, Ave_q_VA_03_04;
struct acc Ave_q_VA_04_05, Ave_q_VA_05_06, Ave_q_VA_06_07, Ave_q_VA_07_08;
struct acc Ave_q_VA_08_09, Ave_q_VA_09_10;
struct acc Ave_q_VD_00_01, Ave_q_VD_01_02, Ave_q_VD_02_03, Ave_q_VD_03_04;
struct acc Ave_q_VD_04_05, Ave_q_VD_05_06, Ave_q_VD_06_07, Ave_q_VD_07_08;
struct acc Ave_q_VD_08_09, Ave_q_VD_09_10;

/* ********************** Distribution on a *********************** */

struct acc a_ng_00_02, a_ng_02_04, a_ng_04_06;
struct acc a_ng_06_08, a_ng_08_10, a_ng_10_12, a_ng_12_14;
struct acc a_ng_14_16, a_ng_16_18, a_ng_18_20, a_ng_20;
struct acc a_s_00_02, a_s_02_04, a_s_04_06;
struct acc a_s_06_08, a_s_08_10, a_s_10_12, a_s_12_14;
struct acc a_s_14_16, a_s_16_18, a_s_18_20, a_s_20;
struct acc a_q_00_02, a_q_02_04, a_q_04_06;
struct acc a_q_06_08, a_q_08_10, a_q_10_12, a_q_12_14;
struct acc a_q_14_16, a_q_16_18, a_q_18_20, a_q_20;
struct acc a_ha_00_02, a_ha_02_04, a_ha_04_06;
struct acc a_ha_06_08, a_ha_08_10, a_ha_10_12, a_ha_12_14;
struct acc a_ha_14_16, a_ha_16_18, a_ha_18_20, a_ha_20;
struct acc a_VA_00_02, a_VA_02_04, a_VA_04_06;
struct acc a_VA_06_08, a_VA_08_10, a_VA_10_12;
struct acc a_VA_12_14, a_VA_14_16, a_VA_16_18;
struct acc a_VA_18_20, a_VA_20;
struct acc a_VD_00_02, a_VD_02_04, a_VD_04_06;
struct acc a_VD_06_08, a_VD_08_10, a_VD_10_12;
struct acc a_VD_12_14, a_VD_14_16, a_VD_16_18;
struct acc a_VD_18_20, a_VD_20;

struct acc a_ng_m00_m02, a_ng_m02_m04, a_ng_m04_m06;
struct acc a_ng_m06_m08, a_ng_m08_m10, a_ng_m10_m12, a_ng_m12_m14;
struct acc a_ng_m14_m16, a_ng_m16_m18, a_ng_m18_m20, a_ng_m20;
struct acc a_s_m00_m02, a_s_m02_m04, a_s_m04_m06;
struct acc a_s_m06_m08, a_s_m08_m10, a_s_m10_m12, a_s_m12_m14;
struct acc a_s_m14_m16, a_s_m16_m18, a_s_m18_m20, a_s_m20;
struct acc a_q_m00_m02, a_q_m02_m04, a_q_m04_m06;
struct acc a_q_m06_m08, a_q_m08_m10, a_q_m10_m12, a_q_m12_m14;
struct acc a_q_m14_m16, a_q_m16_m18, a_q_m18_m20, a_q_m20;
struct acc a_ha_m00_m02, a_ha_m02_m04, a_ha_m04_m06;
struct acc a_ha_m06_m08, a_ha_m08_m10, a_ha_m10_m12, a_ha_m12_m14;
struct acc a_ha_m14_m16, a_ha_m16_m18, a_ha_m18_m20, a_ha_m20;
struct acc a_VA_m00_m02, a_VA_m02_m04, a_VA_m04_m06;
struct acc a_VA_m06_m08, a_VA_m08_m10, a_VA_m10_m12;
struct acc a_VA_m12_m14, a_VA_m14_m16, a_VA_m16_m18;
struct acc a_VA_m18_m20, a_VA_m20;
struct acc a_VD_m00_m02, a_VD_m02_m04, a_VD_m04_m06;
struct acc a_VD_m06_m08, a_VD_m08_m10, a_VD_m10_m12;
struct acc a_VD_m12_m14, a_VD_m14_m16, a_VD_m16_m18;
struct acc a_VD_m18_m20, a_VD_m20;

struct acc Ave_a_ng_00_02, Ave_a_ng_02_04, Ave_a_ng_04_06;
struct acc Ave_a_ng_06_08, Ave_a_ng_08_10, Ave_a_ng_10_12, Ave_a_ng_12_14;
struct acc Ave_a_ng_14_16, Ave_a_ng_16_18, Ave_a_ng_18_20, Ave_a_ng_20;
struct acc Ave_a_s_00_02, Ave_a_s_02_04, Ave_a_s_04_06;
struct acc Ave_a_s_06_08, Ave_a_s_08_10, Ave_a_s_10_12, Ave_a_s_12_14;
struct acc Ave_a_s_14_16, Ave_a_s_16_18, Ave_a_s_18_20, Ave_a_s_20;
struct acc Ave_a_q_00_02, Ave_a_q_02_04, Ave_a_q_04_06;
struct acc Ave_a_q_06_08, Ave_a_q_08_10, Ave_a_q_10_12, Ave_a_q_12_14;
struct acc Ave_a_q_14_16, Ave_a_q_16_18, Ave_a_q_18_20, Ave_a_q_20;
struct acc Ave_a_ha_00_02, Ave_a_ha_02_04, Ave_a_ha_04_06;
struct acc Ave_a_ha_06_08, Ave_a_ha_08_10, Ave_a_ha_10_12, Ave_a_ha_12_14;
struct acc Ave_a_ha_14_16, Ave_a_ha_16_18, Ave_a_ha_18_20, Ave_a_ha_20;
struct acc Ave_a_VA_00_02, Ave_a_VA_02_04, Ave_a_VA_04_06;
struct acc Ave_a_VA_06_08, Ave_a_VA_08_10, Ave_a_VA_10_12;
struct acc Ave_a_VA_12_14, Ave_a_VA_14_16, Ave_a_VA_16_18;
struct acc Ave_a_VA_18_20, Ave_a_VA_20; 
struct acc Ave_a_VD_00_02, Ave_a_VD_02_04, Ave_a_VD_04_06;
struct acc Ave_a_VD_06_08, Ave_a_VD_08_10, Ave_a_VD_10_12;
struct acc Ave_a_VD_12_14, Ave_a_VD_14_16, Ave_a_VD_16_18;
struct acc Ave_a_VD_18_20, Ave_a_VD_20;

struct acc Ave_a_ng_m00_m02, Ave_a_ng_m02_m04, Ave_a_ng_m04_m06;
struct acc Ave_a_ng_m06_m08, Ave_a_ng_m08_m10, Ave_a_ng_m10_m12, Ave_a_ng_m12_m14;
struct acc Ave_a_ng_m14_m16, Ave_a_ng_m16_m18, Ave_a_ng_m18_m20, Ave_a_ng_m20;
struct acc Ave_a_s_m00_m02, Ave_a_s_m02_m04, Ave_a_s_m04_m06;
struct acc Ave_a_s_m06_m08, Ave_a_s_m08_m10, Ave_a_s_m10_m12, Ave_a_s_m12_m14;
struct acc Ave_a_s_m14_m16, Ave_a_s_m16_m18, Ave_a_s_m18_m20, Ave_a_s_m20;
struct acc Ave_a_q_m00_m02, Ave_a_q_m02_m04, Ave_a_q_m04_m06;
struct acc Ave_a_q_m06_m08, Ave_a_q_m08_m10, Ave_a_q_m10_m12, Ave_a_q_m12_m14;
struct acc Ave_a_q_m14_m16, Ave_a_q_m16_m18, Ave_a_q_m18_m20, Ave_a_q_m20;
struct acc Ave_a_ha_m00_m02, Ave_a_ha_m02_m04, Ave_a_ha_m04_m06;
struct acc Ave_a_ha_m06_m08, Ave_a_ha_m08_m10, Ave_a_ha_m10_m12, Ave_a_ha_m12_m14;
struct acc Ave_a_ha_m14_m16, Ave_a_ha_m16_m18, Ave_a_ha_m18_m20, Ave_a_ha_m20;
struct acc Ave_a_VA_m00_m02, Ave_a_VA_m02_m04, Ave_a_VA_m04_m06;
struct acc Ave_a_VA_m06_m08, Ave_a_VA_m08_m10, Ave_a_VA_m10_m12;
struct acc Ave_a_VA_m12_m14, Ave_a_VA_m14_m16, Ave_a_VA_m16_m18;
struct acc Ave_a_VA_m18_m20, Ave_a_VA_m20;
struct acc Ave_a_VD_m00_m02, Ave_a_VD_m02_m04, Ave_a_VD_m04_m06;
struct acc Ave_a_VD_m06_m08, Ave_a_VD_m08_m10, Ave_a_VD_m10_m12;
struct acc Ave_a_VD_m12_m14, Ave_a_VD_m14_m16, Ave_a_VD_m16_m18;
struct acc Ave_a_VD_m18_m20, Ave_a_VD_m20;

/* ********************** Distribution on ha *********************** */

struct acc ha_ng_00_02, ha_ng_02_04, ha_ng_04_06, ha_ng_06_08, ha_ng_08_10;
struct acc ha_q_00_02, ha_q_02_04, ha_q_04_06, ha_q_06_08, ha_q_08_10;
struct acc ha_hs_00_02, ha_hs_02_04, ha_hs_04_06, ha_hs_06_08, ha_hs_08_10;
struct acc ha_VA_00_02, ha_VA_02_04, ha_VA_04_06, ha_VA_06_08, ha_VA_08_10;
struct acc ha_VD_00_02, ha_VD_02_04, ha_VD_04_06, ha_VD_06_08, ha_VD_08_10;

struct acc Ave_ha_ng_00_02, Ave_ha_ng_02_04, Ave_ha_ng_04_06;
struct acc Ave_ha_ng_06_08, Ave_ha_ng_08_10;
struct acc Ave_ha_q_00_02, Ave_ha_q_02_04, Ave_ha_q_04_06;
struct acc Ave_ha_q_06_08, Ave_ha_q_08_10;
struct acc Ave_ha_hs_00_02, Ave_ha_hs_02_04, Ave_ha_hs_04_06;
struct acc Ave_ha_hs_06_08, Ave_ha_hs_08_10;
struct acc Ave_ha_VA_00_02, Ave_ha_VA_02_04, Ave_ha_VA_04_06;
struct acc Ave_ha_VA_06_08, Ave_ha_VA_08_10;
struct acc Ave_ha_VD_00_02, Ave_ha_VD_02_04, Ave_ha_VD_04_06;
struct acc Ave_ha_VD_06_08, Ave_ha_VD_08_10;

/* ********************** Distribution on s *********************** */

struct acc s_ng_0, s_ng_0_106, s_ng_106_104, s_ng_104_102, s_ng_102_101;
struct acc s_ng_01_02, s_ng_02_04, s_ng_04_06, s_ng_06_10, s_ng_10;
struct acc s_q_0, s_q_0_106, s_q_106_104, s_q_104_102, s_q_102_101;
struct acc s_q_01_02, s_q_02_04, s_q_04_06, s_q_06_10, s_q_10;
struct acc s_hs_0_106, s_hs_106_104, s_hs_104_102, s_hs_102_101;
struct acc s_hs_0, s_hs_01_02, s_hs_02_04, s_hs_04_06, s_hs_06_10, s_hs_10;
struct acc s_VA_0, s_VA_0_106, s_VA_106_104, s_VA_104_102, s_VA_102_101;
struct acc s_VA_01_02, s_VA_02_04, s_VA_04_06, s_VA_06_10, s_VA_10;
struct acc s_VD_0, s_VD_0_106, s_VD_106_104, s_VD_104_102, s_VD_102_101;
struct acc s_VD_01_02, s_VD_02_04, s_VD_04_06, s_VD_06_10, s_VD_10;

struct acc Ave_s_ng_0, Ave_s_ng_0_106, Ave_s_ng_106_104, Ave_s_ng_104_102, Ave_s_ng_102_101;
struct acc Ave_s_ng_01_02, Ave_s_ng_02_04, Ave_s_ng_04_06, Ave_s_ng_06_10, Ave_s_ng_10;
struct acc Ave_s_q_0, Ave_s_q_0_106, Ave_s_q_106_104, Ave_s_q_104_102, Ave_s_q_102_101;
struct acc Ave_s_q_01_02, Ave_s_q_02_04, Ave_s_q_04_06, Ave_s_q_06_10, Ave_s_q_10;
struct acc Ave_s_hs_0, Ave_s_hs_0_106, Ave_s_hs_106_104, Ave_s_hs_104_102, Ave_s_hs_102_101;
struct acc Ave_s_hs_01_02, Ave_s_hs_02_04, Ave_s_hs_04_06, Ave_s_hs_06_10, Ave_s_hs_10;
struct acc Ave_s_VA_0, Ave_s_VA_0_106, Ave_s_VA_106_104, Ave_s_VA_104_102, Ave_s_VA_102_101;
struct acc Ave_s_VA_01_02, Ave_s_VA_02_04, Ave_s_VA_04_06, Ave_s_VA_06_10, Ave_s_VA_10;
struct acc Ave_s_VD_0, Ave_s_VD_0_106, Ave_s_VD_106_104, Ave_s_VD_104_102, Ave_s_VD_102_101;
struct acc Ave_s_VD_01_02, Ave_s_VD_02_04, Ave_s_VD_04_06, Ave_s_VD_06_10, Ave_s_VD_10;

/* ************************************************************** */

struct acc SK2[NN], num_SO[NN], b_SS[NN];
struct acc gmean_a[NN], gmean_s[NN], gvar_s[NN], gvar_a[NN], C2[NN], rA[NN];
struct acc VA_s[NN], VD_s[NN];
struct acc rAB[NN];

void initialize(), phenotypeB(), neutral_genes(), disorder_NoSS();
void recombination_masks(), genotypic_values(), regression();
void mutation_neutral(), mutation_selected();
void selected_genes(), mating_binomial(), writelastgeneration();
void fprintmse();

FILE *fptr, *fgen, *fdat, *fpop, *fopen(), *fgens;


/* ***************************************************************** */


main()
{
    fptr = fopen ("dfilename","w");
    fgen = fopen ("genfile","w");
    fpop = fopen ("popfile","w");
    fdat = fopen ("datafile","w");

    getinputs();
    headings();
    recombination_masks(RM);
    initialize (gm,initialgen);
    dumpinitialpopulation();

    for (gen=0; gen<generations; gen++)
    {
	fgens = fopen ("generation.dat","a");
	fprintf (fgens,"generation %d\n", gen);
	fclose(fgens);

	if (tracelevel!=0)     fprintf (fptr,"\n\n ***** generation %d **** \n\n", gen);
	genotypic_values (gm, s, hs, a, ha, genvalm_s, genvalm_a, RM);
               phenotypeB (genvalm_s, genvalm_a, pm_s, pm_a);
	dumpphenotypes();
	regression (pm_s, pm_a);
      	mating_binomial_fecundity (gm, sm, pm_s);
	dumpoffspring();
	neutral_genes (gm, q, initialgen, cumH);
	selected_genes (gm, q, s, hs, a, ha, initialgen, cumH);
	mutation_neutral (gm, q, initialgen, cumH);
	mutation_selected (gm, q, s, hs, a, ha, initialgen, cumH);
	dumpoffspringaftermutation();
	if (gen%block==block-1)
	{
		printout();
		settozero_blocks_on_q();
		settozero_blocks_on_a();
		settozero_blocks_on_ha();
		settozero_blocks_on_s();
	}
    }
    writelastgeneration (gm, s, a, hs, ha, initialgen);
    writeseed();
}


/* ***************************************************************** */


getinputs()
{
	tracestart();
	getseed();
	getintandskip("Number of individuals (max 10000) :",&NIND,2,10000);
	getrealandskip("probability of selfing (Random:99):",&PS,0.0,99.0);
	getrealandskip("Length of genome in Morgans (99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("Number of pieces of chromosome (min 1, max 2000) :",&NCRO,1,2000);
	getintandskip("Number of loci per piece (first is neutral) (min 2, max 30) :",&NLOCI,2,30);
	TOTLOCI = NCRO * NLOCI;
	getrealandskip("Mutation rate per haploid genome for QT :",&Lambda_a,0.0,(double)infinity);
	getrealandskip("Lethal mutation rate (s=1,h=0.02):",&Lambda_L,0.0,(double)infinity);
	getrealandskip("Absolute effect of lethals for QT in heter. Normal (aL,aL/2)(haL=1) :",&ave_aL,0.0,(double)infinity);
	getintandskip ("Random proportion (0) or only large mutants (1) :", &whatproportion,0,1);
	 if(whatproportion==1) getrealandskip ("Values of QTL above threshold :", &threshold,0.0,(double)infinity);
	else		getrealandskip("Proportion of QT-mutants affecting fitness :",&Psi,0.0,1.0);

	getrealandskip("Beta_s :",&beta_s,0.0,(double)infinity);
	getrealandskip("Beta_a :",&beta_a,0.0,(double)infinity);

	getrealandskip("Average |s| :",&ave_s,0.0,1.0);
	alpha_s = beta_s / ave_s;
	epsilon_s = sqrt(beta_s*(1.0+beta_s)/alpha_s/alpha_s);

	getrealandskip("Average |a| :",&ave_a,0.0,(double)infinity);
	alpha_a = beta_a / ave_a;
	epsilon_a = sqrt(beta_a*(1.0+beta_a)/alpha_a/alpha_a);

	getrealandskip("Proportion of positives for fitness :",&P_s,0.0,1.0);
	getrealandskip("Proportion of positives for the trait :",&P_a,0.0,1.0);

	getintandskip("dom. model (0:constant 1:Deng 2:CK94(gamma):", &dom,0,2);

	getrealandskip("hs (mod 0) or k_s (mod 1):", &k_s,0.0,(double)infinity);
	getrealandskip("Average h (mod 2):",&ave_hs,0.0,(double)infinity);
	if (dom==2)	k_s = alpha_s * (pow((2.0*ave_hs),((-1.0)/beta_s))-1.0);

	getrealandskip("ha (mod 0) or k_a (mod 1):", &k_a,0.0,(double)infinity);
	getrealandskip("Average h (mod 2):",&ave_ha,0.0,(double)infinity);
	if (dom==2)	k_a = alpha_a * (pow((2.0*ave_ha),((-1.0)/beta_a))-1.0);

	getrealandskip("rho (99: a=s):",&rho,-1.0,99.0);
	getrealandskip("Stabilizing selection (Vs) :",&Vs,0.0,(double)infinity);
	VE=1.0;
	getintandskip("Model of fitness selection (1:mul; 2: add) :",&model_selection,1,2);
	getintandskip("Number of generations :",&generations,2,infinity);
	getintandskip("Number of generations per block :",&block,1,infinity);
	getintandskip("Generation bottleneck (N=10,t=1; N=10,t=10; N=10,t=20; N=100,t=100) 0,1,10,20,100 :",&GENBOT,0,100);
}


/* **************************************************************** */


headings()
{
	fgen = fopen ("genfile","a");
	fprintf(fgen,"\n\n N=%d    PS=%4.2f    L=%f   N.S-LOCI=%d    N.N-LOCI=%d    gens=%d\n Lambda_a=%6.4f   Lambda_L=%6.4f   ave_aL=%6.4f    Psi=%6.4f   threshold=%6.4f   rho=%6.4f   Vs=%6.4f \n ave_s=%6.4f   PP_s=%6.4f   beta_s=%6.4f   alpha_s=%6.4f  dom=%d   k_s=%6.4f   ave_hs=%6.4f\n ave_a=%6.4f   PP_a=%6.4f   beta_a=%6.4f   alpha_a=%6.4f   k_a=%6.4f   ave_ha=%6.4f\n", NIND, PS, L, TOTLOCI-NCRO, NCRO, generations, Lambda_a, Lambda_L, ave_aL, Psi, threshold, rho, Vs, ave_s, P_s, beta_s, alpha_s, dom, k_s, ave_hs, ave_a, P_a, beta_a, alpha_a, k_a, ave_ha);
    	fclose(fgen);
}


/* ***************************************************************** */


void recombination_masks (RM)
int RM[];
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}


/* ***************************************************************** */


void initialize (gm, initialgen)
int gm[][MMM][2], initialgen[][31];
{
    gen = (-1);

    if(model_selection==1)	    addedfixed_s=1.0;
    else			    addedfixed_s=0.0;
    addedfixed_a=0.0;


    /* NEUTRAL LOCI WITH INITIAL FREQUENCY 0.5 */

    for (k=0; k<NCRO; k++)
    {
	    for (i=0; i<NIND; i++)
	    {
	    	if (uniform() < 0.5)	gm[i][k][0]=(gm[i][k][0] | RM[0]);
	    	if (uniform() < 0.5)	gm[i][k][1]=(gm[i][k][1] | RM[0]);
		initialgen[k][0] = 0;
	    }
    }
    
    /* SELECTED LOCI WITH POISSON (2NL) NEW MUTATIONS */

    if ( (exp(-2.0*(double)NIND*Lambda_a) != 0.0)&&
			(2.0*(double)NIND*Lambda_a < normalthreshold) )
    generatepoissontable(2.0*(double)NIND*Lambda_a, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

    if ( (exp(-2.0*(double)NIND*Lambda_L) != 0.0)&&
			(2.0*(double)NIND*Lambda_L < normalthreshold) )
    generatepoissontable(2.0*(double)NIND*Lambda_L, &lastinmutantspoissontableL, mutantspoissontableL, maxmpt-1);

    for (k=0; k<NCRO; k++)
    {
	for (l=1; l<NLOCI; l++)
	{
	    initialgen[k][l] = (-99);
	    q[k][l] = 0.0;
	    s[k][l] = 0.0; hs[k][l] = 0.0;
	    a[k][l] = 0.0; ha[k][l] = 0.0;
	}
    }

    mutation_selected (gm,q,s,hs,a,ha,initialgen,cumH);
    selected_genes (gm,q,s,hs,a,ha,initialgen,cumH);
}


/* **************************************************************** */


dumpinitialpopulation()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n*****************************************************\n");
	fprintf(fptr,"\n********* Initial population (gm0 gm1) ************\n");
	for (i=0; i<NIND; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}


/* **************************************************************** */


void genotypic_values (gm,s,hs,a,ha,genvalm_s,genvalm_a,RM)
int gm[][MMM][2], RM[];
double genvalm_s[], s[][31], hs[][31];
double genvalm_a[], a[][31], ha[][31];
{
	for (i=0; i<NIND; i++)
	{
		genvalm_s[i]=addedfixed_s;
		genvalm_a[i]=addedfixed_a;

		for (k=0; k<NCRO; k++)
		{
		   for (l=1; l<NLOCI; l++)
		   {
		        if (initialgen[k][l] != (-99))
		        {
	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  				{
				if (model_selection==1)    genvalm_s[i] *= (1.0 + s[k][l]);
				else    genvalm_s[i] += s[k][l];
				genvalm_a[i] += a[k][l];
			}
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 								/* AA */;
			else
			{
			 	if (model_selection==1)     genvalm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
				else    genvalm_s[i] += s[k][l]*hs[k][l];
				genvalm_a[i] += a[k][l]*ha[k][l];
			}
		        }
		   }
		}
		if (tracelevel!=0)	fprintf(fptr," %d    genvalm_a = %f      genvalm_s = %f\n",
			i, genvalm_a[i], genvalm_s[i]);
	}
}


/* **************************************************************** */


void phenotypeB (genvalm_s, genvalm_a, pm_s, pm_a)
double pm_s[], pm_a[];
double genvalm_s[], genvalm_a[];
{
	int ii, it;
	double gsum_s=0.0, gsum2_s=0.0, gsum_sa=0.0;
	double gsum_a=0.0, gsum2_a=0.0;
	double maxpm=0.0, sq_pm_a=0.0;

	for (i=0; i<NIND; i++)
	{	
	    pm_a[i] = genvalm_a[i] + normal(0.0, sqrt(VE));
	    gsum_a += genvalm_a[i];
	    gsum2_a += (genvalm_a[i]* genvalm_a[i]);
	}
	accum (&gmean_a[gen/block], gsum_a/(double)NIND);
	accum (&gvar_a[gen/block], (gsum2_a - (gsum_a*gsum_a / (double)NIND)) / (double)NIND);

	for (i=0; i<NIND; i++)
	{
	    if(Vs==0.0)    pm_s[i]=genvalm_s[i];
	    else
	    {
		sq_pm_a = pow (pm_a[i], 2.0);

	    	if (model_selection==1)    pm_s[i]=genvalm_s[i] * exp(-0.5*(sq_pm_a)/Vs);
	    	else 			pm_s[i]=genvalm_s[i] + exp(-0.5*(sq_pm_a)/Vs);
	    }
	    gsum_s += pm_s[i];
	    gsum2_s += (pm_s[i]*pm_s[i]);
	    gsum_sa += (pm_s[i]*genvalm_a[i]);
	}

	accum (&gmean_s[gen/block], gsum_s/(double)NIND);
	accum (&gvar_s[gen/block], (gsum2_s - (gsum_s*gsum_s / (double)NIND)) / (double)NIND);
	accum (&C2[gen/block], ( (gsum2_s*(double)NIND) / (gsum_s*gsum_s) ) - 1.0);
	accum (&rA[gen/block], (gsum_sa - (gsum_s*gsum_a / (double)NIND))
				/ sqrt(gsum2_a - (gsum_a*gsum_a / (double)NIND)) 
					/ sqrt (gsum2_s - (gsum_s*gsum_s / (double)NIND)));

	if (tracelevel!=0)   fprintf(fptr,"\ngmean_s = %f  gvar_s = %f C2 = %f\n",
		gsum_s/(double)NIND,
			(gsum2_s - (gsum_s*gsum_s / (double)NIND)) / (double)NIND,
				( (gsum2_s*(double)NIND) / (gsum_s*gsum_s) ) - 1.0);

	if (tracelevel!=0)   fprintf(fptr,"\ngmean_a = %f  gvar_a = %f\n",
		gsum_a/(double)NIND,
			(gsum2_a - (gsum_a*gsum_a / (double)NIND)) / ((double)NIND-1.0));

	/* find individual with maximal fitness */
	ii=0;
	for (it=1; it<NIND; it++)		if (pm_s[ii] < pm_s[it])	ii=it;
	maxpm =pm_s[ii];

	for (i=0; i<NIND; i++)	    pm_s[i] = pm_s[i] / maxpm;
}


/* **************************************************************** */


dumpphenotypes()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Phenotypic value and Relative fitnesses \n");
	for (i=0; i<NIND; i++)   fprintf(fptr,"%f    %f\n", pm_a[i], pm_s[i]);
}


/* ***************************************************************** */


void regression (pm_s, pm_a)
double pm_s[], pm_a[];
{
	double phsum_a=0.0, sum_pm_a2=0.0, sum_pm_a4=0.0, sum_pm_s=0.0, sum_pm_a2s=0.0;

	for (i=0; i<NIND; i++)	phsum_a += pm_a[i];

	for (i=0; i<NIND; i++)
	{
		sum_pm_a2 += pow(pm_a[i]-(phsum_a/(double)NIND), 2.0);
		sum_pm_a4 += pow(pm_a[i]-(phsum_a/(double)NIND), 4.0);
		sum_pm_s += pm_s[i];
		sum_pm_a2s += pow(pm_a[i]-(phsum_a/(double)NIND), 2.0) * pm_s[i];
	}
	bSS = ( sum_pm_a2s - (sum_pm_a2 * sum_pm_s / (double)NIND) ) / 
		( sum_pm_a4 - (sum_pm_a2 * sum_pm_a2 / (double)NIND) );

	accum (&b_SS[gen/block], bSS);

	if (tracelevel!=0)	fprintf(fptr," Regression of fitness on squa phen dev from current mean\n");
	if (tracelevel!=0)	fprintf(fptr,"S_a2 = %f   S_a4 = %f   S_s = %f   S_a2s = %f   b_SS = %f\n",
			sum_pm_a2, sum_pm_a4, sum_pm_s, sum_pm_a2s, b_SS[gen/block]);
}


/* ***************************************************************** */

void mating_binomial_fecundity (gm, sm, pm_s)
int gm[][MMM][2], sm[][MMM][2];
double pm_s[];
{
	int EE[MMM], FF[MMM], markerrecom, p1, p2, family[NN];
	double family_sum2=0.0, sq_pm_a, cum=0.0, cum_pm_s[NN], rnd;
		
	for (i=0; i<NIND; i++)
	{
	    for (k=0; k<NCRO; k++)
	    {
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];
	    }
	}

	for (i=0; i<NIND;i++)   family[i]=0;

	/********* Accumulate ps values ***********/

	for (i=0; i<NIND; i++)
	{
		cum += pm_s[i];
		cum_pm_s[i] = cum;
	}
	for (i=0; i<NIND; i++)	cum_pm_s[i] = cum_pm_s[i] / cum;

	if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf (fptr," %d    cum_pm_s = %f\n", i, cum_pm_s[i]);
		
	/******************************************/

	if (tracelevel!=0)	fprintf(fptr,"\n Parents \n");

	for (i=0; i<NIND; i++)
	{
		// ************* BOTTLENECK ***********
	    if (GENBOT==0)		p1 = binarysearchcumulative(cum_pm_s, NIND);
	    else if (GENBOT==1)
	    {
		if (gen == 10000)			p1 = binarysearchcumulative(cum_pm_s, 10);
		else					p1 = binarysearchcumulative(cum_pm_s, NIND);
	    }
	    else if (GENBOT==10)
	    {
		if ((gen >= 9991) && (gen <= 10000)) 	p1 = binarysearchcumulative(cum_pm_s, 10);
		else					p1 = binarysearchcumulative(cum_pm_s, NIND);
	    }
	    else if (GENBOT==20)
	    {
		if ((gen >= 9981) && (gen <= 10000)) 	p1 = binarysearchcumulative(cum_pm_s, 10);
		else					p1 = binarysearchcumulative(cum_pm_s, NIND);
	    }
	    else if (GENBOT==100)
	    {
		if ((gen >= 9901) && (gen <= 10000)) 	p1 = binarysearchcumulative(cum_pm_s, 100);
		else					p1 = binarysearchcumulative(cum_pm_s, NIND);
	    }

	    if (PS==99)
	    {
		// ************* BOTTLENECK ***********
		    if (GENBOT==0)		p2 = binarysearchcumulative(cum_pm_s, NIND);
		    else if (GENBOT==1)
		    {
			if (gen == 10000)			p2 = binarysearchcumulative(cum_pm_s, 10);
			else					p2 = binarysearchcumulative(cum_pm_s, NIND);
		    }
		    else if (GENBOT==10)
		    {
			if ((gen >= 9991) && (gen <= 10000)) 	p2 = binarysearchcumulative(cum_pm_s, 10);
			else					p2 = binarysearchcumulative(cum_pm_s, NIND);
		    }
		    else if (GENBOT==20)
		    {
			if ((gen >= 9981) && (gen <= 10000)) 	p2 = binarysearchcumulative(cum_pm_s, 10);
			else					p2 = binarysearchcumulative(cum_pm_s, NIND);
		    }
		    else if (GENBOT==100)
		    {
			if ((gen >= 9901) && (gen <= 10000)) 	p2 = binarysearchcumulative(cum_pm_s, 100);
			else					p2 = binarysearchcumulative(cum_pm_s, NIND);
		    }
	    }
	    else
	    {
		if (uniform() < PS)		p2=p1;
		else
		{
	    	    do { p2 = binarysearchcumulative(cum_pm_s, NIND); }
	    	    while (p2==p1);
		}
	    }

	    if (tracelevel!=0)   fprintf (fptr,"%d\t%d\n", p1, p2);

	    if (p1==p2)		accum (&num_SO[gen/block], (double)NIND);
	    else		accum (&num_SO[gen/block], 0.0);

	    family[p1]+=1;    family[p2]+=1;

	    if(L==99.0)
	    {	    /* ******************* Free recombination ******************* */

		for (k=0; k<NCRO; k++)
		{
		   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
		   	FF[k] = ~EE[k];
		   	gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
		}
		/* if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);
		*/
		for (k=0; k<NCRO; k++)
		{
		   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
		   	FF[k] = ~EE[k];
		   	gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
		}
		/* if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);
		*/
	    }
	    else
	    {	    /* ************** Restricted recombination ***************** */

		/* ****** Chromosome from father ****** */

		marker = 1;

		for (k=0; k<NCRO; k++)
		{
		      EE[k]=0;

		      for (l=0; l<NLOCI; l++)
		      {
				if (marker==1)
				{
			      if (uniform() < (L / (double)TOTLOCI))
			      {
					EE[k] = EE[k] | RM[l];
					numberofrecom[(gen+1)/block] += 1.0;
					marker = marker * (-1);
			      }
				}
				else
				{
			      if (uniform() < (L / (double)TOTLOCI))
			      {
					numberofrecom[(gen+1)/block] += 1.0;
					marker = marker * (-1);
			      }
			      else   EE[k] = EE[k] | RM[l];
				}
		   	   }
		}
		if (tracelevel!=0)   fprintf (fptr,"\n EE[0]=%d EE[1]=%d \n", EE[0], EE[1]);

		for (k=0; k<NCRO; k++)
		{
			if (uniform() < 0.5) 	EE[k] = ~EE[k];
			FF[k] = ~EE[k];
			gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
		}

		/* ****** Chromosome from mother ****** */

		marker = 1;

		for (k=0; k<NCRO; k++)
		{
		      EE[k]=0;

		      for (l=0; l<NLOCI; l++)
		      {
			if (marker==1)
			{
			      if (uniform() < (L / (double)TOTLOCI))
			      {
				EE[k] = EE[k] | RM[l];
				numberofrecom[(gen+1)/block] += 1.0;
				marker = marker * (-1);
			      }
			}
			else
			{
			      if (uniform() < (L / (double)TOTLOCI))
			      {
				numberofrecom[(gen+1)/block] += 1.0;
				marker = marker * (-1);
			      }
			      else   EE[k] = EE[k] | RM[l];
			}
		      }
		}
		if (tracelevel!=0)   fprintf (fptr,"\n EE[0]=%d EE[1]=%d \n", EE[0], EE[1]);

		for (k=0; k<NCRO; k++)
		{
			if (uniform() < 0.5) 	EE[k] = ~EE[k];
			FF[k] = ~EE[k];
			gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
		}
	    }

		/* *******************PPPPPPPPPPPPPP******************** */

		genvalm_s[i]=addedfixed_s;
		genvalm_a[i]=addedfixed_a;

		for (k=0; k<NCRO; k++)
		{
		   for (l=1; l<NLOCI; l++)
		   {
		        if (initialgen[k][l] != (-99))
		        {
		    		if (((gm[i][k][0] & RM[l])==RM[l])
					&&((gm[i][k][1] & RM[l])==RM[l]))  			
				{
					if (model_selection==1)
						genvalm_s[i] *= (1.0 + s[k][l]);
					else    genvalm_s[i] += s[k][l];
					genvalm_a[i] += a[k][l];
				}
				else    if (((gm[i][k][0] & RM[l])!=RM[l])
						&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
				else
				{
				 	if (model_selection==1)
						genvalm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
					else    genvalm_s[i] += s[k][l]*hs[k][l];
					genvalm_a[i] += a[k][l]*ha[k][l];
				}
		        }
		   }
		}
		if (tracelevel!=0)
			fprintf (fptr," %d    genvalm_a = %f      genvalm_s = %f\n",
			i, genvalm_a[i], genvalm_s[i]);

		if(Vs==0.0)    pm_s[i]=genvalm_s[i];
		else
		{
			sq_pm_a = pow (pm_a[i], 2.0);

			if (model_selection==1)	pm_s[i]=genvalm_s[i] * exp(-0.5*(sq_pm_a)/Vs);
	    		else 			pm_s[i]=genvalm_s[i] + exp(-0.5*(sq_pm_a)/Vs);
		}
	}

	/* ***** VARIANCE OF FAMILY SIZE ***** */

	for (i=0; i<NIND; i++)	family_sum2 += (family[i] * family[i]);
	
	accum(&SK2[gen/block], ((family_sum2-(4.0*NIND))/(NIND-1.0)) );
	if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", ((family_sum2-(4.0*NIND))/(NIND-1.0)));
}


/* **************************************************************** */


dumpoffspring()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Offspring before mutation (gm0 gm1)\n");	
	for (i=0; i<NIND; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}


/* ***************************************************************** */


void neutral_genes (gm,q,initialgen,cumH)
int gm[][MMM][2], initialgen[][31];
double q[][31], cumH[][31];
{
    if (tracelevel!=0)
	fprintf(fptr,"\n neutral genes\npiece\tlocus\tAA\tAa\taa\tfreqAa\t            q\t     inigen\t    cumH\n");

    for (k=0; k<NCRO; k++)
    {
	    AA=0.0; Aa=0.0; aa=0.0;

	    for (i=0; i<NIND; i++)
	    {
	        if (((gm[i][k][0] & RM[0])==RM[0])&&((gm[i][k][1] & RM[0])==RM[0]))          aa+=1.0;
	        else    if (((gm[i][k][0] & RM[0])!=RM[0])&&((gm[i][k][1] & RM[0])!=RM[0])) AA+=1.0;
	        else								            Aa+=1.0;
	    }

	    q[k][0] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
	    cumH[k][0] += Aa/(double)NIND;
	    freqn[(gen+1)/block] += q[k][0];

	    if (q[k][0]==0.0)
	    {
		lostn[(gen+1)/block] += 1.0;
		CHn[(gen+1)/block] += cumH[k][0];
		timetolossn[(gen+1)/block] += (double)((gen+1)-initialgen[k][0]);
	    }
	    else if (q[k][0]==1.0)
	    {
		fixn[(gen+1)/block] += 1.0;
		CHn[(gen+1)/block] += cumH[k][0];
		timetofixn[(gen+1)/block] += (double)((gen+1)-initialgen[k][0]);
		/* make it wild type */
		for (i=0; i<NIND; i++)
		{
		    gm[i][k][0] = ( gm[i][k][0] & (~RM[0]) );
		    gm[i][k][1] = ( gm[i][k][1] & (~RM[0]) );
		}
		cumH[k][0]=0.0;
		q[k][0]=0.0;
	    }
	    else
	    {
		othern[(gen+1)/block] += 1.0;
	    }

	    if (tracelevel!=0)    fprintf(fptr,"%d\t0\t%1.0f\t%1.0f\t%1.0f\t%f\t%f\t%d\t%f\n"
		,k,AA,Aa,aa,Aa/(double)NIND,q[k][0],initialgen[k][0],cumH[k][0]);
}
}


/* ***************************************************************** */


void selected_genes (gm,q,s,hs,a,ha,initialgen,cumH)
int gm[][MMM][2], initialgen[][31];
double q[][31], s[][31], hs[][31], a[][31], ha[][31], cumH[][31];
{
    double VAs=0.0, VDs=0.0, AB=0.0, ab=0.0, Ab=0.0, aB=0.0;

    if (tracelevel!=0)
	fprintf(fptr,"\n selected genes\npiece\tlocus\tAA\tAa\taa\tfreqAa\t            q\t     inigen\t    cumH\t    d_a\t    d_s\t    alfa_a\t    alfa_s\t    VAa\t    VAs\n");

    for (k=0; k<NCRO; k++)
    {
	for (l=1; l<NLOCI; l++)
	{
	    AA=0.0; Aa=0.0; aa=0.0;

	    for (i=0; i<NIND; i++)
	    {
	         if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))           aa+=1.0;
	         else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))   AA+=1.0;
	         else								           Aa+=1.0;
	    }

	    q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
	    cumH[k][l] += Aa/(double)NIND;

	    if (initialgen[k][l]==(-99)) /* gene non segregating */;
	    else if (q[k][l]==0.0)
	    {
		lost[(gen+1)/block] += 1.0;
		CH[(gen+1)/block] += cumH[k][l];
		timetoloss[(gen+1)/block] += (double)((gen+1)-initialgen[k][l]);
		initialgen[k][l]=(-99);
		cumH[k][l]=0.0;
		s[k][l]=0.0; hs[k][l]=0.0;
		a[k][l]=0.0; ha[k][l]=0.0;
	    }
	    else if (q[k][l]==1.0)
	    {
		fix[(gen+1)/block] += 1.0;
		CH[(gen+1)/block] += cumH[k][l];
		timetofix[(gen+1)/block] += (double)((gen+1)-initialgen[k][l]);
		/* make it wild type */
		for (i=0; i<NIND; i++)
		{
		    gm[i][k][0] = ( gm[i][k][0] & (~RM[l]) );
		    gm[i][k][1] = ( gm[i][k][1] & (~RM[l]) );
		}
		initialgen[k][l]=(-99);
		if(model_selection==1)    addedfixed_s *= (1.0 + s[k][l]);
		else 			addedfixed_s += s[k][l];
		addedfixed_a += a[k][l];
		cumH[k][l]=0.0;
		q[k][l]=0.0;
		s[k][l]=0.0; hs[k][l]=0.0;
		a[k][l]=0.0; ha[k][l]=0.0;
	    }
	    else
	    {
		other[(gen+1)/block] += 1.0;
		freq[(gen+1)/block] += q[k][l];

		d_a = (a[k][l]/2.0) * (2.0*ha[k][l] - 1.0);
		if (a[k][l] >= 0.0)		alfa_a = (a[k][l]/2.0) + ( d_a * (1.0 - 2.0*q[k][l]) );
		else			alfa_a = (-a[k][l]/2.0) + ( d_a * (2.0*q[k][l] - 1.0) );
		
		d_s = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
		if (s[k][l] >= 0.0)		alfa_s = (s[k][l]/2.0) + ( d_s * (1.0 - 2.0*q[k][l]) );
		else			alfa_s = (-s[k][l]/2.0) + ( d_s * (2.0*q[k][l] - 1.0) );

		VAs += 2.0 * alfa_s * alfa_s * q[k][l] * (1.0 - q[k][l]);
		VDs +=  pow(2.0 * d_s * q[k][l] * (1.0 - q[k][l]), 2.0);

		distribution_gens_on_q();
		distribution_gens_on_a();
		if(s[k][l]!=(-1.0))	distribution_gens_on_ha();
		distribution_gens_on_s();
	    }

	    if (tracelevel!=0)    fprintf(fptr,"%d\t%d\t%1.0f\t%1.0f\t%1.0f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",k,l,AA,Aa,aa,Aa/(double)NIND,q[k][l],initialgen[k][l],cumH[k][l],d_a,d_s,alfa_a,alfa_s,2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]),2.0 * alfa_s * alfa_s * q[k][l] * (1.0 - q[k][l]));
	}
    }

    accum (&VA_s[(gen+1)/block], VAs);
    accum (&VD_s[(gen+1)/block], VDs);
    distribution_blocks_on_q();
    distribution_blocks_on_a();
    distribution_blocks_on_ha();
    distribution_blocks_on_s();
    settozero_gens_on_q();
    settozero_gens_on_a();
    settozero_gens_on_ha();
    settozero_gens_on_s();

    /* LINKAGE DISEQUILIBRIUM */
    for (k=0; k<NCRO; k++)
    for (l=1; l<NLOCI; l++)
    {
	AB=0.0; ab=0.0; Ab=0.0; aB=0.0;

	if ((q[k][l]>0) && (q[k][l+1]>0))
	{
		for (i=0; i<NIND; i++)
		{
			if ((((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][0] & RM[l+1])!=RM[l+1])) || (((gm[i][k][1] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l+1])!=RM[l+1])))	AB += 1.0;
			else if ((((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][0] & RM[l+1])==RM[l+1])) || (((gm[i][k][1] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l+1])==RM[l+1])))	ab += 1.0;
			else if ((((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][0] & RM[l+1])==RM[l+1])) || (((gm[i][k][1] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l+1])==RM[l+1])))	Ab += 1.0;
			else if ((((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][0] & RM[l+1])!=RM[l+1])) || (((gm[i][k][1] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l+1])!=RM[l+1])))	aB += 1.0;
		}
		accum (&rAB[gen/block],  (((AB/(2.0*(double)NIND))*(ab/(2.0*(double)NIND))) - ((Ab/(2.0*(double)NIND))*(aB/(2.0*(double)NIND))))/sqrt(q[k][l]*(1-q[k][l])*q[k][l+1]*(1-q[k][l+1])) );

		if (tracelevel!=0)    fprintf(fptr,"LD k=%d l=%d l+1=%d AB=%f ab=%f Ab=%f aB=%f ql=%f ql+1=%f	rAB=%f\n",k,l,l+1,(AB/(2.0*(double)NIND)),(ab/(2.0*(double)NIND)),(Ab/(2.0*(double)NIND)),(aB/(2.0*(double)NIND)),q[k][l],q[k][l+1],(((AB/(2.0*(double)NIND))*(ab/(2.0*(double)NIND))) - ((Ab/(2.0*(double)NIND))*(aB/(2.0*(double)NIND))))/sqrt(q[k][l]*(1-q[k][l])*q[k][l+1]*(1-q[k][l+1])));
	}
    }
}


/* ***************************************************************** */


distribution_gens_on_q()
{
	accum (&ng_all,  1.0);
	accum (&s_all,  s[k][l]);
	if (a[k][l] < 0.0)	accum (&a_all,  (-a[k][l]));
	else		accum (&a_all,  a[k][l]);
	if (s[k][l] != 0.0)   accum (&hs_all,  hs[k][l]);
	accum (&ha_all,  ha[k][l]);
	accum (&VA_all,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
	accum (&VD_all,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));

	if ( (q[k][l] > 0.0) && (q[k][l] <= 0.1) )
	{
		accum (&q_ng_00_01,  1.0);
		accum (&q_s_00_01,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_00_01,  (-a[k][l]));
		else		accum (&q_a_00_01,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_00_01,  hs[k][l]);
		accum (&q_ha_00_01,  ha[k][l]);
		accum (&q_VA_00_01,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_00_01,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.1) && (q[k][l] <= 0.2) )
	{
		accum (&q_ng_01_02,  1.0);
		accum (&q_s_01_02,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_01_02,  (-a[k][l]));
		else		accum (&q_a_01_02,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_01_02,  hs[k][l]);
		accum (&q_ha_01_02,  ha[k][l]);
		accum (&q_VA_01_02,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_01_02,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.2) && (q[k][l] <= 0.3) )
	{
		accum (&q_ng_02_03,  1.0);
		accum (&q_s_02_03,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_02_03,  (-a[k][l]));
		else		accum (&q_a_02_03,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_02_03,  hs[k][l]);
		accum (&q_ha_02_03,  ha[k][l]);
		accum (&q_VA_02_03,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_02_03,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.3) && (q[k][l] <= 0.4) )
	{
		accum (&q_ng_03_04,  1.0);
		accum (&q_s_03_04,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_03_04,  (-a[k][l]));
		else		accum (&q_a_03_04,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_03_04,  hs[k][l]);
		accum (&q_ha_03_04,  ha[k][l]);
		accum (&q_VA_03_04,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_03_04,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.4) && (q[k][l] <= 0.5) )
	{
		accum (&q_ng_04_05,  1.0);
		accum (&q_s_04_05,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_04_05,  (-a[k][l]));
		else		accum (&q_a_04_05,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_04_05,  hs[k][l]);
		accum (&q_ha_04_05,  ha[k][l]); 
		accum (&q_VA_04_05,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_04_05,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));

	}
	else if ( (q[k][l] > 0.5) && (q[k][l] <= 0.6) )
	{
		accum (&q_ng_05_06,  1.0);
		accum (&q_s_05_06,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_05_06,  (-a[k][l]));
		else		accum (&q_a_05_06,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_05_06,  hs[k][l]);
		accum (&q_ha_05_06,  ha[k][l]);
		accum (&q_VA_05_06,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_05_06,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.6) && (q[k][l] <= 0.7) )
	{
		accum (&q_ng_06_07,  1.0);
		accum (&q_s_06_07,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_06_07,  (-a[k][l]));
		else		accum (&q_a_06_07,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_06_07,  hs[k][l]);
		accum (&q_ha_06_07,  ha[k][l]);
		accum (&q_VA_06_07,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_06_07,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.7) && (q[k][l] <= 0.8) )
	{
		accum (&q_ng_07_08,  1.0);
		accum (&q_s_07_08,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_07_08,  (-a[k][l]));
		else		accum (&q_a_07_08,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_07_08,  hs[k][l]);
		accum (&q_ha_07_08,  ha[k][l]);
		accum (&q_VA_07_08,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_07_08,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.8) && (q[k][l] <= 0.9) )
	{
		accum (&q_ng_08_09,  1.0);
		accum (&q_s_08_09,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_08_09,  (-a[k][l]));
		else		accum (&q_a_08_09,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_08_09,  hs[k][l]);
		accum (&q_ha_08_09,  ha[k][l]);
		accum (&q_VA_08_09,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_08_09,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (q[k][l] > 0.9) && (q[k][l] < 1.0) )
	{
		accum (&q_ng_09_10,  1.0);
		accum (&q_s_09_10,  s[k][l]);
		if (a[k][l] < 0.0)	accum (&q_a_09_10,  (-a[k][l]));
		else		accum (&q_a_09_10,  a[k][l]);
		if (s[k][l] != 0.0)   accum (&q_hs_09_10,  hs[k][l]);
		accum (&q_ha_09_10,  ha[k][l]);
		accum (&q_VA_09_10,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&q_VD_09_10,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
}


/* **************************************************************** */


distribution_gens_on_a()
{
	/* ***** negative values ***** */

	if ( (a[k][l] < 0.0) && (a[k][l] >= (-0.2)) )
	{
		accum (&a_ng_m00_m02,  1.0);
		accum (&a_s_m00_m02,  s[k][l]);
		accum (&a_q_m00_m02,  q[k][l]);
		accum (&a_ha_m00_m02,  ha[k][l]);
		accum (&a_VA_m00_m02,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m00_m02,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-0.2)) && (a[k][l] > (-0.4)) )
	{
		accum (&a_ng_m02_m04,  1.0);
		accum (&a_s_m02_m04,  s[k][l]);
		accum (&a_q_m02_m04,  q[k][l]);
		accum (&a_ha_m02_m04,  ha[k][l]);
		accum (&a_VA_m02_m04,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m02_m04,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-0.4)) && (a[k][l] > (-0.6)) )
	{
		accum (&a_ng_m04_m06,  1.0);
		accum (&a_s_m04_m06,  s[k][l]);
		accum (&a_q_m04_m06,  q[k][l]);
		accum (&a_ha_m04_m06,  ha[k][l]);
		accum (&a_VA_m04_m06,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m04_m06,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-0.6)) && (a[k][l] > (-0.8)) )
	{
		accum (&a_ng_m06_m08,  1.0);
		accum (&a_s_m06_m08,  s[k][l]);
		accum (&a_q_m06_m08,  q[k][l]);
		accum (&a_ha_m06_m08,  ha[k][l]);
		accum (&a_VA_m06_m08,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m06_m08,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-0.8)) && (a[k][l] > (-1.0)) )
	{
		accum (&a_ng_m08_m10,  1.0);
		accum (&a_s_m08_m10,  s[k][l]);
		accum (&a_q_m08_m10,  q[k][l]);
		accum (&a_ha_m08_m10,  ha[k][l]);
		accum (&a_VA_m08_m10,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m08_m10,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-1.0)) && (a[k][l] > (-1.2)) )
	{
		accum (&a_ng_m10_m12,  1.0);
		accum (&a_s_m10_m12,  s[k][l]);
		accum (&a_q_m10_m12,  q[k][l]);
		accum (&a_ha_m10_m12,  ha[k][l]);
		accum (&a_VA_m10_m12,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m10_m12,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-1.2)) && (a[k][l] > (-1.4)) )
	{
		accum (&a_ng_m12_m14,  1.0);
		accum (&a_s_m12_m14,  s[k][l]);
		accum (&a_q_m12_m14,  q[k][l]);
		accum (&a_ha_m12_m14,  ha[k][l]);
		accum (&a_VA_m12_m14,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m12_m14,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-1.4)) && (a[k][l] > (-1.6)) )
	{
		accum (&a_ng_m14_m16,  1.0);
		accum (&a_s_m14_m16,  s[k][l]);
		accum (&a_q_m14_m16,  q[k][l]);
		accum (&a_ha_m14_m16,  ha[k][l]);
		accum (&a_VA_m14_m16,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m14_m16,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-1.6)) && (a[k][l] > (-1.8)) )
	{
		accum (&a_ng_m16_m18,  1.0);
		accum (&a_s_m16_m18,  s[k][l]);
		accum (&a_q_m16_m18,  q[k][l]);
		accum (&a_ha_m16_m18,  ha[k][l]);
		accum (&a_VA_m16_m18,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m16_m18,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] <= (-1.8)) && (a[k][l] > (-2.0)) )
	{
		accum (&a_ng_m18_m20,  1.0);
		accum (&a_s_m18_m20,  s[k][l]);
		accum (&a_q_m18_m20,  q[k][l]);
		accum (&a_ha_m18_m20,  ha[k][l]);
		accum (&a_VA_m18_m20,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m18_m20,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if (a[k][l] <= (-2.0))
	{
		accum (&a_ng_m20,  1.0);
		accum (&a_s_m20,  s[k][l]);
		accum (&a_q_m20,  q[k][l]);
		accum (&a_ha_m20,  ha[k][l]);
		accum (&a_VA_m20,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_m20,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}

	/* ***** positive values ***** */

	if ( (a[k][l] >= 0.0) && (a[k][l] < 0.2) )
	{
		accum (&a_ng_00_02,  1.0);
		accum (&a_s_00_02,  s[k][l]);
		accum (&a_q_00_02,  q[k][l]);
		accum (&a_ha_00_02,  ha[k][l]);
		accum (&a_VA_00_02,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_00_02,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 0.2) && (a[k][l] < 0.4) )
	{
		accum (&a_ng_02_04,  1.0);
		accum (&a_s_02_04,  s[k][l]);
		accum (&a_q_02_04,  q[k][l]);
		accum (&a_ha_02_04,  ha[k][l]);
		accum (&a_VA_02_04,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_02_04,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 0.4) && (a[k][l] < 0.6) )
	{
		accum (&a_ng_04_06,  1.0);
		accum (&a_s_04_06,  s[k][l]);
		accum (&a_q_04_06,  q[k][l]);
		accum (&a_ha_04_06,  ha[k][l]);
		accum (&a_VA_04_06,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_04_06,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 0.6) && (a[k][l] < 0.8) )
	{
		accum (&a_ng_06_08,  1.0);
		accum (&a_s_06_08,  s[k][l]);
		accum (&a_q_06_08,  q[k][l]);
		accum (&a_ha_06_08,  ha[k][l]);
		accum (&a_VA_06_08,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_06_08,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 0.8) && (a[k][l] < 1.0) )
	{
		accum (&a_ng_08_10,  1.0);
		accum (&a_s_08_10,  s[k][l]);
		accum (&a_q_08_10,  q[k][l]);
		accum (&a_ha_08_10,  ha[k][l]);
		accum (&a_VA_08_10,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_08_10,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 1.0) && (a[k][l] < 1.2) )
	{
		accum (&a_ng_10_12,  1.0);
		accum (&a_s_10_12,  s[k][l]);
		accum (&a_q_10_12,  q[k][l]);
		accum (&a_ha_10_12,  ha[k][l]);
		accum (&a_VA_10_12,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_10_12,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 1.2) && (a[k][l] < 1.4) )
	{
		accum (&a_ng_12_14,  1.0);
		accum (&a_s_12_14,  s[k][l]);
		accum (&a_q_12_14,  q[k][l]);
		accum (&a_ha_12_14,  ha[k][l]);
		accum (&a_VA_12_14,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_12_14,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 1.4) && (a[k][l] < 1.6) )
	{
		accum (&a_ng_14_16,  1.0);
		accum (&a_s_14_16,  s[k][l]);
		accum (&a_q_14_16,  q[k][l]);
		accum (&a_ha_14_16,  ha[k][l]);
		accum (&a_VA_14_16,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_14_16,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 1.6) && (a[k][l] < 1.8) )
	{
		accum (&a_ng_16_18,  1.0);
		accum (&a_s_16_18,  s[k][l]);
		accum (&a_q_16_18,  q[k][l]);
		accum (&a_ha_16_18,  ha[k][l]);
		accum (&a_VA_16_18,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_16_18,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (a[k][l] >= 1.8) && (a[k][l] < 2.0) )
	{
		accum (&a_ng_18_20,  1.0);
		accum (&a_s_18_20,  s[k][l]);
		accum (&a_q_18_20,  q[k][l]);
		accum (&a_ha_18_20,  ha[k][l]);
		accum (&a_VA_18_20,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_18_20,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if (a[k][l] >= 2.0)
	{
		accum (&a_ng_20,  1.0);
		accum (&a_s_20,  s[k][l]);
		accum (&a_q_20,  q[k][l]);
		accum (&a_ha_20,  ha[k][l]);
		accum (&a_VA_20,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&a_VD_20,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
}


/* ***************************************************************** */


distribution_gens_on_ha()
{
	if ( (ha[k][l] >= 0.0) && (ha[k][l] < 0.2) )
	{
		accum (&ha_ng_00_02,  1.0);
		accum (&ha_q_00_02,  q[k][l]);
		if (s[k][l] != 0.0)   accum (&ha_hs_00_02,  hs[k][l]);
		accum (&ha_VA_00_02,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&ha_VD_00_02,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (ha[k][l] >= 0.2) && (ha[k][l] < 0.4) )
	{
		accum (&ha_ng_02_04,  1.0);
		accum (&ha_q_02_04,  q[k][l]);
		if (s[k][l] != 0.0)   accum (&ha_hs_02_04,  hs[k][l]);
		accum (&ha_VA_02_04,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&ha_VD_02_04,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (ha[k][l] >= 0.4) && (ha[k][l] < 0.6) )
	{
		accum (&ha_ng_04_06,  1.0);
		accum (&ha_q_04_06,  q[k][l]);
		if (s[k][l] != 0.0)   accum (&ha_hs_04_06,  hs[k][l]);
		accum (&ha_VA_04_06,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&ha_VD_04_06,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (ha[k][l] >= 0.6) && (ha[k][l] < 0.8) )
	{
		accum (&ha_ng_06_08,  1.0);
		accum (&ha_q_06_08,  q[k][l]);
		if (s[k][l] != 0.0)   accum (&ha_hs_06_08,  hs[k][l]);
		accum (&ha_VA_06_08,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&ha_VD_06_08,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (ha[k][l] >= 0.8) && (ha[k][l] <= 1.0) )
	{
		accum (&ha_ng_08_10,  1.0);
		accum (&ha_q_08_10,  q[k][l]);
		if (s[k][l] != 0.0)   accum (&ha_hs_08_10,  hs[k][l]);
		accum (&ha_VA_08_10,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&ha_VD_08_10,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
}


/* ***************************************************************** */


distribution_gens_on_s()
{
	if (s[k][l] == 0.0)
	{
		accum (&s_ng_0,  1.0);
		accum (&s_q_0,  q[k][l]);
		accum (&s_hs_0,  hs[k][l]);
		accum (&s_VA_0,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_0,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )
	{
		accum (&s_ng_0_106,  1.0);
		accum (&s_q_0_106,  q[k][l]);
		accum (&s_hs_0_106,  hs[k][l]);
		accum (&s_VA_0_106,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_0_106,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )
	{
		accum (&s_ng_106_104,  1.0);
		accum (&s_q_106_104,  q[k][l]);
		accum (&s_hs_106_104,  hs[k][l]);
		accum (&s_VA_106_104,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_106_104,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )
	{
		accum (&s_ng_104_102,  1.0);
		accum (&s_q_104_102,  q[k][l]);
		accum (&s_hs_104_102,  hs[k][l]);
		accum (&s_VA_104_102,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_104_102,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )
	{
		accum (&s_ng_102_101,  1.0);
		accum (&s_q_102_101,  q[k][l]);
		accum (&s_hs_102_101,  hs[k][l]);
		accum (&s_VA_102_101,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_102_101,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )
	{
		accum (&s_ng_01_02,  1.0);
		accum (&s_q_01_02,  q[k][l]);
		accum (&s_hs_01_02,  hs[k][l]);
		accum (&s_VA_01_02,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_01_02,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )
	{
		accum (&s_ng_02_04,  1.0);
		accum (&s_q_02_04,  q[k][l]);
		accum (&s_hs_02_04,  hs[k][l]);
		accum (&s_VA_02_04,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_02_04,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )
	{
		accum (&s_ng_04_06,  1.0);
		accum (&s_q_04_06,  q[k][l]);
		accum (&s_hs_04_06,  hs[k][l]);
		accum (&s_VA_04_06,  2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_04_06,  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )
	{
		accum (&s_ng_06_10,  1.0);
		accum (&s_q_06_10,  q[k][l]);
		accum (&s_hs_06_10,  hs[k][l]);
		accum (&s_VA_06_10,
				2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_06_10,
				pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}
	else if (s[k][l] == (-1.0))
	{
		accum (&s_ng_10,  1.0);
		accum (&s_q_10,  q[k][l]);
		accum (&s_hs_10,  hs[k][l]);
		accum (&s_VA_10,
				2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
		accum (&s_VD_10,
				pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0));
	}

}


/* **************************************************************** */


distribution_blocks_on_q()
{
		accum (&Ave_ng_all,  accsum(&ng_all));
		accum (&Ave_VA_all,  accsum(&VA_all));
		accum (&Ave_VD_all,  accsum(&VD_all));

if (tracelevel!=0)    fprintf(fptr,"escribe VA_all=%f  Ave_VA_all=%f \n", accsum(&VA_all), accsum(&Ave_VA_all));

		if (accsum(&ng_all) != 0.0)
		{
			accum (&Ave_s_all,  accmean(&s_all));
			accum (&Ave_a_all,  accmean(&a_all));
			if (accmean(&s_all) != 0.0)   accum (&Ave_hs_all,  accmean(&hs_all));
			accum (&Ave_ha_all,  accmean(&ha_all));
		}

		accum (&Ave_q_ng_00_01,  accsum(&q_ng_00_01));
		accum (&Ave_q_VA_00_01,  accsum(&q_VA_00_01));
		accum (&Ave_q_VD_00_01,  accsum(&q_VD_00_01));
		if (accsum(&q_ng_00_01) != 0.0)
		{
			accum (&Ave_q_s_00_01,  accmean(&q_s_00_01));
			accum (&Ave_q_a_00_01,  accmean(&q_a_00_01));
			if (accmean(&q_s_00_01) != 0.0)
				accum (&Ave_q_hs_00_01, accmean(&q_hs_00_01));
			accum (&Ave_q_ha_00_01,  accmean(&q_ha_00_01));
		}

		accum (&Ave_q_ng_01_02,  accsum(&q_ng_01_02));
		accum (&Ave_q_VA_01_02,  accsum(&q_VA_01_02));
		accum (&Ave_q_VD_01_02,  accsum(&q_VD_01_02));
		if (accsum(&q_ng_01_02) != 0.0)
		{
			accum (&Ave_q_s_01_02,  accmean(&q_s_01_02));
			accum (&Ave_q_a_01_02,  accmean(&q_a_01_02));
			if (accmean(&q_s_01_02) != 0.0)
				accum (&Ave_q_hs_01_02, accmean(&q_hs_01_02));
			accum (&Ave_q_ha_01_02,  accmean(&q_ha_01_02));
		}

		accum (&Ave_q_ng_02_03,  accsum(&q_ng_02_03));
		accum (&Ave_q_VA_02_03,  accsum(&q_VA_02_03));
		accum (&Ave_q_VD_02_03,  accsum(&q_VD_02_03));
		if (accsum(&q_ng_02_03) != 0.0)
		{
			accum (&Ave_q_s_02_03,  accmean(&q_s_02_03));
			accum (&Ave_q_a_02_03,  accmean(&q_a_02_03));
			if (accmean(&q_s_02_03) != 0.0)
				accum (&Ave_q_hs_02_03, accmean(&q_hs_02_03));
			accum (&Ave_q_ha_02_03,  accmean(&q_ha_02_03));
		}

		accum (&Ave_q_ng_03_04,  accsum(&q_ng_03_04));
		accum (&Ave_q_VA_03_04,  accsum(&q_VA_03_04));
		accum (&Ave_q_VD_03_04,  accsum(&q_VD_03_04));
		if (accsum(&q_ng_03_04) != 0.0)
		{
			accum (&Ave_q_s_03_04,  accmean(&q_s_03_04));
			accum (&Ave_q_a_03_04,  accmean(&q_a_03_04));
			if (accmean(&q_s_03_04) != 0.0)
				accum (&Ave_q_hs_03_04, accmean(&q_hs_03_04));
			accum (&Ave_q_ha_03_04,  accmean(&q_ha_03_04));
		}

		accum (&Ave_q_ng_04_05,  accsum(&q_ng_04_05));
		accum (&Ave_q_VA_04_05,  accsum(&q_VA_04_05));
		accum (&Ave_q_VD_04_05,  accsum(&q_VD_04_05));
		if (accsum(&q_ng_04_05) != 0.0)
		{
			accum (&Ave_q_s_04_05,  accmean(&q_s_04_05));
			accum (&Ave_q_a_04_05,  accmean(&q_a_04_05));
			if (accmean(&q_s_04_05) != 0.0)
				accum (&Ave_q_hs_04_05, accmean(&q_hs_04_05));
			accum (&Ave_q_ha_04_05,  accmean(&q_ha_04_05)); 
		}

		accum (&Ave_q_ng_05_06,  accsum(&q_ng_05_06));
		accum (&Ave_q_VA_05_06,  accsum(&q_VA_05_06));
		accum (&Ave_q_VD_05_06,  accsum(&q_VD_05_06));
		if (accsum(&q_ng_05_06) != 0.0)
		{
			accum (&Ave_q_s_05_06,  accmean(&q_s_05_06));
			accum (&Ave_q_a_05_06,  accmean(&q_a_05_06));
			if (accmean(&q_s_05_06) != 0.0)
				accum (&Ave_q_hs_05_06, accmean(&q_hs_05_06));
			accum (&Ave_q_ha_05_06,  accmean(&q_ha_05_06));
		}

		accum (&Ave_q_ng_06_07,  accsum(&q_ng_06_07));
		accum (&Ave_q_VA_06_07,  accsum(&q_VA_06_07));
		accum (&Ave_q_VD_06_07,  accsum(&q_VD_06_07));
		if (accsum(&q_ng_06_07) != 0.0)
		{
			accum (&Ave_q_s_06_07,  accmean(&q_s_06_07));
			accum (&Ave_q_a_06_07,  accmean(&q_a_06_07));
			if (accmean(&q_s_06_07) != 0.0)
				accum (&Ave_q_hs_06_07, accmean(&q_hs_06_07));
			accum (&Ave_q_ha_06_07,  accmean(&q_ha_06_07));
		}

		accum (&Ave_q_ng_07_08,  accsum(&q_ng_07_08));
		accum (&Ave_q_VA_07_08,  accsum(&q_VA_07_08));
		accum (&Ave_q_VD_07_08,  accsum(&q_VD_07_08));
		if (accsum(&q_ng_07_08) != 0.0)
		{
			accum (&Ave_q_s_07_08,  accmean(&q_s_07_08));
			accum (&Ave_q_a_07_08,  accmean(&q_a_07_08));
			if (accmean(&q_s_07_08) != 0.0)
				accum (&Ave_q_hs_07_08, accmean(&q_hs_07_08));
			accum (&Ave_q_ha_07_08,  accmean(&q_ha_07_08));
		}

		accum (&Ave_q_ng_08_09,  accsum(&q_ng_08_09));
		accum (&Ave_q_VA_08_09,  accsum(&q_VA_08_09));
		accum (&Ave_q_VD_08_09,  accsum(&q_VD_08_09));
		if (accsum(&q_ng_08_09) != 0.0)
		{
			accum (&Ave_q_s_08_09,  accmean(&q_s_08_09));
			accum (&Ave_q_a_08_09,  accmean(&q_a_08_09));
			if (accmean(&q_s_08_09) != 0.0)
				accum (&Ave_q_hs_08_09, accmean(&q_hs_08_09));
			accum (&Ave_q_ha_08_09,  accmean(&q_ha_08_09));
		}

		accum (&Ave_q_ng_09_10,  accsum(&q_ng_09_10));
		accum (&Ave_q_VA_09_10,  accsum(&q_VA_09_10));
		accum (&Ave_q_VD_09_10,  accsum(&q_VD_09_10));
		if (accsum(&q_ng_09_10) != 0.0)
		{
			accum (&Ave_q_s_09_10,  accmean(&q_s_09_10));
			accum (&Ave_q_a_09_10,  accmean(&q_a_09_10));
			if (accmean(&q_s_09_10) != 0.0)
				accum (&Ave_q_hs_09_10, accmean(&q_hs_09_10));
			accum (&Ave_q_ha_09_10,  accmean(&q_ha_09_10));
		}
}


/* **************************************************************** */


distribution_blocks_on_a()
{
	/* ***** negative values ***** */

		accum (&Ave_a_ng_m00_m02,  accsum(&a_ng_m00_m02));
		if (accsum(&a_ng_m00_m02) != 0.0)
		{
		accum (&Ave_a_s_m00_m02,  accmean(&a_s_m00_m02));
		accum (&Ave_a_q_m00_m02,  accmean(&a_q_m00_m02));
		accum (&Ave_a_ha_m00_m02,  accmean(&a_ha_m00_m02));
		}
		accum (&Ave_a_VA_m00_m02,  accsum(&a_VA_m00_m02));
		accum (&Ave_a_VD_m00_m02,  accsum(&a_VD_m00_m02));

		accum (&Ave_a_ng_m02_m04,  accsum(&a_ng_m02_m04));
		if (accsum(&a_ng_m02_m04) != 0.0)
		{
		accum (&Ave_a_s_m02_m04,  accmean(&a_s_m02_m04));
		accum (&Ave_a_q_m02_m04,  accmean(&a_q_m02_m04));
		accum (&Ave_a_ha_m02_m04,  accmean(&a_ha_m02_m04));
		}
		accum (&Ave_a_VA_m02_m04,  accsum(&a_VA_m02_m04));
		accum (&Ave_a_VD_m02_m04,  accsum(&a_VD_m02_m04));

		accum (&Ave_a_ng_m04_m06,  accsum(&a_ng_m04_m06));
		if (accsum(&a_ng_m04_m06) != 0.0)
		{
		accum (&Ave_a_s_m04_m06,  accmean(&a_s_m04_m06));
		accum (&Ave_a_q_m04_m06,  accmean(&a_q_m04_m06));
		accum (&Ave_a_ha_m04_m06,  accmean(&a_ha_m04_m06));
		}
		accum (&Ave_a_VA_m04_m06,  accsum(&a_VA_m04_m06));
		accum (&Ave_a_VD_m04_m06,  accsum(&a_VD_m04_m06));

		accum (&Ave_a_ng_m06_m08,  accsum(&a_ng_m06_m08));
		if (accsum(&a_ng_m06_m08) != 0.0)
		{
		accum (&Ave_a_s_m06_m08,  accmean(&a_s_m06_m08));
		accum (&Ave_a_q_m06_m08,  accmean(&a_q_m06_m08));
		accum (&Ave_a_ha_m06_m08,  accmean(&a_ha_m06_m08));
		}
		accum (&Ave_a_VA_m06_m08,  accsum(&a_VA_m06_m08));
		accum (&Ave_a_VD_m06_m08,  accsum(&a_VD_m06_m08));

		accum (&Ave_a_ng_m08_m10,  accsum(&a_ng_m08_m10));
		if (accsum(&a_ng_m08_m10) != 0.0)
		{
		accum (&Ave_a_s_m08_m10,  accmean(&a_s_m08_m10));
		accum (&Ave_a_q_m08_m10,  accmean(&a_q_m08_m10));
		accum (&Ave_a_ha_m08_m10,  accmean(&a_ha_m08_m10));
		}
		accum (&Ave_a_VA_m08_m10,  accsum(&a_VA_m08_m10));
		accum (&Ave_a_VD_m08_m10,  accsum(&a_VD_m08_m10));

		accum (&Ave_a_ng_m10_m12,  accsum(&a_ng_m10_m12));
		if (accsum(&a_ng_m10_m12) != 0.0)
		{
		accum (&Ave_a_s_m10_m12,  accmean(&a_s_m10_m12));
		accum (&Ave_a_q_m10_m12,  accmean(&a_q_m10_m12));
		accum (&Ave_a_ha_m10_m12,  accmean(&a_ha_m10_m12));
		}
		accum (&Ave_a_VA_m10_m12,  accsum(&a_VA_m10_m12));
		accum (&Ave_a_VD_m10_m12,  accsum(&a_VD_m10_m12));

		accum (&Ave_a_ng_m12_m14,  accsum(&a_ng_m12_m14));
		if (accsum(&a_ng_m12_m14) != 0.0)
		{
		accum (&Ave_a_s_m12_m14,  accmean(&a_s_m12_m14));
		accum (&Ave_a_q_m12_m14,  accmean(&a_q_m12_m14));
		accum (&Ave_a_ha_m12_m14,  accmean(&a_ha_m12_m14));
		}
		accum (&Ave_a_VA_m12_m14,  accsum(&a_VA_m12_m14));
		accum (&Ave_a_VD_m12_m14,  accsum(&a_VD_m12_m14));

		accum (&Ave_a_ng_m14_m16,  accsum(&a_ng_m14_m16));
		if (accsum(&a_ng_m14_m16) != 0.0)
		{
		accum (&Ave_a_s_m14_m16,  accmean(&a_s_m14_m16));
		accum (&Ave_a_q_m14_m16,  accmean(&a_q_m14_m16));
		accum (&Ave_a_ha_m14_m16,  accmean(&a_ha_m14_m16));
		}
		accum (&Ave_a_VA_m14_m16,  accsum(&a_VA_m14_m16));
		accum (&Ave_a_VD_m14_m16,  accsum(&a_VD_m14_m16));

		accum (&Ave_a_ng_m16_m18,  accsum(&a_ng_m16_m18));
		if (accsum(&a_ng_m16_m18) != 0.0)
		{
		accum (&Ave_a_s_m16_m18,  accmean(&a_s_m16_m18));
		accum (&Ave_a_q_m16_m18,  accmean(&a_q_m16_m18));
		accum (&Ave_a_ha_m16_m18,  accmean(&a_ha_m16_m18));
		}
		accum (&Ave_a_VA_m16_m18,  accsum(&a_VA_m16_m18));
		accum (&Ave_a_VD_m16_m18,  accsum(&a_VD_m16_m18));

		accum (&Ave_a_ng_m18_m20,  accsum(&a_ng_m18_m20));
		if (accsum(&a_ng_m18_m20) != 0.0)
		{
		accum (&Ave_a_s_m18_m20,  accmean(&a_s_m18_m20));
		accum (&Ave_a_q_m18_m20,  accmean(&a_q_m18_m20));
		accum (&Ave_a_ha_m18_m20,  accmean(&a_ha_m18_m20));
		}
		accum (&Ave_a_VA_m18_m20,  accsum(&a_VA_m18_m20));
		accum (&Ave_a_VD_m18_m20,  accsum(&a_VD_m18_m20));

		accum (&Ave_a_ng_m20,  accsum(&a_ng_m20));
		if (accsum(&a_ng_m20) != 0.0)
		{
		accum (&Ave_a_s_m20,  accmean(&a_s_m20));
		accum (&Ave_a_q_m20,  accmean(&a_q_m20));
		accum (&Ave_a_ha_m20,  accmean(&a_ha_m20));
		}
		accum (&Ave_a_VA_m20,  accsum(&a_VA_m20));
		accum (&Ave_a_VD_m20,  accsum(&a_VD_m20));

		/* ***** positive values ***** */

		accum (&Ave_a_ng_00_02,  accsum(&a_ng_00_02));
		if (accsum(&a_ng_00_02) != 0.0)
		{
		accum (&Ave_a_s_00_02,  accmean(&a_s_00_02));
		accum (&Ave_a_q_00_02,  accmean(&a_q_00_02));
		accum (&Ave_a_ha_00_02,  accmean(&a_ha_00_02));
		}
		accum (&Ave_a_VA_00_02,  accsum(&a_VA_00_02));
		accum (&Ave_a_VD_00_02,  accsum(&a_VD_00_02));

		accum (&Ave_a_ng_02_04,  accsum(&a_ng_02_04));
		if (accsum(&a_ng_02_04) != 0.0)
		{
		accum (&Ave_a_s_02_04,  accmean(&a_s_02_04));
		accum (&Ave_a_q_02_04,  accmean(&a_q_02_04));
		accum (&Ave_a_ha_02_04,  accmean(&a_ha_02_04));
		}
		accum (&Ave_a_VA_02_04,  accsum(&a_VA_02_04));
		accum (&Ave_a_VD_02_04,  accsum(&a_VD_02_04));

		accum (&Ave_a_ng_04_06,  accsum(&a_ng_04_06));
		if (accsum(&a_ng_04_06) != 0.0)
		{
		accum (&Ave_a_s_04_06,  accmean(&a_s_04_06));
		accum (&Ave_a_q_04_06,  accmean(&a_q_04_06));
		accum (&Ave_a_ha_04_06,  accmean(&a_ha_04_06));
		}
		accum (&Ave_a_VA_04_06,  accsum(&a_VA_04_06));
		accum (&Ave_a_VD_04_06,  accsum(&a_VD_04_06));

		accum (&Ave_a_ng_06_08,  accsum(&a_ng_06_08));
		if (accsum(&a_ng_06_08) != 0.0)
		{
		accum (&Ave_a_s_06_08,  accmean(&a_s_06_08));
		accum (&Ave_a_q_06_08,  accmean(&a_q_06_08));
		accum (&Ave_a_ha_06_08,  accmean(&a_ha_06_08));
		}
		accum (&Ave_a_VA_06_08,  accsum(&a_VA_06_08));
		accum (&Ave_a_VD_06_08,  accsum(&a_VD_06_08));

		accum (&Ave_a_ng_08_10,  accsum(&a_ng_08_10));
		if (accsum(&a_ng_08_10) != 0.0)
		{
		accum (&Ave_a_s_08_10,  accmean(&a_s_08_10));
		accum (&Ave_a_q_08_10,  accmean(&a_q_08_10));
		accum (&Ave_a_ha_08_10,  accmean(&a_ha_08_10));
		}
		accum (&Ave_a_VA_08_10,  accsum(&a_VA_08_10));
		accum (&Ave_a_VD_08_10,  accsum(&a_VD_08_10));

		accum (&Ave_a_ng_10_12,  accsum(&a_ng_10_12));
		if (accsum(&a_ng_10_12) != 0.0)
		{
		accum (&Ave_a_s_10_12,  accmean(&a_s_10_12));
		accum (&Ave_a_q_10_12,  accmean(&a_q_10_12));
		accum (&Ave_a_ha_10_12,  accmean(&a_ha_10_12));
		}
		accum (&Ave_a_VA_10_12,  accsum(&a_VA_10_12));
		accum (&Ave_a_VD_10_12,  accsum(&a_VD_10_12));

		accum (&Ave_a_ng_12_14,  accsum(&a_ng_12_14));
		if (accsum(&a_ng_12_14) != 0.0)
		{
		accum (&Ave_a_s_12_14,  accmean(&a_s_12_14));
		accum (&Ave_a_q_12_14,  accmean(&a_q_12_14));
		accum (&Ave_a_ha_12_14,  accmean(&a_ha_12_14));
		}
		accum (&Ave_a_VA_12_14,  accsum(&a_VA_12_14));
		accum (&Ave_a_VD_12_14,  accsum(&a_VD_12_14));

		accum (&Ave_a_ng_14_16,  accsum(&a_ng_14_16));
		if (accsum(&a_ng_14_16) != 0.0)
		{
		accum (&Ave_a_s_14_16,  accmean(&a_s_14_16));
		accum (&Ave_a_q_14_16,  accmean(&a_q_14_16));
		accum (&Ave_a_ha_14_16,  accmean(&a_ha_14_16));
		}
		accum (&Ave_a_VA_14_16,  accsum(&a_VA_14_16));
		accum (&Ave_a_VD_14_16,  accsum(&a_VD_14_16));

		accum (&Ave_a_ng_16_18,  accsum(&a_ng_16_18));
		if (accsum(&a_ng_16_18) != 0.0)
		{
		accum (&Ave_a_s_16_18,  accmean(&a_s_16_18));
		accum (&Ave_a_q_16_18,  accmean(&a_q_16_18));
		accum (&Ave_a_ha_16_18,  accmean(&a_ha_16_18));
		}
		accum (&Ave_a_VA_16_18,  accsum(&a_VA_16_18));
		accum (&Ave_a_VD_16_18,  accsum(&a_VD_16_18));

		accum (&Ave_a_ng_18_20,  accsum(&a_ng_18_20));
		if (accsum(&a_ng_18_20) != 0.0)
		{
		accum (&Ave_a_s_18_20,  accmean(&a_s_18_20));
		accum (&Ave_a_q_18_20,  accmean(&a_q_18_20));
		accum (&Ave_a_ha_18_20,  accmean(&a_ha_18_20));
		}
		accum (&Ave_a_VA_18_20,  accsum(&a_VA_18_20));
		accum (&Ave_a_VD_18_20,  accsum(&a_VD_18_20));

		accum (&Ave_a_ng_20,  accsum(&a_ng_20));
		if (accsum(&a_ng_20) != 0.0)
		{
		accum (&Ave_a_s_20,  accmean(&a_s_20));
		accum (&Ave_a_q_20,  accmean(&a_q_20));
		accum (&Ave_a_ha_20,  accmean(&a_ha_20));
		}
		accum (&Ave_a_VA_20,  accsum(&a_VA_20));
		accum (&Ave_a_VD_20,  accsum(&a_VD_20));
}


/* ***************************************************************** */


distribution_blocks_on_ha()
{
	accum (&Ave_ha_ng_00_02,  accsum(&ha_ng_00_02));
	if (accsum(&ha_ng_00_02) != 0.0)
	{
	accum (&Ave_ha_q_00_02,  accmean(&ha_q_00_02));
		accum (&Ave_ha_hs_00_02,  accmean(&ha_hs_00_02));
	}
accum (&Ave_ha_VA_00_02,  accsum(&ha_VA_00_02));
	accum (&Ave_ha_VD_00_02,  accsum(&ha_VD_00_02));

	accum (&Ave_ha_ng_02_04,  accsum(&ha_ng_02_04));
	if (accsum(&ha_ng_02_04) != 0.0)
	{
		accum (&Ave_ha_q_02_04,  accmean(&ha_q_02_04));
		accum (&Ave_ha_hs_02_04,  accmean(&ha_hs_02_04));
	}
accum (&Ave_ha_VA_02_04,  accsum(&ha_VA_02_04));
	accum (&Ave_ha_VD_02_04,  accsum(&ha_VD_02_04));

	accum (&Ave_ha_ng_04_06,  accsum(&ha_ng_04_06));
	if (accsum(&ha_ng_04_06) != 0.0)
	{
		accum (&Ave_ha_q_04_06,  accmean(&ha_q_04_06));
	accum (&Ave_ha_hs_04_06,  accmean(&ha_hs_04_06));
	}
accum (&Ave_ha_VA_04_06,  accsum(&ha_VA_04_06));
	accum (&Ave_ha_VD_04_06,  accsum(&ha_VD_04_06));

	accum (&Ave_ha_ng_06_08,  accsum(&ha_ng_06_08));
	if (accsum(&ha_ng_06_08) != 0.0)
	{
	accum (&Ave_ha_q_06_08,  accmean(&ha_q_06_08));
	accum (&Ave_ha_hs_06_08,  accmean(&ha_hs_06_08));
	}
accum (&Ave_ha_VA_06_08,  accsum(&ha_VA_06_08));
	accum (&Ave_ha_VD_06_08,  accsum(&ha_VD_06_08));

	accum (&Ave_ha_ng_08_10,  accsum(&ha_ng_08_10));
	if (accsum(&ha_ng_08_10) != 0.0)
	{
		accum (&Ave_ha_q_08_10,  accmean(&ha_q_08_10));
	accum (&Ave_ha_hs_08_10,  accmean(&ha_hs_08_10));
	}
	accum (&Ave_ha_VA_08_10,  accsum(&ha_VA_08_10));
	accum (&Ave_ha_VD_08_10,  accsum(&ha_VD_08_10));
}


/* ***************************************************************** */


distribution_blocks_on_s()
{
	accum (&Ave_s_ng_0,  accsum(&s_ng_0));
	if (accsum(&s_ng_0) != 0.0)
	{
		accum (&Ave_s_q_0,  accmean(&s_q_0));
		accum (&Ave_s_hs_0,  accmean(&s_hs_0));
	}
	accum (&Ave_s_VA_0,  accsum(&s_VA_0));
	accum (&Ave_s_VD_0,  accsum(&s_VD_0));

	accum (&Ave_s_ng_0_106,  accsum(&s_ng_0_106));
	if (accsum(&s_ng_0_106) != 0.0)
	{
		accum (&Ave_s_q_0_106,  accmean(&s_q_0_106));
		accum (&Ave_s_hs_0_106,  accmean(&s_hs_0_106));
	}
	accum (&Ave_s_VA_0_106,  accsum(&s_VA_0_106));
	accum (&Ave_s_VD_0_106,  accsum(&s_VD_0_106));

	accum (&Ave_s_ng_106_104,  accsum(&s_ng_106_104));
	if (accsum(&s_ng_106_104) != 0.0)
	{
		accum (&Ave_s_q_106_104,  accmean(&s_q_106_104));
		accum (&Ave_s_hs_106_104,  accmean(&s_hs_106_104));
	}
	accum (&Ave_s_VA_106_104,  accsum(&s_VA_106_104));
	accum (&Ave_s_VD_106_104,  accsum(&s_VD_106_104));

	accum (&Ave_s_ng_104_102,  accsum(&s_ng_104_102));
	if (accsum(&s_ng_104_102) != 0.0)
	{
		accum (&Ave_s_q_104_102,  accmean(&s_q_104_102));
		accum (&Ave_s_hs_104_102,  accmean(&s_hs_104_102));
	}
	accum (&Ave_s_VA_104_102,  accsum(&s_VA_104_102));
	accum (&Ave_s_VD_104_102,  accsum(&s_VD_104_102));

	accum (&Ave_s_ng_102_101,  accsum(&s_ng_102_101));
	if (accsum(&s_ng_102_101) != 0.0)
	{
		accum (&Ave_s_q_102_101,  accmean(&s_q_102_101));
		accum (&Ave_s_hs_102_101,  accmean(&s_hs_102_101));
	}
	accum (&Ave_s_VA_102_101,  accsum(&s_VA_102_101));
	accum (&Ave_s_VD_102_101,  accsum(&s_VD_102_101));

	accum (&Ave_s_ng_01_02,  accsum(&s_ng_01_02));
	if (accsum(&s_ng_01_02) != 0.0)
	{
		accum (&Ave_s_q_01_02,  accmean(&s_q_01_02));
		accum (&Ave_s_hs_01_02,  accmean(&s_hs_01_02));
	}
	accum (&Ave_s_VA_01_02,  accsum(&s_VA_01_02));
	accum (&Ave_s_VD_01_02,  accsum(&s_VD_01_02));

	accum (&Ave_s_ng_02_04,  accsum(&s_ng_02_04));
	if (accsum(&s_ng_02_04) != 0.0)
	{
		accum (&Ave_s_q_02_04,  accmean(&s_q_02_04));
		accum (&Ave_s_hs_02_04,  accmean(&s_hs_02_04));
	}
	accum (&Ave_s_VA_02_04,  accsum(&s_VA_02_04));
	accum (&Ave_s_VD_02_04,  accsum(&s_VD_02_04));

	accum (&Ave_s_ng_04_06,  accsum(&s_ng_04_06));
	if (accsum(&s_ng_04_06) != 0.0)
	{
		accum (&Ave_s_q_04_06,  accmean(&s_q_04_06));
		accum (&Ave_s_hs_04_06,  accmean(&s_hs_04_06));
	}
	accum (&Ave_s_VA_04_06,  accsum(&s_VA_04_06));
	accum (&Ave_s_VD_04_06,  accsum(&s_VD_04_06));

	accum (&Ave_s_ng_06_10,  accsum(&s_ng_06_10));
	if (accsum(&s_ng_06_10) != 0.0)
	{
		accum (&Ave_s_q_06_10,  accmean(&s_q_06_10));
		accum (&Ave_s_hs_06_10,  accmean(&s_hs_06_10));
	}
	accum (&Ave_s_VA_06_10,  accsum(&s_VA_06_10));
	accum (&Ave_s_VD_06_10,  accsum(&s_VD_06_10));

	accum (&Ave_s_ng_10,  accsum(&s_ng_10));
	if (accsum(&s_ng_10) != 0.0)
	{
		accum (&Ave_s_q_10,  accmean(&s_q_10));
		accum (&Ave_s_hs_10,  accmean(&s_hs_10));
	}
	accum (&Ave_s_VA_10,  accsum(&s_VA_10));
	accum (&Ave_s_VD_10,  accsum(&s_VD_10));
}


/* **************************************************************** */


void mutation_neutral (gm,q,initialgen,cumH)
int gm[][MMM][2], initialgen[][31];
double q[][31], cumH[][31];
{
    /* NEUTRAL GENES: (POISSON) 2N(Lambda_a) NEW MUTATIONS (TWO DIRECTIONS) */

    if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants\n");

    muts = mutationnumber();

    for (m=0; m<muts; m++)
    {
	ran_i = (int)(uniform()*NIND);
	ran_k = (int)(uniform()*NCRO);
	ran_h = (int)(uniform()*2.0);
	if ( (gm[ran_i][ran_k][ran_h] & RM[0])==RM[0] )
	    gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] & (~RM[0]));
	else
	    gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[0]);
    }
}


/* **************************************************************** */


void mutation_selected (gm,q,s,hs,a,ha,initialgen,cumH)
int gm[][MMM][2], initialgen[][31];
double q[][31], s[][31], hs[][31], a[][31], ha[][31], cumH[][31];
{
    /* SELECTED GENES: (POISSON) 2N(Lambda_a + Lambda_L) NEW MUTATIONS */
    
    muts = mutationnumber();
    mutsL = mutationnumberL();

    if ( (muts + mutsL) == 0 ) goto label;

    if (tracelevel!=0)    fprintf(fptr,"\n New selected mutants\n");

    countNoSS = 0;

    for (k=0; k<NCRO; k++)
    {
	for (l=1; l<NLOCI; l++)
	{
	    if (q[k][l]==0.0)
	    {
		countNoSS += 1;
		NoSS_k[countNoSS-1] = k;
		NoSS_l[countNoSS-1] = l;
	    }
	}
    }
    
    mutants_ocurred = 0;

    if (countNoSS != 0)
    {
	disorder_NoSS (NoSS_k,NoSS_l);

	for (m=0; m<countNoSS; m++)
	{
	    if (mutants_ocurred==(muts+mutsL))    goto label;

	    ran_i = (int)(uniform()*NIND);
   	    ran_h = (int)(uniform()*2.0);
    	    gm[ran_i][NoSS_k[m]][ran_h]=
				(gm[ran_i][NoSS_k[m]][ran_h] | RM[NoSS_l[m]]);
	    mutants_ocurred += 1;
	    initialgen[NoSS_k[m]][NoSS_l[m]] = gen+1;
	    cumH[NoSS_k[m]][NoSS_l[m]] = 0.0;

	 /* ****** Lethal mutations with effect on the trait ****** */

	 if(mutants_ocurred <= mutsL)
	 {
		zs = 1.0;
		if(rho==99.0)	za = zs;
		else
		{
			za = normal(ave_aL, ave_aL/2.0);
			if (uniform() < 0.5) za = (-za);
		}
		s[NoSS_k[m]][NoSS_l[m]] = (-zs);
		meanmutsL_s[(gen+1)/block] +=  1.0;
		hs[NoSS_k[m]][NoSS_l[m]] = 0.02;
	 }
	 else
	 {
	    /* ****** values of s, hs, a and ha for the new mutation ****** */

	   if(whatproportion==0)
	   {
	    	if(uniform() < Psi)
	   	{
		       if(rho==0.0)
		       {
		    		za = gengam (alpha_a, beta_a);
		    		zs = gengam (alpha_s, beta_s);
				if (zs > 1.0)		zs = 1.0;	
		       }
		       else if(rho==99.0)
		       {
		    		za = gengam (alpha_a, beta_a);
		    		zs = za;
				if (zs > 1.0)		zs = 1.0;
		       }
		       else
		       {
				bvgam(alpha_s, alpha_a, beta_s, beta_a, rho, &zs, &za);
		       }
	    	       if (zs > 1.0)		zs = 1.0;
			meanmuts_s[(gen+1)/block] +=  1.0;
	   	       if (uniform()<P_s)	s[NoSS_k[m]][NoSS_l[m]] = zs;
	    	       else				s[NoSS_k[m]][NoSS_l[m]] = (-zs);

		       if(dom==0)		hs[NoSS_k[m]][NoSS_l[m]] = k_s;
		       else if(dom==1)
		       {
	    	       	if (s[NoSS_k[m]][NoSS_l[m]] < 0.0)
				hs[NoSS_k[m]][NoSS_l[m]] =
				0.5 * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
	    	       	else	hs[NoSS_k[m]][NoSS_l[m]] =
				0.5 * exp(-k_s*s[NoSS_k[m]][NoSS_l[m]]);
		       }
		       else
		       {
	    	       	if (s[NoSS_k[m]][NoSS_l[m]] < 0.0)
				hs[NoSS_k[m]][NoSS_l[m]] =
				uniform() * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
	    	       	else	hs[NoSS_k[m]][NoSS_l[m]] =
				uniform() * exp(-k_s*s[NoSS_k[m]][NoSS_l[m]]);
		       }

	   	       mean_s[(gen+1)/block] += zs;
	  	       mean_s2[(gen+1)/block] += zs*zs;
	 	       mean_s3[(gen+1)/block] +=  zs*zs*zs;
	   	       mean_s4[(gen+1)/block] +=  zs*zs*zs*zs;
	  	       mean_hs[(gen+1)/block] += hs[NoSS_k[m]][NoSS_l[m]];
	   	       mean_a_s[(gen+1)/block] += za;
	  	       mean_a2_s[(gen+1)/block] += za*za;
	  	       mean_sa[(gen+1)/block] += zs*za;
	   	}
	   	else
	   	{
	   		za = gengam (alpha_a, beta_a);
			zs = 0.0;
	   	}
	   }
	   else
	   {
		 if(rho==0)
		{
	    		za = gengam (alpha_a, beta_a);
	    		if(za<threshold)		zs = 0.0;
			else
			{
	    			zs = gengam (alpha_s, beta_s);
				if (zs > 1.0)		zs = 1.0;
	    			meanmuts_s[(gen+1)/block] +=  1.0;
	    			if (uniform()<P_s)	s[NoSS_k[m]][NoSS_l[m]] = zs;
				else			s[NoSS_k[m]][NoSS_l[m]] = (-zs);

		       		if(dom==0)		hs[NoSS_k[m]][NoSS_l[m]] = k_s;
			       	else if(dom==1)
			       	{
	    	       			if (s[NoSS_k[m]][NoSS_l[m]] < 0.0)
						hs[NoSS_k[m]][NoSS_l[m]] = 0.5 * 									exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
	    	       			else	hs[NoSS_k[m]][NoSS_l[m]] = 0.5 *
						exp(-k_s*s[NoSS_k[m]][NoSS_l[m]]);
		       		}
			       	else
			       	{
	    	       			if (s[NoSS_k[m]][NoSS_l[m]] < 0.0)
						hs[NoSS_k[m]][NoSS_l[m]] = uniform() * 								exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
	    	       			else	hs[NoSS_k[m]][NoSS_l[m]] = uniform() *
						exp(-k_s*s[NoSS_k[m]][NoSS_l[m]]);
		       		}

				mean_s[(gen+1)/block] += zs;
				mean_s2[(gen+1)/block] += zs*zs;
				mean_s3[(gen+1)/block] +=  zs*zs*zs;
				mean_s4[(gen+1)/block] +=  zs*zs*zs*zs;
				mean_hs[(gen+1)/block] += hs[NoSS_k[m]][NoSS_l[m]];
			   	mean_a_s[(gen+1)/block] += za;
		  	       	mean_a2_s[(gen+1)/block] += za*za;
				mean_sa[(gen+1)/block] += zs*za;
			}
		}
		else
		{
			bvgam(alpha_s, alpha_a, beta_s, beta_a, rho, &zs, &za);
	    		if(za<threshold)		zs = 0.0;
			else
			{
	    			meanmuts_s[(gen+1)/block] +=  1.0;
	    			if (zs > 1.0)		zs = 1.0;
				if (uniform()<P_s)	s[NoSS_k[m]][NoSS_l[m]] = zs;
	    			else			s[NoSS_k[m]][NoSS_l[m]] = (-zs);

				if(dom==0)		hs[NoSS_k[m]][NoSS_l[m]] = k_s;
				else if(dom==1)
		       	{
	    	       		if (s[NoSS_k[m]][NoSS_l[m]] < 0.0)
					hs[NoSS_k[m]][NoSS_l[m]] = 0.5 * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
	    	       		else	hs[NoSS_k[m]][NoSS_l[m]] = 0.5 * exp(-k_s*s[NoSS_k[m]][NoSS_l[m]]);
		       	}
				else
		       	{
	    	       		if (s[NoSS_k[m]][NoSS_l[m]] < 0.0)
					hs[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
	    	       		else	hs[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(-k_s*s[NoSS_k[m]][NoSS_l[m]]);
		       	}

				mean_s[(gen+1)/block] += zs;
				mean_s2[(gen+1)/block] += zs*zs;
				mean_s3[(gen+1)/block] +=  zs*zs*zs;
				mean_s4[(gen+1)/block] +=  zs*zs*zs*zs;
				mean_hs[(gen+1)/block] += hs[NoSS_k[m]][NoSS_l[m]];
	   			mean_a_s[(gen+1)/block] += za;
			  	mean_a2_s[(gen+1)/block] += za*za;
				mean_sa[(gen+1)/block] += zs*za;
			}
		}
	   }
	 }
	    totalmutant[(gen+1)/block] +=  1;
	    meanmuts_a[(gen+1)/block] +=  1.0;
	    if (uniform()<P_a)	a[NoSS_k[m]][NoSS_l[m]] = za;
	    else			a[NoSS_k[m]][NoSS_l[m]] = (-za);
	    mean_real_a[(gen+1)/block] += a[NoSS_k[m]][NoSS_l[m]];

	    if(dom==0)	ha[NoSS_k[m]][NoSS_l[m]] = k_a;
	    else if(dom==1)
	    {
	       	if (a[NoSS_k[m]][NoSS_l[m]] < 0.0)
			ha[NoSS_k[m]][NoSS_l[m]] = 0.5 * 
			exp(k_a*a[NoSS_k[m]][NoSS_l[m]]);
	    		else
			ha[NoSS_k[m]][NoSS_l[m]] = 0.5 *
			exp(-k_a*a[NoSS_k[m]][NoSS_l[m]]);
    	    }
	    else
	    {
	       	if (a[NoSS_k[m]][NoSS_l[m]] < 0.0)
			ha[NoSS_k[m]][NoSS_l[m]] = uniform() * 
			exp(k_a*a[NoSS_k[m]][NoSS_l[m]]);
	    		else
			ha[NoSS_k[m]][NoSS_l[m]] = uniform() *
			exp(-k_a*a[NoSS_k[m]][NoSS_l[m]]);
    	    }
	    if(rho==99.0)	ha[NoSS_k[m]][NoSS_l[m]] = hs[NoSS_k[m]][NoSS_l[m]]; 

	    if(zs==1.0)    ha[NoSS_k[m]][NoSS_l[m]] = 1.0;

	    mean_a[(gen+1)/block] += za;
	    mean_a2[(gen+1)/block] += za*za;
	    mean_a3[(gen+1)/block] +=  za*za*za;
	    mean_a4[(gen+1)/block] +=  za*za*za*za;
	    if(zs!=1.0) 	mean_ha[(gen+1)/block] += ha[NoSS_k[m]][NoSS_l[m]]; 
	    else    	mean_ha[(gen+1)/block] += ave_ha; 

	 

	   if (tracelevel!=0)    fprintf(fptr,"s=%f  a=%f\n", s[NoSS_k[m]][NoSS_l[m]], a[NoSS_k[m]][NoSS_l[m]]);
	   if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  k=%d  l=%d   ran_h=%d\n", ran_i, NoSS_k[m],                                        	   NoSS_l[m], ran_h);
	}
    }

    for (m=mutants_ocurred; m<muts; m++)
    {
	totalmutant[(gen+1)/block] +=  1;
	ran_i = (int)(uniform()*NIND);
	ran_k = (int)(uniform()*NCRO);
	do {ran_l = (int)(uniform()*NLOCI);}   while (ran_l== 0);
	ran_h = (int)(uniform()*2.0);
	if ( (gm[ran_i][ran_k][ran_h] & RM[ran_l])==RM[ran_l] )
	{
	    doublemutant[(gen+1)/block] += 1;
	}
	else
	{
	    gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[ran_l]);
	    recurrentmutant[(gen+1)/block] += 1;
	}
	if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  ran_k=%d  ran_l=%d  ran_h=%d\n", ran_i, ran_k, ran_l, ran_h);
    }
    label: /* end of mutations */;
    if (tracelevel!=0)    fprintf(fptr,"doublemutants = %d   recurrentmutants = %d\n"
		, doublemutant[(gen+1)/block], recurrentmutant[(gen+1)/block]);
}


/* ***************************************************************** */


int mutationnumber ()
{
	int r;
	if ((2.0*(double)NIND*Lambda_a < normalthreshold) && (exp(-2.0*(double)NIND*Lambda_a) != 0.0) )
	{
		r= poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_a, sqrt(2.0*(double)NIND*Lambda_a)) );
	return(r);
}


/* ************************************************************** */


int mutationnumberL ()
{
	int r;
	if ((2.0*(double)NIND*Lambda_L < normalthreshold) &&
(exp(-2.0*(double)NIND*Lambda_L) != 0.0) )
	{
		r= poisson(lastinmutantspoissontableL, mutantspoissontableL);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_L, 
sqrt(2.0*(double)NIND*Lambda_L)) );
	return(r);
}


/* *************************************************************** */


void disorder_NoSS (NoSS_k,NoSS_l)
int NoSS_k[], NoSS_l[];
{
	int a, b, rnd;
	
	for (i=0; i<countNoSS-1; i++)
	{
	   rnd=(int)(uniform()*(countNoSS-i));
	   a=NoSS_k[countNoSS-1-i]; NoSS_k[countNoSS-1-i]=NoSS_k[rnd]; NoSS_k[rnd]=a;
	   b=NoSS_l[countNoSS-1-i]; NoSS_l[countNoSS-1-i]=NoSS_l[rnd]; NoSS_l[rnd]=b;
	}
}


/* ************************************************************** */


dumpoffspringaftermutation()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Offspring after mutation (gm0 gm1)\n");	
	for (i=0; i<NIND; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}


/* **************************************************************** */


printout()
{ 
    if ( (fixn[gen/block]==0.0) && (lostn[gen/block]==0.0) )
    {
	Ave_CHn = 0.0;
	Ave_NSSn = 0.0;
    }
    else
    {
	Ave_CHn = CHn[gen/block]/(lostn[gen/block]+fixn[gen/block]);
	Ave_NSSn = (timetolossn[gen/block]+timetofixn[gen/block]) / (lostn[gen/block] + 				fixn[gen/block]);
    }

    if (lostn[gen/block]==0.0)    Ave_timetolossn = 0.0;
    else    Ave_timetolossn = timetolossn[gen/block]/lostn[gen/block];

    if (fixn[gen/block]==0.0)    Ave_timetofixn = 0.0;
    else    Ave_timetofixn = timetofixn[gen/block]/fixn[gen/block];

    if (other[gen/block]==0.0)    	Ave_freq = 0.0;
    else    Ave_freq = freq[gen/block] / other[gen/block];

    if ( (fix[gen/block]==0.0) && (lost[gen/block]==0.0) )
    {
	Ave_CH = 0.0;
	Ave_NSS = 0.0;
    }
    else
    {
	Ave_CH = CH[gen/block]/(lost[gen/block]+fix[gen/block]);
	Ave_NSS = (timetoloss[gen/block]+timetofix[gen/block]) / (lost[gen/block]+fix[gen/block]);
    }

    if (lost[gen/block]==0.0)    Ave_timetoloss = 0.0;
    else    Ave_timetoloss = timetoloss[gen/block]/lost[gen/block];

    if (fix[gen/block]==0.0)    Ave_timetofix = 0.0;
    else    Ave_timetofix = timetofix[gen/block]/fix[gen/block];

    fgen = fopen ("genfile","a");

    fprintf(fgen,"\n\n\n gen=%d  N_SO=%4.2f+-%4.2f     SK2=%4.2f+-%4.2f    N_Recom=%4.2f\n", gen+1, accmean(&num_SO[gen/block]), se(&num_SO[gen/block]), accmean(&SK2[gen/block]), se(&SK2[gen/block]), (double)numberofrecom[gen/block]/(2.0*(double)NIND*(double)block));

fprintf(fgen, "         gmean_s=%6.4f+-%6.4f      gvar_s=%6.4f+-%6.4f    C2=%6.4f+-%6.4f \n         gmean_a=%6.4f+-%6.4f     gvar_a=%6.4f+-%6.4f    H2_a=%6.4f \n", accmean(&gmean_s[gen/block]), se(&gmean_s[gen/block]), accmean(&gvar_s[gen/block]), se(&gvar_s[gen/block]), accmean(&C2[gen/block]), se(&C2[gen/block]), accmean(&gmean_a[gen/block]), se(&gmean_a[gen/block]), accmean(&gvar_a[gen/block]), se(&gvar_a[gen/block]), accmean(&gvar_a[gen/block]) / (accmean(&gvar_a[gen/block]) + VE)); 

fprintf(fgen, "         VA_s=%6.4f+-%6.4f      VD_s=%6.4f+-%6.4f\n", accmean(&VA_s[gen/block]), se(&VA_s[gen/block]), accmean(&VD_s[gen/block]), se(&VD_s[gen/block])); 

h2a = accmean(&Ave_VA_all) /
	(  accmean(&Ave_VA_all)  +  accmean(&Ave_VD_all)  +  VE ) ;

fprintf(fgen, "         VA_a=%6.4f+-%6.4f      VD_a=%6.4f+-%6.4f         h2a=%6.4f\n         rA=%6.4f+-%6.4f      b_SS=%f+-%f     Vs(obs)=%6.4f \n", accmean(&Ave_VA_all), se(&Ave_VA_all), accmean(&Ave_VD_all), se(&Ave_VD_all), h2a, accmean(&rA[gen/block]), se(&rA[gen/block]), accmean(&b_SS[gen/block]), se(&b_SS[gen/block]), (-1.0) / (2.0*accmean(&b_SS[gen/block]))); 

    fprintf(fgen, "\n         NEU  lost=%4.2f   fix=%4.2f   other=%4.2f   NSS=%4.2f   H=%f\n         t_loss=%4.2f  t_fix=%4.2f   q=%f\n"
	,lostn[gen/block]/(double)block,fixn[gen/block]/(double)block,othern[gen/block]/(double)block
		,Ave_NSSn,Ave_CHn,Ave_timetolossn,Ave_timetofixn
			,freqn[gen/block]/((double)block*(double)NCRO));

    fprintf(fgen, "\n         SEL  lost=%4.2f   fix=%4.2f   other=%4.2f   NSS=%4.2f   H=%f\n         t_loss=%4.2f   t_fix=%4.2f   q=%f   rAB=%f\n\n", 	lost[gen/block]/(double)block,fix[gen/block]/(double)block,other[gen/block]/(double)block
		,Ave_NSS,Ave_CH,Ave_timetoloss,Ave_timetofix,Ave_freq,accmean(&rAB[gen/block]));

    if(meanmuts_s[gen/block] != 0.0)
    {
   	h1_s=mean_s[gen/block] / meanmuts_s[gen/block];
    	h2_s=mean_s2[gen/block] / meanmuts_s[gen/block];
    	h3_s=mean_s3[gen/block] / meanmuts_s[gen/block];
    	h4_s=mean_s4[gen/block] / meanmuts_s[gen/block];
	h1_a_s=mean_a_s[gen/block] / meanmuts_s[gen/block];
    	h2_a_s=mean_a2_s[gen/block] / meanmuts_s[gen/block];
	h_sa=mean_sa[gen/block] / meanmuts_s[gen/block];
    }
    if(meanmuts_a[gen/block] != 0.0)
    {
	h1_a=mean_a[gen/block] / meanmuts_a[gen/block];
    	h2_a=mean_a2[gen/block] / meanmuts_a[gen/block];
    	h3_a=mean_a3[gen/block] / meanmuts_a[gen/block];
    	h4_a=mean_a4[gen/block] / meanmuts_a[gen/block];
    }
    fprintf (fgen, "         Ave. number of mut. per gen.(QT): %f\n", 			(double)totalmutant[gen/block]/(double)block);
    fprintf (fgen, "         Ave. number of recurrent mutants : %f\n", 	(double)recurrentmutant[gen/block]/(double)block);
    fprintf (fgen, "         Ave. number of double mutants : %f\n", 	
	(double)doublemutant[gen/block]/(double)block);
    fprintf (fgen, "         Ave. number of lethal mutants : %f\n", 	(double)meanmutsL_s[gen/block]/(double)block);

    fprintf (fgen, "         Mean  a : %f\n", mean_real_a[gen/block] / meanmuts_a[gen/block]);
    fprintf (fgen, "         Mean |a| : %f\n", h1_a);
    fprintf (fgen, "         Variance |a| : %f\n", h2_a - (h1_a*h1_a));
    fprintf (fgen, "         Kurtosis |a| : %f\n", ((h4_a - 4.0*h1_a*h3_a + 6.0*h1_a*h1_a*h2_a - 	3.0*h1_a*h1_a*h1_a*h1_a) / ((h2_a - h1_a*h1_a) * (h2_a - h1_a*h1_a)))-3.0);
    if(meanmuts_a[gen/block] != 0.0)
          fprintf (fgen, "         Mean ha : %f\n", mean_ha[gen/block] / meanmuts_a[gen/block]);

    fprintf (fgen, "         Ave. number of mut. per gen.(on fitness, excl. lethals): %f\n", 	meanmuts_s[gen/block]/(double)block);
    fprintf (fgen, "         Mean |s| : %f\n", h1_s);
    fprintf (fgen, "         Variance |s| : %f\n", h2_s - (h1_s*h1_s));
    fprintf (fgen, "         Kurtosis |s| : %f\n", ((h4_s - 4.0*h1_s*h3_s + 6.0*h1_s*h1_s*h2_s - 	3.0*h1_s*h1_s*h1_s*h1_s) / ((h2_s - h1_s*h1_s) * (h2_s - h1_s*h1_s)))-3.0);
    if(meanmuts_s[gen/block] != 0.0)
          fprintf (fgen, "         Mean hs : %f\n", mean_hs[gen/block] / meanmuts_s[gen/block]);

    fprintf (fgen, "         Cov (|a|,|s|) : %f\n", h_sa-(h1_s*h1_a_s));
    if ( (h2_s-(h1_s*h1_s)==0) || (h2_a_s-(h1_a_s*h1_a_s)==0) )
		fprintf (fgen, "         Correlation (|a|,|s|) : 0\n\n");
    else		fprintf (fgen, "         Correlation (|a|,|s|) : %f\n\n", (h_sa-(h1_s*h1_a_s)) / sqrt(h2_s-(			h1_s*h1_s)) / sqrt(h2_a_s-(h1_a_s*h1_a_s)));


    /* ****************** TOTAL VALUES OF DISTRIBUTION  ************************ */

fprintf (fgen, "\n             n.genes       s                    |a|             hs                  ha\n");

    if(accsum(&Ave_ng_all) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "all values %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n\n", accmean(&Ave_ng_all), se(&Ave_ng_all), accmean(&Ave_s_all), se(&Ave_s_all), accmean(&Ave_a_all), se(&Ave_a_all), accmean(&Ave_hs_all), se(&Ave_hs_all), accmean(&Ave_ha_all), se(&Ave_ha_all));

/* ************************************************************ */

fprintf (fgen, "\n             n.genes       s                    |a|             hs                  ha\n");

/* **************** DISTRIBUTION FOR GENE FREQUENCY ********************** */

if(accsum(&Ave_q_ng_00_01) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.0-0.1 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_00_01), se(&Ave_q_ng_00_01), accmean(&Ave_q_s_00_01), se(&Ave_q_s_00_01), accmean(&Ave_q_a_00_01), se(&Ave_q_a_00_01), accmean(&Ave_q_hs_00_01), se(&Ave_q_hs_00_01), accmean(&Ave_q_ha_00_01), se(&Ave_q_ha_00_01));

if(accsum(&Ave_q_ng_01_02) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.1-0.2 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_01_02), se(&Ave_q_ng_01_02), accmean(&Ave_q_s_01_02), se(&Ave_q_s_01_02), accmean(&Ave_q_a_01_02), se(&Ave_q_a_01_02), accmean(&Ave_q_hs_01_02), se(&Ave_q_hs_01_02), accmean(&Ave_q_ha_01_02), se(&Ave_q_ha_01_02));

if(accsum(&Ave_q_ng_02_03) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.2-0.3 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_02_03), se(&Ave_q_ng_02_03), accmean(&Ave_q_s_02_03), se(&Ave_q_s_02_03), accmean(&Ave_q_a_02_03), se(&Ave_q_a_02_03), accmean(&Ave_q_hs_02_03), se(&Ave_q_hs_02_03), accmean(&Ave_q_ha_02_03), se(&Ave_q_ha_02_03));

if(accsum(&Ave_q_ng_03_04) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.3-0.4 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_03_04), se(&Ave_q_ng_03_04), accmean(&Ave_q_s_03_04), se(&Ave_q_s_03_04), accmean(&Ave_q_a_03_04), se(&Ave_q_a_03_04), accmean(&Ave_q_hs_03_04), se(&Ave_q_hs_03_04), accmean(&Ave_q_ha_03_04), se(&Ave_q_ha_03_04));

if(accsum(&Ave_q_ng_04_05) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.4-0.5 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_04_05), se(&Ave_q_ng_04_05), accmean(&Ave_q_s_04_05), se(&Ave_q_s_04_05), accmean(&Ave_q_a_04_05), se(&Ave_q_a_04_05), accmean(&Ave_q_hs_04_05), se(&Ave_q_hs_04_05), accmean(&Ave_q_ha_04_05), se(&Ave_q_ha_04_05));

if(accsum(&Ave_q_ng_05_06) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.5-0.6 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_05_06), se(&Ave_q_ng_05_06), accmean(&Ave_q_s_05_06), se(&Ave_q_s_05_06), accmean(&Ave_q_a_05_06), se(&Ave_q_a_05_06), accmean(&Ave_q_hs_05_06), se(&Ave_q_hs_05_06), accmean(&Ave_q_ha_05_06), se(&Ave_q_ha_05_06));

if(accsum(&Ave_q_ng_06_07) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.6-0.7 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_06_07), se(&Ave_q_ng_06_07), accmean(&Ave_q_s_06_07), se(&Ave_q_s_06_07), accmean(&Ave_q_a_06_07), se(&Ave_q_a_06_07), accmean(&Ave_q_hs_06_07), se(&Ave_q_hs_06_07), accmean(&Ave_q_ha_06_07), se(&Ave_q_ha_06_07));

if(accsum(&Ave_q_ng_07_08) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.7-0.8 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_07_08), se(&Ave_q_ng_07_08), accmean(&Ave_q_s_07_08), se(&Ave_q_s_07_08), accmean(&Ave_q_a_07_08), se(&Ave_q_a_07_08), accmean(&Ave_q_hs_07_08), se(&Ave_q_hs_07_08), accmean(&Ave_q_ha_07_08), se(&Ave_q_ha_07_08));

if(accsum(&Ave_q_ng_08_09) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.8-0.9 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_08_09), se(&Ave_q_ng_08_09), accmean(&Ave_q_s_08_09), se(&Ave_q_s_08_09), accmean(&Ave_q_a_08_09), se(&Ave_q_a_08_09), accmean(&Ave_q_hs_08_09), se(&Ave_q_hs_08_09), accmean(&Ave_q_ha_08_09), se(&Ave_q_ha_08_09));

if(accsum(&Ave_q_ng_09_10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.9-1.0 %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_q_ng_09_10), se(&Ave_q_ng_09_10), accmean(&Ave_q_s_09_10), se(&Ave_q_s_09_10), accmean(&Ave_q_a_09_10), se(&Ave_q_a_09_10), accmean(&Ave_q_hs_09_10), se(&Ave_q_hs_09_10), accmean(&Ave_q_ha_09_10), se(&Ave_q_ha_09_10));

/* ********************* DISTRIBUTION FOR GENE EFFECTS  ******************** */

    fprintf (fgen, "\n               n.genes       s                    q                   ha\n");
   
    /* ************************** negative values ************************* */

if(accmean(&Ave_a_ng_m20) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-inf-2.0  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m20), se(&Ave_a_ng_m20), accmean(&Ave_a_s_m20), se(&Ave_a_s_m20), accmean(&Ave_a_q_m20), se(&Ave_a_q_m20), accmean(&Ave_a_ha_m20), se(&Ave_a_ha_m20));

if(accmean(&Ave_a_ng_m18_m20) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-2.0-1.8  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m18_m20), se(&Ave_a_ng_m18_m20), accmean(&Ave_a_s_m18_m20), se(&Ave_a_s_m18_m20), accmean(&Ave_a_q_m18_m20), se(&Ave_a_q_m18_m20), accmean(&Ave_a_ha_m18_m20), se(&Ave_a_ha_m18_m20));

if(accmean(&Ave_a_ng_m16_m18) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-1.8-1.6  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m16_m18), se(&Ave_a_ng_m16_m18), accmean(&Ave_a_s_m16_m18), se(&Ave_a_s_m16_m18), accmean(&Ave_a_q_m16_m18), se(&Ave_a_q_m16_m18), accmean(&Ave_a_ha_m16_m18), se(&Ave_a_ha_m16_m18));

if(accmean(&Ave_a_ng_m14_m16) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-1.6-1.4  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m14_m16), se(&Ave_a_ng_m14_m16), accmean(&Ave_a_s_m14_m16), se(&Ave_a_s_m14_m16), accmean(&Ave_a_q_m14_m16), se(&Ave_a_q_m14_m16), accmean(&Ave_a_ha_m14_m16), se(&Ave_a_ha_m14_m16));

if(accmean(&Ave_a_ng_m12_m14) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-1.4-1.2  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m12_m14), se(&Ave_a_ng_m12_m14), accmean(&Ave_a_s_m12_m14), se(&Ave_a_s_m12_m14), accmean(&Ave_a_q_m12_m14), se(&Ave_a_q_m12_m14), accmean(&Ave_a_ha_m12_m14), se(&Ave_a_ha_m12_m14));

if(accmean(&Ave_a_ng_m10_m12) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-1.2-1.0  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m10_m12), se(&Ave_a_ng_m10_m12), accmean(&Ave_a_s_m10_m12), se(&Ave_a_s_m10_m12), accmean(&Ave_a_q_m10_m12), se(&Ave_a_q_m10_m12), accmean(&Ave_a_ha_m10_m12), se(&Ave_a_ha_m10_m12));

if(accmean(&Ave_a_ng_m08_m10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-1.0-0.8  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m08_m10), se(&Ave_a_ng_m08_m10), accmean(&Ave_a_s_m08_m10), se(&Ave_a_s_m08_m10), accmean(&Ave_a_q_m08_m10), se(&Ave_a_q_m08_m10), accmean(&Ave_a_ha_m08_m10), se(&Ave_a_ha_m08_m10));

if(accmean(&Ave_a_ng_m06_m08) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-0.8-0.6  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m06_m08), se(&Ave_a_ng_m06_m08), accmean(&Ave_a_s_m06_m08), se(&Ave_a_s_m06_m08), accmean(&Ave_a_q_m06_m08), se(&Ave_a_q_m06_m08), accmean(&Ave_a_ha_m06_m08), se(&Ave_a_ha_m06_m08));

if(accmean(&Ave_a_ng_m04_m06) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-0.6-0.4  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m04_m06), se(&Ave_a_ng_m04_m06), accmean(&Ave_a_s_m04_m06), se(&Ave_a_s_m04_m06), accmean(&Ave_a_q_m04_m06), se(&Ave_a_q_m04_m06), accmean(&Ave_a_ha_m04_m06), se(&Ave_a_ha_m04_m06));

if(accmean(&Ave_a_ng_m02_m04) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-0.4-0.2  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m02_m04), se(&Ave_a_ng_m02_m04), accmean(&Ave_a_s_m02_m04), se(&Ave_a_s_m02_m04), accmean(&Ave_a_q_m02_m04), se(&Ave_a_q_m02_m04), accmean(&Ave_a_ha_m02_m04), se(&Ave_a_ha_m02_m04));

if(accmean(&Ave_a_ng_m00_m02) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a=-0.2-0.0  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_m00_m02), se(&Ave_a_ng_m00_m02), accmean(&Ave_a_s_m00_m02), se(&Ave_a_s_m00_m02), accmean(&Ave_a_q_m00_m02), se(&Ave_a_q_m00_m02), accmean(&Ave_a_ha_m00_m02), se(&Ave_a_ha_m00_m02));

/* *********************** positive values ************************ */

if(accmean(&Ave_a_ng_00_02) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 0.0-0.2  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_00_02), se(&Ave_a_ng_00_02), accmean(&Ave_a_s_00_02), se(&Ave_a_s_00_02), accmean(&Ave_a_q_00_02), se(&Ave_a_q_00_02), accmean(&Ave_a_ha_00_02), se(&Ave_a_ha_00_02));

if(accmean(&Ave_a_ng_02_04) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 0.2-0.4  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_02_04), se(&Ave_a_ng_02_04), accmean(&Ave_a_s_02_04), se(&Ave_a_s_02_04), accmean(&Ave_a_q_02_04), se(&Ave_a_q_02_04), accmean(&Ave_a_ha_02_04), se(&Ave_a_ha_02_04));

if(accmean(&Ave_a_ng_04_06) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 0.4-0.6  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_04_06), se(&Ave_a_ng_04_06), accmean(&Ave_a_s_04_06), se(&Ave_a_s_04_06), accmean(&Ave_a_q_04_06), se(&Ave_a_q_04_06), accmean(&Ave_a_ha_04_06), se(&Ave_a_ha_04_06));

if(accmean(&Ave_a_ng_06_08) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 0.6-0.8  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_06_08), se(&Ave_a_ng_06_08), accmean(&Ave_a_s_06_08), se(&Ave_a_s_06_08), accmean(&Ave_a_q_06_08), se(&Ave_a_q_06_08), accmean(&Ave_a_ha_06_08), se(&Ave_a_ha_06_08));

if(accmean(&Ave_a_ng_08_10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 0.8-1.0  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_08_10), se(&Ave_a_ng_08_10), accmean(&Ave_a_s_08_10), se(&Ave_a_s_08_10), accmean(&Ave_a_q_08_10), se(&Ave_a_q_08_10), accmean(&Ave_a_ha_08_10), se(&Ave_a_ha_08_10));

if(accmean(&Ave_a_ng_10_12) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 1.0-1.2  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_10_12), se(&Ave_a_ng_10_12), accmean(&Ave_a_s_10_12), se(&Ave_a_s_10_12), accmean(&Ave_a_q_10_12), se(&Ave_a_q_10_12), accmean(&Ave_a_ha_10_12), se(&Ave_a_ha_10_12));

if(accmean(&Ave_a_ng_12_14) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 1.2-1.4  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_12_14), se(&Ave_a_ng_12_14), accmean(&Ave_a_s_12_14), se(&Ave_a_s_12_14), accmean(&Ave_a_q_12_14), se(&Ave_a_q_12_14), accmean(&Ave_a_ha_12_14), se(&Ave_a_ha_12_14));

if(accmean(&Ave_a_ng_14_16) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 1.4-1.6  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_14_16), se(&Ave_a_ng_14_16), accmean(&Ave_a_s_14_16), se(&Ave_a_s_14_16), accmean(&Ave_a_q_14_16), se(&Ave_a_q_14_16), accmean(&Ave_a_ha_14_16), se(&Ave_a_ha_14_16));

if(accmean(&Ave_a_ng_16_18) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 1.6-1.8  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_16_18), se(&Ave_a_ng_16_18), accmean(&Ave_a_s_16_18), se(&Ave_a_s_16_18), accmean(&Ave_a_q_16_18), se(&Ave_a_q_16_18), accmean(&Ave_a_ha_16_18), se(&Ave_a_ha_16_18));

if(accmean(&Ave_a_ng_18_20) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 1.8-2.0  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_18_20), se(&Ave_a_ng_18_20), accmean(&Ave_a_s_18_20), se(&Ave_a_s_18_20), accmean(&Ave_a_q_18_20), se(&Ave_a_q_18_20), accmean(&Ave_a_ha_18_20), se(&Ave_a_ha_18_20));

if(accmean(&Ave_a_ng_20) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "a= 2.0-inf  %6.2f+-%4.2f  %f+-%f  %f+-%f  %f+-%f\n", accmean(&Ave_a_ng_20), se(&Ave_a_ng_20), accmean(&Ave_a_s_20), se(&Ave_a_s_20), accmean(&Ave_a_q_20), se(&Ave_a_q_20), accmean(&Ave_a_ha_20), se(&Ave_a_ha_20));

/* ********************* DISTRIBUTION FOR DOMINANCE  ********************** */

    fprintf (fgen, "\n               n.genes      q                  hs\n");
   
if(accmean(&Ave_ha_ng_00_02) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "ha=0.0-0.2  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_ha_ng_00_02), se(&Ave_ha_ng_00_02), accmean(&Ave_ha_q_00_02), se(&Ave_ha_q_00_02), accmean(&Ave_ha_hs_00_02), se(&Ave_ha_hs_00_02));

if(accmean(&Ave_ha_ng_02_04) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "ha=0.2-0.4  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_ha_ng_02_04), se(&Ave_ha_ng_02_04), accmean(&Ave_ha_q_02_04), se(&Ave_ha_q_02_04), accmean(&Ave_ha_hs_02_04), se(&Ave_ha_hs_02_04));

if(accmean(&Ave_ha_ng_04_06) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "ha=0.4-0.6  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_ha_ng_04_06), se(&Ave_ha_ng_04_06), accmean(&Ave_ha_q_04_06), se(&Ave_ha_q_04_06), accmean(&Ave_ha_hs_04_06), se(&Ave_ha_hs_04_06));

if(accmean(&Ave_ha_ng_06_08) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "ha=0.6-0.8  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_ha_ng_06_08), se(&Ave_ha_ng_06_08), accmean(&Ave_ha_q_06_08), se(&Ave_ha_q_06_08), accmean(&Ave_ha_hs_06_08), se(&Ave_ha_hs_06_08));

if(accmean(&Ave_ha_ng_08_10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "ha=0.8-1.0  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_ha_ng_08_10), se(&Ave_ha_ng_08_10), accmean(&Ave_ha_q_08_10), se(&Ave_ha_q_08_10), accmean(&Ave_ha_hs_08_10), se(&Ave_ha_hs_08_10));

/* ************************ DISTRIBUTION FOR FITNESS  ********************** */

    fprintf (fgen, "\n                 n.genes     q                   hs\n");
   
if(accmean(&Ave_s_ng_0) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=0.0         %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_0), se(&Ave_s_ng_0), accmean(&Ave_s_q_0), se(&Ave_s_q_0), accmean(&Ave_s_hs_0), se(&Ave_s_hs_0));

if(accmean(&Ave_s_ng_0_106) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=0.0-10-6    %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_0_106), se(&Ave_s_ng_0_106), accmean(&Ave_s_q_0_106), se(&Ave_s_q_0_106), accmean(&Ave_s_hs_0_106), se(&Ave_s_hs_0_106));

if(accmean(&Ave_s_ng_106_104) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-10-6-10-4  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_106_104), se(&Ave_s_ng_106_104), accmean(&Ave_s_q_106_104), se(&Ave_s_q_106_104), accmean(&Ave_s_hs_106_104), se(&Ave_s_hs_106_104));

if(accmean(&Ave_s_ng_104_102) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-10-4-10-2  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_104_102), se(&Ave_s_ng_104_102), accmean(&Ave_s_q_104_102), se(&Ave_s_q_104_102), accmean(&Ave_s_hs_104_102), se(&Ave_s_hs_104_102));

if(accmean(&Ave_s_ng_102_101) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-10-2-10-1  %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_102_101), se(&Ave_s_ng_102_101), accmean(&Ave_s_q_102_101), se(&Ave_s_q_102_101), accmean(&Ave_s_hs_102_101), se(&Ave_s_hs_102_101));

if(accmean(&Ave_s_ng_01_02) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-0.1-0.2    %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_01_02), se(&Ave_s_ng_01_02), accmean(&Ave_s_q_01_02), se(&Ave_s_q_01_02), accmean(&Ave_s_hs_01_02), se(&Ave_s_hs_01_02));

if(accmean(&Ave_s_ng_02_04) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-0.2-0.4    %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_02_04), se(&Ave_s_ng_02_04), accmean(&Ave_s_q_02_04), se(&Ave_s_q_02_04), accmean(&Ave_s_hs_02_04), se(&Ave_s_hs_02_04));

if(accmean(&Ave_s_ng_04_06) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-0.4-0.6    %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_04_06), se(&Ave_s_ng_04_06), accmean(&Ave_s_q_04_06), se(&Ave_s_q_04_06), accmean(&Ave_s_hs_04_06), se(&Ave_s_hs_04_06));

if(accmean(&Ave_s_ng_06_10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-0.6-1.0    %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_06_10), se(&Ave_s_ng_06_10), accmean(&Ave_s_q_06_10), se(&Ave_s_q_06_10), accmean(&Ave_s_hs_06_10), se(&Ave_s_hs_06_10));

if(accmean(&Ave_s_ng_10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "s=-1.0        %6.2f+-%4.2f  %f+-%f  %f+-%f\n", accmean(&Ave_s_ng_10), se(&Ave_s_ng_10), accmean(&Ave_s_q_10), se(&Ave_s_q_10), accmean(&Ave_s_hs_10), se(&Ave_s_hs_10));

/* ****************** TOTAL VALUES OF DISTRIBUTION  ************************ */

fprintf (fgen, "\n               n.genes              VA                           VD\n");

    if(accmean(&Ave_ng_all) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "all values  %6.2f+-%4.2f           %f+-%f           %f+-%f\n\n", accmean(&Ave_ng_all), se(&Ave_ng_all), accmean(&Ave_VA_all), se(&Ave_VA_all), accmean(&Ave_VD_all), se(&Ave_VD_all));

/* **************************************************************** */

if(accmean(&Ave_q_ng_00_01) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.0-0.1  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_00_01), se(&Ave_q_ng_00_01), 100.0 * accmean(&Ave_q_ng_00_01) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_00_01), se(&Ave_q_VA_00_01), 100.0 * accmean(&Ave_q_VA_00_01) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_00_01), se(&Ave_q_VD_00_01), 100.0 * accmean(&Ave_q_VD_00_01) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_01_02) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.1-0.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_01_02), se(&Ave_q_ng_01_02), 100.0 * accmean(&Ave_q_ng_01_02) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_01_02), se(&Ave_q_VA_01_02), 100.0 * accmean(&Ave_q_VA_01_02) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_01_02), se(&Ave_q_VD_01_02), 100.0 * accmean(&Ave_q_VD_01_02) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_02_03) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.2-0.3  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_02_03), se(&Ave_q_ng_02_03), 100.0 * accmean(&Ave_q_ng_02_03) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_02_03), se(&Ave_q_VA_02_03), 100.0 * accmean(&Ave_q_VA_02_03) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_02_03), se(&Ave_q_VD_02_03), 100.0 * accmean(&Ave_q_VD_02_03) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_03_04) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.3-0.4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_03_04), se(&Ave_q_ng_03_04), 100.0 * accmean(&Ave_q_ng_03_04) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_03_04), se(&Ave_q_VA_03_04), 100.0 * accmean(&Ave_q_VA_03_04) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_03_04), se(&Ave_q_VD_03_04), 100.0 * accmean(&Ave_q_VD_03_04) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_04_05) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.4-0.5  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_04_05), se(&Ave_q_ng_04_05), 100.0 * accmean(&Ave_q_ng_04_05) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_04_05), se(&Ave_q_VA_04_05), 100.0 * accmean(&Ave_q_VA_04_05) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_04_05), se(&Ave_q_VD_04_05), 100.0 * accmean(&Ave_q_VD_04_05) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_05_06) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.5-0.6  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_05_06), se(&Ave_q_ng_05_06), 100.0 * accmean(&Ave_q_ng_05_06) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_05_06), se(&Ave_q_VA_05_06), 100.0 * accmean(&Ave_q_VA_05_06) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_05_06), se(&Ave_q_VD_05_06), 100.0 * accmean(&Ave_q_VD_05_06) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_06_07) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.6-0.7  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_06_07), se(&Ave_q_ng_06_07), 100.0 * accmean(&Ave_q_ng_06_07) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_06_07), se(&Ave_q_VA_06_07), 100.0 * accmean(&Ave_q_VA_06_07) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_06_07), se(&Ave_q_VD_06_07), 100.0 * accmean(&Ave_q_VD_06_07) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_07_08) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.7-0.8  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_07_08), se(&Ave_q_ng_07_08), 100.0 * accmean(&Ave_q_ng_07_08) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_07_08), se(&Ave_q_VA_07_08), 100.0 * accmean(&Ave_q_VA_07_08) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_07_08), se(&Ave_q_VD_07_08), 100.0 * accmean(&Ave_q_VD_07_08) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_08_09) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.8-0.9  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_08_09), se(&Ave_q_ng_08_09), 100.0 * accmean(&Ave_q_ng_08_09) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_08_09), se(&Ave_q_VA_08_09), 100.0 * accmean(&Ave_q_VA_08_09) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_08_09), se(&Ave_q_VD_08_09), 100.0 * accmean(&Ave_q_VD_08_09) / accmean(&Ave_VD_all));

if(accmean(&Ave_q_ng_09_10) == 0.0)		/* do not print */;
    else 	fprintf (fgen, "q=0.9-1.0  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_q_ng_09_10), se(&Ave_q_ng_09_10), 100.0 * accmean(&Ave_q_ng_09_10) / accmean(&Ave_ng_all), accmean(&Ave_q_VA_09_10), se(&Ave_q_VA_09_10), 100.0 * accmean(&Ave_q_VA_09_10) / accmean(&Ave_VA_all), accmean(&Ave_q_VD_09_10), se(&Ave_q_VD_09_10), 100.0 * accmean(&Ave_q_VD_09_10) / accmean(&Ave_VD_all));

fprintf (fgen, "\n");

/* ********************* DISTRIBUTION FOR GENE EFFECTS  ******************** */

/* ************************** negative values ************************* */

if(accmean(&Ave_a_ng_m20) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-inf-2.0  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m20), se(&Ave_a_ng_m20), 100.0 * accmean(&Ave_a_ng_m20) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m20), se(&Ave_a_VA_m20), 100.0 * accmean(&Ave_a_VA_m20) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m20), se(&Ave_a_VD_m20), 100.0 * accmean(&Ave_a_VD_m20) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m18_m20) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-2.0-1.8  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m18_m20), se(&Ave_a_ng_m18_m20), 100.0 * accmean(&Ave_a_ng_m18_m20) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m18_m20), se(&Ave_a_VA_m18_m20), 100.0 * accmean(&Ave_a_VA_m18_m20) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m18_m20), se(&Ave_a_VD_m18_m20), 100.0 * accmean(&Ave_a_VD_m18_m20) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m16_m18) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-1.8-1.6  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m16_m18), se(&Ave_a_ng_m16_m18), 100.0 * accmean(&Ave_a_ng_m16_m18) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m16_m18), se(&Ave_a_VA_m16_m18), 100.0 * accmean(&Ave_a_VA_m16_m18) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m16_m18), se(&Ave_a_VD_m16_m18), 100.0 * accmean(&Ave_a_VD_m16_m18) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m14_m16) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-1.6-1.4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m14_m16), se(&Ave_a_ng_m14_m16), 100.0 * accmean(&Ave_a_ng_m14_m16) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m14_m16), se(&Ave_a_VA_m14_m16), 100.0 * accmean(&Ave_a_VA_m14_m16) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m14_m16), se(&Ave_a_VD_m14_m16), 100.0 * accmean(&Ave_a_VD_m14_m16) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m12_m14) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-1.4-1.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m12_m14), se(&Ave_a_ng_m12_m14), 100.0 * accmean(&Ave_a_ng_m12_m14) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m12_m14), se(&Ave_a_VA_m12_m14), 100.0 * accmean(&Ave_a_VA_m12_m14) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m12_m14), se(&Ave_a_VD_m12_m14), 100.0 * accmean(&Ave_a_VD_m12_m14) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m10_m12) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-1.2-1.0  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m10_m12), se(&Ave_a_ng_m10_m12), 100.0 * accmean(&Ave_a_ng_m10_m12) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m10_m12), se(&Ave_a_VA_m10_m12), 100.0 * accmean(&Ave_a_VA_m10_m12) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m10_m12), se(&Ave_a_VD_m10_m12), 100.0 * accmean(&Ave_a_VD_m10_m12) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m08_m10) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-1.0-0.8  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m08_m10), se(&Ave_a_ng_m08_m10), 100.0 * accmean(&Ave_a_ng_m08_m10) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m08_m10), se(&Ave_a_VA_m08_m10), 100.0 * accmean(&Ave_a_VA_m08_m10) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m08_m10), se(&Ave_a_VD_m08_m10), 100.0 * accmean(&Ave_a_VD_m08_m10) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m06_m08) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-0.8-0.6  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m06_m08), se(&Ave_a_ng_m06_m08), 100.0 * accmean(&Ave_a_ng_m06_m08) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m06_m08), se(&Ave_a_VA_m06_m08), 100.0 * accmean(&Ave_a_VA_m06_m08) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m06_m08), se(&Ave_a_VD_m06_m08), 100.0 * accmean(&Ave_a_VD_m06_m08) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m04_m06) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-0.6-0.4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m04_m06), se(&Ave_a_ng_m04_m06), 100.0 * accmean(&Ave_a_ng_m04_m06) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m04_m06), se(&Ave_a_VA_m04_m06), 100.0 * accmean(&Ave_a_VA_m04_m06) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m04_m06), se(&Ave_a_VD_m04_m06), 100.0 * accmean(&Ave_a_VD_m04_m06) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m02_m04) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-0.4-0.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m02_m04), se(&Ave_a_ng_m02_m04), 100.0 * accmean(&Ave_a_ng_m02_m04) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m02_m04), se(&Ave_a_VA_m02_m04), 100.0 * accmean(&Ave_a_VA_m02_m04) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m02_m04), se(&Ave_a_VD_m02_m04), 100.0 * accmean(&Ave_a_VD_m02_m04) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_m00_m02) == 0.0)		/* do not print */;
    else fprintf (fgen, "a=-0.0-0.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_m00_m02), se(&Ave_a_ng_m00_m02), 100.0 * accmean(&Ave_a_ng_m00_m02) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_m00_m02), se(&Ave_a_VA_m00_m02), 100.0 * accmean(&Ave_a_VA_m00_m02) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_m00_m02), se(&Ave_a_VD_m00_m02), 100.0 * accmean(&Ave_a_VD_m00_m02) / accmean(&Ave_VD_all)); 

/* *********************** positive values ************************ */

if(accmean(&Ave_a_ng_00_02) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 0.0-0.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_00_02), se(&Ave_a_ng_00_02), 100.0 * accmean(&Ave_a_ng_00_02) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_00_02), se(&Ave_a_VA_00_02), 100.0 * accmean(&Ave_a_VA_00_02) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_00_02), se(&Ave_a_VD_00_02), 100.0 * accmean(&Ave_a_VD_00_02) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_02_04) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 0.2-0.4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_02_04), se(&Ave_a_ng_02_04), 100.0 * accmean(&Ave_a_ng_02_04) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_02_04), se(&Ave_a_VA_02_04), 100.0 * accmean(&Ave_a_VA_02_04) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_02_04), se(&Ave_a_VD_02_04), 100.0 * accmean(&Ave_a_VD_02_04) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_04_06) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 0.4-0.6  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_04_06), se(&Ave_a_ng_04_06), 100.0 * accmean(&Ave_a_ng_04_06) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_04_06), se(&Ave_a_VA_04_06), 100.0 * accmean(&Ave_a_VA_04_06) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_04_06), se(&Ave_a_VD_04_06), 100.0 * accmean(&Ave_a_VD_04_06) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_06_08) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 0.6-0.8  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_06_08), se(&Ave_a_ng_06_08), 100.0 * accmean(&Ave_a_ng_06_08) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_06_08), se(&Ave_a_VA_06_08), 100.0 * accmean(&Ave_a_VA_06_08) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_06_08), se(&Ave_a_VD_06_08), 100.0 * accmean(&Ave_a_VD_06_08) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_08_10) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 0.8-1.0  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_08_10), se(&Ave_a_ng_08_10), 100.0 * accmean(&Ave_a_ng_08_10) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_08_10), se(&Ave_a_VA_08_10), 100.0 * accmean(&Ave_a_VA_08_10) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_08_10), se(&Ave_a_VD_08_10), 100.0 * accmean(&Ave_a_VD_08_10) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_10_12) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 1.0-1.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_10_12), se(&Ave_a_ng_10_12), 100.0 * accmean(&Ave_a_ng_10_12) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_10_12), se(&Ave_a_VA_10_12), 100.0 * accmean(&Ave_a_VA_10_12) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_10_12), se(&Ave_a_VD_10_12), 100.0 * accmean(&Ave_a_VD_10_12) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_12_14) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 1.2-1.4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_12_14), se(&Ave_a_ng_12_14), 100.0 * accmean(&Ave_a_ng_12_14) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_12_14), se(&Ave_a_VA_12_14), 100.0 * accmean(&Ave_a_VA_12_14) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_12_14), se(&Ave_a_VD_12_14), 100.0 * accmean(&Ave_a_VD_12_14) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_14_16) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 1.4-1.6  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_14_16), se(&Ave_a_ng_14_16), 100.0 * accmean(&Ave_a_ng_14_16) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_14_16), se(&Ave_a_VA_14_16), 100.0 * accmean(&Ave_a_VA_14_16) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_14_16), se(&Ave_a_VD_14_16), 100.0 * accmean(&Ave_a_VD_14_16) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_16_18) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 1.6-1.8  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_16_18), se(&Ave_a_ng_16_18), 100.0 * accmean(&Ave_a_ng_16_18) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_16_18), se(&Ave_a_VA_16_18), 100.0 * accmean(&Ave_a_VA_16_18) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_16_18), se(&Ave_a_VD_16_18), 100.0 * accmean(&Ave_a_VD_16_18) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_18_20) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 1.8-2.0  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_18_20), se(&Ave_a_ng_18_20), 100.0 * accmean(&Ave_a_ng_18_20) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_18_20), se(&Ave_a_VA_18_20), 100.0 * accmean(&Ave_a_VA_18_20) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_18_20), se(&Ave_a_VD_18_20), 100.0 * accmean(&Ave_a_VD_18_20) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_a_ng_20) == 0.0)		/* do not print */;
    else fprintf (fgen, "a= 2.0-inf  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_a_ng_20), se(&Ave_a_ng_20), 100.0 * accmean(&Ave_a_ng_20) / accmean(&Ave_ng_all), accmean(&Ave_a_VA_20), se(&Ave_a_VA_20), 100.0 * accmean(&Ave_a_VA_20) / accmean(&Ave_VA_all), accmean(&Ave_a_VD_20), se(&Ave_a_VD_20), 100.0 * accmean(&Ave_a_VD_20) / accmean(&Ave_VD_all)); 

fprintf (fgen, "\n");

/* ********************* DISTRIBUTION FOR DOMINANCE  ********************** */

if(accmean(&Ave_ha_ng_00_02) == 0.0)		/* do not print */;
    else fprintf (fgen, "ha=0.0-0.2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_ha_ng_00_02), se(&Ave_ha_ng_00_02), 100.0 * accmean(&Ave_ha_ng_00_02) / accmean(&Ave_ng_all), accmean(&Ave_ha_VA_00_02), se(&Ave_ha_VA_00_02), 100.0 * accmean(&Ave_ha_VA_00_02) / accmean(&Ave_VA_all), accmean(&Ave_ha_VD_00_02), se(&Ave_ha_VD_00_02), 100.0 * accmean(&Ave_ha_VD_00_02) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_ha_ng_02_04) == 0.0)		/* do not print */;
    else fprintf (fgen, "ha=0.2-0.4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_ha_ng_02_04), se(&Ave_ha_ng_02_04), 100.0 * accmean(&Ave_ha_ng_02_04) / accmean(&Ave_ng_all), accmean(&Ave_ha_VA_02_04), se(&Ave_ha_VA_02_04), 100.0 * accmean(&Ave_ha_VA_02_04) / accmean(&Ave_VA_all), accmean(&Ave_ha_VD_02_04), se(&Ave_ha_VD_02_04), 100.0 * accmean(&Ave_ha_VD_02_04) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_ha_ng_04_06) == 0.0)		/* do not print */;
    else fprintf (fgen, "ha=0.4-0.6  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_ha_ng_04_06), se(&Ave_ha_ng_04_06), 100.0 * accmean(&Ave_ha_ng_04_06) / accmean(&Ave_ng_all), accmean(&Ave_ha_VA_04_06), se(&Ave_ha_VA_04_06), 100.0 * accmean(&Ave_ha_VA_04_06) / accmean(&Ave_VA_all), accmean(&Ave_ha_VD_04_06), se(&Ave_ha_VD_04_06), 100.0 * accmean(&Ave_ha_VD_04_06) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_ha_ng_06_08) == 0.0)		/* do not print */;
    else fprintf (fgen, "ha=0.6-0.8  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_ha_ng_06_08), se(&Ave_ha_ng_06_08), 100.0 * accmean(&Ave_ha_ng_06_08) / accmean(&Ave_ng_all), accmean(&Ave_ha_VA_06_08), se(&Ave_ha_VA_06_08), 100.0 * accmean(&Ave_ha_VA_06_08) / accmean(&Ave_VA_all), accmean(&Ave_ha_VD_06_08), se(&Ave_ha_VD_06_08), 100.0 * accmean(&Ave_ha_VD_06_08) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_ha_ng_08_10) == 0.0)		/* do not print */;
    else fprintf (fgen, "ha=0.8-1.0  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_ha_ng_08_10), se(&Ave_ha_ng_08_10), 100.0 * accmean(&Ave_ha_ng_08_10) / accmean(&Ave_ng_all), accmean(&Ave_ha_VA_08_10), se(&Ave_ha_VA_08_10), 100.0 * accmean(&Ave_ha_VA_08_10) / accmean(&Ave_VA_all), accmean(&Ave_ha_VD_08_10), se(&Ave_ha_VD_08_10), 100.0 * accmean(&Ave_ha_VD_08_10) / accmean(&Ave_VD_all)); 

fprintf (fgen, "\n");

/* ********************* DISTRIBUTION FOR FITNESS  ********************** */

if(accmean(&Ave_s_ng_0) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=0.0         %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_0), se(&Ave_s_ng_0), 100.0 * accmean(&Ave_s_ng_0) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_0), se(&Ave_s_VA_0), 100.0 * accmean(&Ave_s_VA_0) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_0), se(&Ave_s_VD_0), 100.0 * accmean(&Ave_s_VD_0) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_0_106) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=0.0-10-6    %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_0_106), se(&Ave_s_ng_0_106), 100.0 * accmean(&Ave_s_ng_0_106) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_0_106), se(&Ave_s_VA_0_106), 100.0 * accmean(&Ave_s_VA_0_106) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_0_106), se(&Ave_s_VD_0_106), 100.0 * accmean(&Ave_s_VD_0_106) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_106_104) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-10-6-10-4  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_106_104), se(&Ave_s_ng_106_104), 100.0 * accmean(&Ave_s_ng_106_104) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_106_104), se(&Ave_s_VA_106_104), 100.0 * accmean(&Ave_s_VA_106_104) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_106_104), se(&Ave_s_VD_106_104), 100.0 * accmean(&Ave_s_VD_106_104) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_104_102) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-10-4-10-2  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_104_102), se(&Ave_s_ng_104_102), 100.0 * accmean(&Ave_s_ng_104_102) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_104_102), se(&Ave_s_VA_104_102), 100.0 * accmean(&Ave_s_VA_104_102) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_104_102), se(&Ave_s_VD_104_102), 100.0 * accmean(&Ave_s_VD_104_102) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_102_101) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-10-2-10-1  %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_102_101), se(&Ave_s_ng_102_101), 100.0 * accmean(&Ave_s_ng_102_101) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_102_101), se(&Ave_s_VA_102_101), 100.0 * accmean(&Ave_s_VA_102_101) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_102_101), se(&Ave_s_VD_102_101), 100.0 * accmean(&Ave_s_VD_102_101) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_01_02) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-0.1-0.2    %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_01_02), se(&Ave_s_ng_01_02), 100.0 * accmean(&Ave_s_ng_01_02) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_01_02), se(&Ave_s_VA_01_02), 100.0 * accmean(&Ave_s_VA_01_02) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_01_02), se(&Ave_s_VD_01_02), 100.0 * accmean(&Ave_s_VD_01_02) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_02_04) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-0.2-0.4    %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_02_04), se(&Ave_s_ng_02_04), 100.0 * accmean(&Ave_s_ng_02_04) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_02_04), se(&Ave_s_VA_02_04), 100.0 * accmean(&Ave_s_VA_02_04) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_02_04), se(&Ave_s_VD_02_04), 100.0 * accmean(&Ave_s_VD_02_04) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_04_06) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-0.4-0.6    %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_04_06), se(&Ave_s_ng_04_06), 100.0 * accmean(&Ave_s_ng_04_06) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_04_06), se(&Ave_s_VA_04_06), 100.0 * accmean(&Ave_s_VA_04_06) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_04_06), se(&Ave_s_VD_04_06), 100.0 * accmean(&Ave_s_VD_04_06) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_06_10) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-0.6-1.0    %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_06_10), se(&Ave_s_ng_06_10), 100.0 * accmean(&Ave_s_ng_06_10) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_06_10), se(&Ave_s_VA_06_10), 100.0 * accmean(&Ave_s_VA_06_10) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_06_10), se(&Ave_s_VD_06_10), 100.0 * accmean(&Ave_s_VD_06_10) / accmean(&Ave_VD_all)); 

if(accmean(&Ave_s_ng_10) == 0.0)		/* do not print */;
    else fprintf (fgen, "s=-1.0        %6.2f+-%4.2f  (%5.2f)  %f+-%f  (%5.2f)  %f+-%f  (%5.2f)\n", accmean(&Ave_s_ng_10), se(&Ave_s_ng_10), 100.0 * accmean(&Ave_s_ng_10) / accmean(&Ave_ng_all), accmean(&Ave_s_VA_10), se(&Ave_s_VA_10), 100.0 * accmean(&Ave_s_VA_10) / accmean(&Ave_VA_all), accmean(&Ave_s_VD_10), se(&Ave_s_VD_10), 100.0 * accmean(&Ave_s_VD_10) / accmean(&Ave_VD_all)); 

/* ************************************************************* */

    fclose(fgen);
}


/* ************************************************************ */


settozero_gens_on_q()
{
    initacc(&ng_all); initacc(&s_all); initacc(&a_all); initacc(&hs_all); initacc(&ha_all); initacc(&VA_all); initacc(&VD_all);

	/* ****************************************************** */

    initacc(&q_ng_00_01); initacc(&q_s_00_01); initacc(&q_a_00_01); initacc(&q_hs_00_01); initacc(&q_ha_00_01); initacc(&q_VA_00_01); initacc(&q_VD_00_01);
	
    initacc(&q_ng_01_02); initacc(&q_s_01_02); initacc(&q_a_01_02); initacc(&q_hs_01_02); initacc(&q_ha_01_02); initacc(&q_VA_01_02); initacc(&q_VD_01_02);

    initacc(&q_ng_02_03); initacc(&q_s_02_03); initacc(&q_a_02_03); initacc(&q_hs_02_03); initacc(&q_ha_02_03); initacc(&q_VA_02_03); initacc(&q_VD_02_03);

    initacc(&q_ng_03_04); initacc(&q_s_03_04); initacc(&q_a_03_04); initacc(&q_hs_03_04); initacc(&q_ha_03_04); initacc(&q_VA_03_04); initacc(&q_VD_03_04);

    initacc(&q_ng_04_05); initacc(&q_s_04_05); initacc(&q_a_04_05); initacc(&q_hs_04_05); initacc(&q_ha_04_05); initacc(&q_VA_04_05); initacc(&q_VD_04_05);

    initacc(&q_ng_05_06); initacc(&q_s_05_06); initacc(&q_a_05_06); initacc(&q_hs_05_06); initacc(&q_ha_05_06); initacc(&q_VA_05_06); initacc(&q_VD_05_06);

    initacc(&q_ng_06_07); initacc(&q_s_06_07); initacc(&q_a_06_07); initacc(&q_hs_06_07); initacc(&q_ha_06_07); initacc(&q_VA_06_07); initacc(&q_VD_06_07);

    initacc(&q_ng_07_08); initacc(&q_s_07_08); initacc(&q_a_07_08); initacc(&q_hs_07_08); initacc(&q_ha_07_08); initacc(&q_VA_07_08); initacc(&q_VD_07_08);

    initacc(&q_ng_08_09); initacc(&q_s_08_09); initacc(&q_a_08_09); initacc(&q_hs_08_09); initacc(&q_ha_08_09); initacc(&q_VA_08_09); initacc(&q_VD_08_09);

    initacc(&q_ng_09_10); initacc(&q_s_09_10); initacc(&q_a_09_10); initacc(&q_hs_09_10); initacc(&q_ha_09_10); initacc(&q_VA_09_10); initacc(&q_VD_09_10);
}


/* ************************************************************** */


settozero_gens_on_a()
{
	initacc(&a_ng_00_02); initacc(&a_s_00_02); initacc(&a_q_00_02);
	initacc(&a_ha_00_02); initacc(&a_VA_00_02); initacc(&a_VD_00_02);

	initacc(&a_ng_02_04); initacc(&a_s_02_04); initacc(&a_q_02_04);
	initacc(&a_ha_02_04); initacc(&a_VA_02_04); initacc(&a_VD_02_04);

	initacc(&a_ng_04_06); initacc(&a_s_04_06); initacc(&a_q_04_06);
	initacc(&a_ha_04_06); initacc(&a_VA_04_06); initacc(&a_VD_04_06);

	initacc(&a_ng_06_08); initacc(&a_s_06_08); initacc(&a_q_06_08);
	initacc(&a_ha_06_08); initacc(&a_VA_06_08); initacc(&a_VD_06_08);

	initacc(&a_ng_08_10); initacc(&a_s_08_10); initacc(&a_q_08_10);
	initacc(&a_ha_08_10); initacc(&a_VA_08_10); initacc(&a_VD_08_10);

	initacc(&a_ng_10_12); initacc(&a_s_10_12); initacc(&a_q_10_12);
	initacc(&a_ha_10_12); initacc(&a_VA_10_12); initacc(&a_VD_10_12);

	initacc(&a_ng_12_14); initacc(&a_s_12_14); initacc(&a_q_12_14);
	initacc(&a_ha_12_14); initacc(&a_VA_12_14); initacc(&a_VD_12_14);

	initacc(&a_ng_14_16); initacc(&a_s_14_16); initacc(&a_q_14_16);
	initacc(&a_ha_14_16); initacc(&a_VA_14_16); initacc(&a_VD_14_16);

	initacc(&a_ng_16_18); initacc(&a_s_16_18); initacc(&a_q_16_18);
	initacc(&a_ha_16_18); initacc(&a_VA_16_18); initacc(&a_VD_16_18);

	initacc(&a_ng_18_20); initacc(&a_s_18_20); initacc(&a_q_18_20);
	initacc(&a_ha_18_20); initacc(&a_VA_18_20); initacc(&a_VD_18_20);

	initacc(&a_ng_20); initacc(&a_s_20); initacc(&a_q_20);
	initacc(&a_ha_20); initacc(&a_VA_20); initacc(&a_VD_20);

	/* ****************************************************** */

	initacc(&a_ng_m00_m02); initacc(&a_s_m00_m02); initacc(&a_q_m00_m02);
	initacc(&a_ha_m00_m02); initacc(&a_VA_m00_m02); initacc(&a_VD_m00_m02);

	initacc(&a_ng_m02_m04); initacc(&a_s_m02_m04); initacc(&a_q_m02_m04);
	initacc(&a_ha_m02_m04); initacc(&a_VA_m02_m04); initacc(&a_VD_m02_m04);

	initacc(&a_ng_m04_m06); initacc(&a_s_m04_m06); initacc(&a_q_m04_m06);
	initacc(&a_ha_m04_m06); initacc(&a_VA_m04_m06); initacc(&a_VD_m04_m06);

	initacc(&a_ng_m06_m08); initacc(&a_s_m06_m08); initacc(&a_q_m06_m08);
	initacc(&a_ha_m06_m08); initacc(&a_VA_m06_m08); initacc(&a_VD_m06_m08);

	initacc(&a_ng_m08_m10); initacc(&a_s_m08_m10); initacc(&a_q_m08_m10);
	initacc(&a_ha_m08_m10); initacc(&a_VA_m08_m10); initacc(&a_VD_m08_m10);

	initacc(&a_ng_m10_m12); initacc(&a_s_m10_m12); initacc(&a_q_m10_m12);
	initacc(&a_ha_m10_m12); initacc(&a_VA_m10_m12); initacc(&a_VD_m10_m12);

	initacc(&a_ng_m12_m14); initacc(&a_s_m12_m14); initacc(&a_q_m12_m14);
	initacc(&a_ha_m12_m14); initacc(&a_VA_m12_m14); initacc(&a_VD_m12_m14);

	initacc(&a_ng_m14_m16); initacc(&a_s_m14_m16); initacc(&a_q_m14_m16);
	initacc(&a_ha_m14_m16); initacc(&a_VA_m14_m16); initacc(&a_VD_m14_m16);

	initacc(&a_ng_m16_m18); initacc(&a_s_m16_m18); initacc(&a_q_m16_m18);
	initacc(&a_ha_m16_m18); initacc(&a_VA_m16_m18); initacc(&a_VD_m16_m18);

	initacc(&a_ng_m18_m20); initacc(&a_s_m18_m20); initacc(&a_q_m18_m20);
	initacc(&a_ha_m18_m20); initacc(&a_VA_m18_m20); initacc(&a_VD_m18_m20);

	initacc(&a_ng_m20); initacc(&a_s_m20); initacc(&a_q_m20);
	initacc(&a_ha_m20); initacc(&a_VA_m20); initacc(&a_VD_m20);
}


/* *************************************************************** */


settozero_gens_on_ha()
{
	initacc(&ha_ng_00_02); initacc(&ha_q_00_02);
	initacc(&ha_hs_00_02); initacc(&ha_VA_00_02); initacc(&ha_VD_00_02);

	initacc(&ha_ng_02_04); initacc(&ha_q_02_04);
	initacc(&ha_hs_02_04); initacc(&ha_VA_02_04); initacc(&ha_VD_02_04);

	initacc(&ha_ng_04_06); initacc(&ha_q_04_06);
	initacc(&ha_hs_04_06); initacc(&ha_VA_04_06); initacc(&ha_VD_04_06);

	initacc(&ha_ng_06_08); initacc(&ha_q_06_08);
	initacc(&ha_hs_06_08); initacc(&ha_VA_06_08); initacc(&ha_VD_06_08);

	initacc(&ha_ng_08_10); initacc(&ha_q_08_10);
	initacc(&ha_hs_08_10); initacc(&ha_VA_08_10); initacc(&ha_VD_08_10);
}


/* *************************************************************** */


settozero_gens_on_s()
{
	initacc(&s_ng_0); initacc(&s_q_0);
	initacc(&s_hs_0); initacc(&s_VA_0); initacc(&s_VD_0);

	initacc(&s_ng_0_106); initacc(&s_q_0_106);
	initacc(&s_hs_0_106); initacc(&s_VA_0_106); initacc(&s_VD_0_106);

	initacc(&s_ng_106_104); initacc(&s_q_106_104);
	initacc(&s_hs_106_104); initacc(&s_VA_106_104); initacc(&s_VD_106_104);

	initacc(&s_ng_104_102); initacc(&s_q_104_102);
	initacc(&s_hs_104_102); initacc(&s_VA_104_102); initacc(&s_VD_104_102);

	initacc(&s_ng_102_101); initacc(&s_q_102_101);
	initacc(&s_hs_102_101); initacc(&s_VA_102_101); initacc(&s_VD_102_101);

	initacc(&s_ng_01_02); initacc(&s_q_01_02);
	initacc(&s_hs_01_02); initacc(&s_VA_01_02); initacc(&s_VD_01_02);

	initacc(&s_ng_02_04); initacc(&s_q_02_04);
	initacc(&s_hs_02_04); initacc(&s_VA_02_04); initacc(&s_VD_02_04);

	initacc(&s_ng_04_06); initacc(&s_q_04_06);
	initacc(&s_hs_04_06); initacc(&s_VA_04_06); initacc(&s_VD_04_06);

	initacc(&s_ng_06_10); initacc(&s_q_06_10);
	initacc(&s_hs_06_10); initacc(&s_VA_06_10); initacc(&s_VD_06_10);

	initacc(&s_ng_10); initacc(&s_q_10);
	initacc(&s_hs_10); initacc(&s_VA_10); initacc(&s_VD_10);

}


/* ************************************************************** */


settozero_blocks_on_q()
{
    initacc(&Ave_ng_all); initacc(&Ave_s_all); initacc(&Ave_a_all); initacc(&Ave_hs_all); initacc(&Ave_ha_all); initacc(&Ave_VA_all); initacc(&Ave_VD_all);

    /* ************************************************************* */

    initacc(&Ave_q_ng_00_01); initacc(&Ave_q_s_00_01); initacc(&Ave_q_a_00_01); initacc(&Ave_q_hs_00_01); initacc(&Ave_q_ha_00_01); initacc(&Ave_q_VA_00_01); initacc(&Ave_q_VD_00_01);
	
    initacc(&Ave_q_ng_01_02); initacc(&Ave_q_s_01_02); initacc(&Ave_q_a_01_02); initacc(&Ave_q_hs_01_02); initacc(&Ave_q_ha_01_02); initacc(&Ave_q_VA_01_02); initacc(&Ave_q_VD_01_02);

    initacc(&Ave_q_ng_02_03); initacc(&Ave_q_s_02_03); initacc(&Ave_q_a_02_03); initacc(&Ave_q_hs_02_03); initacc(&Ave_q_ha_02_03); initacc(&Ave_q_VA_02_03); initacc(&Ave_q_VD_02_03);

    initacc(&Ave_q_ng_03_04); initacc(&Ave_q_s_03_04); initacc(&Ave_q_a_03_04); initacc(&Ave_q_hs_03_04); initacc(&Ave_q_ha_03_04); initacc(&Ave_q_VA_03_04); initacc(&Ave_q_VD_03_04);

    initacc(&Ave_q_ng_04_05); initacc(&Ave_q_s_04_05); initacc(&Ave_q_a_04_05); initacc(&Ave_q_hs_04_05); initacc(&Ave_q_ha_04_05); initacc(&Ave_q_VA_04_05); initacc(&Ave_q_VD_04_05);

    initacc(&Ave_q_ng_05_06); initacc(&Ave_q_s_05_06); initacc(&Ave_q_a_05_06); initacc(&Ave_q_hs_05_06); initacc(&Ave_q_ha_05_06); initacc(&Ave_q_VA_05_06); initacc(&Ave_q_VD_05_06);

    initacc(&Ave_q_ng_06_07); initacc(&Ave_q_s_06_07); initacc(&Ave_q_a_06_07); initacc(&Ave_q_hs_06_07); initacc(&Ave_q_ha_06_07); initacc(&Ave_q_VA_06_07); initacc(&Ave_q_VD_06_07);

    initacc(&Ave_q_ng_07_08); initacc(&Ave_q_s_07_08); initacc(&Ave_q_a_07_08); initacc(&Ave_q_hs_07_08); initacc(&Ave_q_ha_07_08); initacc(&Ave_q_VA_07_08); initacc(&Ave_q_VD_07_08);

    initacc(&Ave_q_ng_08_09); initacc(&Ave_q_s_08_09); initacc(&Ave_q_a_08_09); initacc(&Ave_q_hs_08_09); initacc(&Ave_q_ha_08_09); initacc(&Ave_q_VA_08_09); initacc(&Ave_q_VD_08_09);

    initacc(&Ave_q_ng_09_10); initacc(&Ave_q_s_09_10); initacc(&Ave_q_a_09_10); initacc(&Ave_q_hs_09_10); initacc(&Ave_q_ha_09_10); initacc(&Ave_q_VA_09_10); initacc(&Ave_q_VD_09_10);
}


/* *************************************************************** */


settozero_blocks_on_a()
{
	initacc(&Ave_a_ng_00_02); initacc(&Ave_a_s_00_02); initacc(&Ave_a_q_00_02);
	initacc(&Ave_a_ha_00_02); initacc(&Ave_a_VA_00_02); initacc(&Ave_a_VD_00_02);

	initacc(&Ave_a_ng_02_04); initacc(&Ave_a_s_02_04); initacc(&Ave_a_q_02_04);
	initacc(&Ave_a_ha_02_04); initacc(&Ave_a_VA_02_04); initacc(&Ave_a_VD_02_04);

	initacc(&Ave_a_ng_04_06); initacc(&Ave_a_s_04_06); initacc(&Ave_a_q_04_06);
	initacc(&Ave_a_ha_04_06); initacc(&Ave_a_VA_04_06); initacc(&Ave_a_VD_04_06);

	initacc(&Ave_a_ng_06_08); initacc(&Ave_a_s_06_08); initacc(&Ave_a_q_06_08);
	initacc(&Ave_a_ha_06_08); initacc(&Ave_a_VA_06_08); initacc(&Ave_a_VD_06_08);

	initacc(&Ave_a_ng_08_10); initacc(&Ave_a_s_08_10); initacc(&Ave_a_q_08_10);
	initacc(&Ave_a_ha_08_10); initacc(&Ave_a_VA_08_10); initacc(&Ave_a_VD_08_10);

	initacc(&Ave_a_ng_10_12); initacc(&Ave_a_s_10_12); initacc(&Ave_a_q_10_12);
	initacc(&Ave_a_ha_10_12); initacc(&Ave_a_VA_10_12); initacc(&Ave_a_VD_10_12);

	initacc(&Ave_a_ng_12_14); initacc(&Ave_a_s_12_14); initacc(&Ave_a_q_12_14);
	initacc(&Ave_a_ha_12_14); initacc(&Ave_a_VA_12_14); initacc(&Ave_a_VD_12_14);

	initacc(&Ave_a_ng_14_16); initacc(&Ave_a_s_14_16); initacc(&Ave_a_q_14_16);
	initacc(&Ave_a_ha_14_16); initacc(&Ave_a_VA_14_16); initacc(&Ave_a_VD_14_16);

	initacc(&Ave_a_ng_16_18); initacc(&Ave_a_s_16_18); initacc(&Ave_a_q_16_18);
	initacc(&Ave_a_ha_16_18); initacc(&Ave_a_VA_16_18); initacc(&Ave_a_VD_16_18);

	initacc(&Ave_a_ng_18_20); initacc(&Ave_a_s_18_20); initacc(&Ave_a_q_18_20);
	initacc(&Ave_a_ha_18_20); initacc(&Ave_a_VA_18_20); initacc(&Ave_a_VD_18_20);

	initacc(&Ave_a_ng_20); initacc(&Ave_a_s_20); initacc(&Ave_a_q_20);
	initacc(&Ave_a_ha_20); initacc(&Ave_a_VA_20); initacc(&Ave_a_VD_20);

	/* ****************************************************** */

	initacc(&Ave_a_ng_m00_m02); initacc(&Ave_a_s_m00_m02); initacc(&Ave_a_q_m00_m02);
initacc(&Ave_a_ha_m00_m02); initacc(&Ave_a_VA_m00_m02); initacc(&Ave_a_VD_m00_m02);

	initacc(&Ave_a_ng_m02_m04); initacc(&Ave_a_s_m02_m04); initacc(&Ave_a_q_m02_m04);
initacc(&Ave_a_ha_m02_m04); initacc(&Ave_a_VA_m02_m04); initacc(&Ave_a_VD_m02_m04);

	initacc(&Ave_a_ng_m04_m06); initacc(&Ave_a_s_m04_m06); initacc(&Ave_a_q_m04_m06);
initacc(&Ave_a_ha_m04_m06); initacc(&Ave_a_VA_m04_m06); initacc(&Ave_a_VD_m04_m06);

	initacc(&Ave_a_ng_m06_m08); initacc(&Ave_a_s_m06_m08); initacc(&Ave_a_q_m06_m08);
initacc(&Ave_a_ha_m06_m08); initacc(&Ave_a_VA_m06_m08); initacc(&Ave_a_VD_m06_m08);

	initacc(&Ave_a_ng_m08_m10); initacc(&Ave_a_s_m08_m10); initacc(&Ave_a_q_m08_m10);
initacc(&Ave_a_ha_m08_m10); initacc(&Ave_a_VA_m08_m10); initacc(&Ave_a_VD_m08_m10);

	initacc(&Ave_a_ng_m10_m12); initacc(&Ave_a_s_m10_m12); initacc(&Ave_a_q_m10_m12);
initacc(&Ave_a_ha_m10_m12); initacc(&Ave_a_VA_m10_m12); initacc(&Ave_a_VD_m10_m12);

	initacc(&Ave_a_ng_m12_m14); initacc(&Ave_a_s_m12_m14); initacc(&Ave_a_q_m12_m14);
initacc(&Ave_a_ha_m12_m14); initacc(&Ave_a_VA_m12_m14); initacc(&Ave_a_VD_m12_m14);

	initacc(&Ave_a_ng_m14_m16); initacc(&Ave_a_s_m14_m16); initacc(&Ave_a_q_m14_m16);
initacc(&Ave_a_ha_m14_m16); initacc(&Ave_a_VA_m14_m16); initacc(&Ave_a_VD_m14_m16);

	initacc(&Ave_a_ng_m16_m18); initacc(&Ave_a_s_m16_m18); initacc(&Ave_a_q_m16_m18);
initacc(&Ave_a_ha_m16_m18); initacc(&Ave_a_VA_m16_m18); initacc(&Ave_a_VD_m16_m18);

	initacc(&Ave_a_ng_m18_m20); initacc(&Ave_a_s_m18_m20); initacc(&Ave_a_q_m18_m20);
initacc(&Ave_a_ha_m18_m20); initacc(&Ave_a_VA_m18_m20); initacc(&Ave_a_VD_m18_m20);

	initacc(&Ave_a_ng_m20); initacc(&Ave_a_s_m20); initacc(&Ave_a_q_m20);
	initacc(&Ave_a_ha_m20); initacc(&Ave_a_VA_m20); initacc(&Ave_a_VD_m20);
}


/* **************************************************************** */


settozero_blocks_on_ha()
{
	initacc(&Ave_ha_ng_00_02); initacc(&Ave_ha_q_00_02);
	initacc(&Ave_ha_hs_00_02); initacc(&Ave_ha_VA_00_02); initacc(&Ave_ha_VD_00_02);

	initacc(&Ave_ha_ng_02_04); initacc(&Ave_ha_q_02_04);
	initacc(&Ave_ha_hs_02_04); initacc(&Ave_ha_VA_02_04); initacc(&Ave_ha_VD_02_04);

	initacc(&Ave_ha_ng_04_06); initacc(&Ave_ha_q_04_06);
	initacc(&Ave_ha_hs_04_06); initacc(&Ave_ha_VA_04_06); initacc(&Ave_ha_VD_04_06);

	initacc(&Ave_ha_ng_06_08); initacc(&Ave_ha_q_06_08);
	initacc(&Ave_ha_hs_06_08); initacc(&Ave_ha_VA_06_08); initacc(&Ave_ha_VD_06_08);

	initacc(&Ave_ha_ng_08_10); initacc(&Ave_ha_q_08_10);
	initacc(&Ave_ha_hs_08_10); initacc(&Ave_ha_VA_08_10); initacc(&Ave_ha_VD_08_10);
}


/* *************************************************************** */


settozero_blocks_on_s()
{
	initacc(&Ave_s_ng_0); initacc(&Ave_s_q_0);
	initacc(&Ave_s_hs_0); initacc(&Ave_s_VA_0); initacc(&Ave_s_VD_0);

	initacc(&Ave_s_ng_0_106); initacc(&Ave_s_q_0_106);
	initacc(&Ave_s_hs_0_106); initacc(&Ave_s_VA_0_106); initacc(&Ave_s_VD_0_106);

	initacc(&Ave_s_ng_106_104); initacc(&Ave_s_q_106_104);
	initacc(&Ave_s_hs_106_104); initacc(&Ave_s_VA_106_104); initacc(&Ave_s_VD_106_104);

	initacc(&Ave_s_ng_104_102); initacc(&Ave_s_q_104_102);
	initacc(&Ave_s_hs_104_102); initacc(&Ave_s_VA_104_102); initacc(&Ave_s_VD_104_102);

	initacc(&Ave_s_ng_102_101); initacc(&Ave_s_q_102_101);
	initacc(&Ave_s_hs_102_101); initacc(&Ave_s_VA_102_101); initacc(&Ave_s_VD_102_101);

	initacc(&Ave_s_ng_01_02); initacc(&Ave_s_q_01_02);
	initacc(&Ave_s_hs_01_02); initacc(&Ave_s_VA_01_02); initacc(&Ave_s_VD_01_02);

	initacc(&Ave_s_ng_02_04); initacc(&Ave_s_q_02_04);
	initacc(&Ave_s_hs_02_04); initacc(&Ave_s_VA_02_04); initacc(&Ave_s_VD_02_04);

	initacc(&Ave_s_ng_04_06); initacc(&Ave_s_q_04_06);
	initacc(&Ave_s_hs_04_06); initacc(&Ave_s_VA_04_06); initacc(&Ave_s_VD_04_06);

	initacc(&Ave_s_ng_06_10); initacc(&Ave_s_q_06_10);
	initacc(&Ave_s_hs_06_10); initacc(&Ave_s_VA_06_10); initacc(&Ave_s_VD_06_10);

	initacc(&Ave_s_ng_10); initacc(&Ave_s_q_10);
	initacc(&Ave_s_hs_10); initacc(&Ave_s_VA_10); initacc(&Ave_s_VD_10);
}


/* *************************************************************** */


void writelastgeneration (gm, s, a, hs, ha, initialgen)
int gm[][MMM][2], initialgen[][31];
double s[][31], a[][31], hs[][31], ha[][31];
{
    fprintf(fpop,"%d\n", NIND);

    for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
		fprintf(fpop,"%d  %d\n", gm[i][k][0], gm[i][k][1]);

    fprintf(fdat,"%42.40f  %42.40f\n", addedfixed_s, addedfixed_a);

    for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
		fprintf(fdat,"%42.40f  %42.40f  %f  %f  %d\n", s[k][l], a[k][l], hs[k][l], ha[k][l], initialgen[k][l]);
}


/* ************************************************************ */


/*generate_bivariate_gammas()*/

/* GTVR algorithm given in Schmeiser and Lal (1982).  Allowable values of
   the correlation are restricted to be positive and 
   <= min (beta1, beta2)/(sqrt(beta1*beta2))   */

double bvgam(double alpha1, double alpha2, double beta1, double beta2, double rho, double *g1, double *g2)
{
   double betay1, betay2, betay3;
   double y1, y2, y3;
   double x1, x2;

/* betay1 and betay2 should be larger than zero */

   betay3 = rho*(sqrt(beta1*beta2));
   betay1 = beta1 - betay3;
   betay2 = beta2 - betay3;
/*   printf("bvgam: betay3 %lf betay1 %lf betay2 %lf\n", betay3, betay1, betay2);
   printf("bvgam: alpha1 %lf alpha2 %lf beta1 %lf beta2 %lf rho %lf\n", alpha1,alpha2,beta1,beta2,rho);
*/
/*Check that rho is within bounds (see Figs 1 and 2 of Schmeiser and Lal (1982)
   */

   if ((rho > beta1/betay3) || (rho > beta2/betay3 )) 
   {
      printf("bvgam: Gabort condition 1\n");
      gabort ("rho is too large for this algorithm");
   }
   if (rho <= 0.0 )
   {
      printf("bvgam: Gabort condition 2\n");
      gabort ("This algorithm requires positive correlations.");
   }
   y1 = (double)gengam(1.0, (float)betay1);
   y2 = (double)gengam(1.0, (float)betay2);
   y3 = (double)gengam(1.0, (float)betay3);
/*   printf("bvgam: y1 %lf, y2 %lf y3 %lf\n", y1, y2, y3);*/

   *g1 = (y1 + y3)/alpha1;
   *g2 = (y2 + y3)/alpha2;
/*   accum(&g1acc, *g1);
   accum(&g2acc, *g2);
   covaccum(&c1acc, *g1, *g2);
   printf("betay1 %lf, betay2 %lf, betay3 %lf\n", betay1, betay2, betay3);
   printf("y1 %lf, y2 %lf, y3 %lf\n", y1, y2, y3);
*/
/*   printf("bvgam: g1 %lf g2 %lf\n", *g1, *g2);*/
}


/* ********************************************************** */

