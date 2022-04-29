/* extinctiongwasPFS.c (only gives the output GENFILE, the other outputs are commented) */
/* New factors: */
/*   - Mating_RCW() added */
/*   - Environmental source of mortality (A = 0.3) */
/*   - Natural catastrophes (beta distribution with mean P = 0.058) */
/*   - Quantitative trait (all sites affect fitness, regardless the value of a; i.e. whatproportion = 0 and Psi = 1 in naturalfBOT.c) */

/* ***************************************************** */

#include "libhdr"
#include "ranlib.h"
#define NN 2001  /* max number of NIND and gen */
#define NO 10001 /* max born offspring */
#define MM 2001  /* max number of NCRO */
#define maxmpt 5001
#define normalthreshold 30

int muts, mutsL, lastinmutantspoissontable, lastinmutantspoissontableL, lastinrecombinantpoissontable;
int NMAX, NIND[NN], NCRO, NLOCI, TOTLOCI, reps_alive_gen[1001], NSEGLOCNP;
int gen, tgen, generations, i, j, k, l, m, rep, replicates, born, dead;
int mutants_ocurred, RM[NN], initialgenNP[MM][31], initialgen[MM][31], p1, p2, compl_fail, parc_fail;
int ran_i, ran_k, ran_l, ran_h, countNoSS, NoSS_k[60001], NoSS_l[60001];
int gm[NN][MM][2], gmo[NO][MM][2], gmb[NN][MM][2], sm[NN][MM][2];
int NINDNP, gmNP[10001][MM][2], type, currenttype, CHT, CHTg, dom, C, scaling, ECrm;
int a, b, n, nn[NN], father[NN], mother[NN], FS, assign;
int family[NN], maxfec[NN], count, neutral, fec_via, selfing, maternalf,zygotes;
int tagPS[NN], tagFS[NN], tagHS[NN], failed[NN];
int mark[NO], mark_zyg[NN], mark_zygb[NN], mark_zygo[NO];
int fathero[NO], mothero[NO], tagPSo[NO], tagFSo[NO], tagHSo[NO];
int ss, SNP[MM][31], pps;
int envmort, catastrophe;

double genval_sf[NN], genval_sv[NN], pm_sf[NN], pm_sv[NN], pm_sfb[NN], pm_svb[NN];
double genval_at[NN], pm_at[NN], pm_atb[NN], VE;
double Lambda_s, Lambda_L, MAXFEC;
double L, w, AA, Aa, aa, q[MM][31], qNP[MM][31], LEQ_NP, LEQ_NPf, LEQ_NPv, TLEQ_NP;
double ave_s, addedfixed_sNP, addedfixed_sf, addedfixed_sv, addedfixed_at, alpha_s, beta_s, k_s, ave_hs, zs, za;
double mutantspoissontable[maxmpt], mutantspoissontableL[maxmpt], recombinantpoissontable[maxmpt];
double rho, beta_a, alpha_a, ave_a, ave_ha, k_a, Vs, P_a, ave_aL;

double sNP[MM][31], atNP[MM][31], hsNP[MM][31], hatNP[MM][31], s[MM][31], at[MM][31], hs[MM][31], hat[MM][31];
double H, PFSM, leqLv[NN], leqNLv[NN], leqLf[NN], leqNLf[NN], leqT[NN];
double smat[NN][NN], mat[NN][NN];
double sk2, ne, alpha, ne_alpha, fvk, fvk_alpha, relaxedsel, OPT;
double cum_pm_sf[NN], Wfo, Wvo;
double numPSborn, numFSborn, numHSborn;

double bvgam();
int fecunditysearch();
void disorder_NoSS();

struct acc gmean_sf_zygotes[NN], gmean_sv_zygotes[NN], gmean_a_zygotes[NN], DHWE[NN],Ho[NN], Hw[NN],Ho_del[NN], mean_pcat[NN];
struct acc mean_N[NN], mean_FS[NN], gmean_sf[NN], gmean_sv[NN], gmean_at[NN], gvar_at[NN], SK2[NN], cfail[NN], pfail[NN], fp[NN], Fp[NN];
struct acc LEQ_Lf[NN], LEQ_Lv[NN], LEQ_NLf[NN], LEQ_NLv[NN], LEQ_T[NN], TIME_EXT, P_PSBORN[NN], P_FSBORN[NN], P_HSBORN[NN];
struct acc NeVk[NN], Nevk_alpha[NN], Fvk[NN], Fvk_alpha[NN], PROD[NN];

/* ********************** Distribution of genes  *********************** */

struct acc q_ng_00[1001], q_ng_00_01[1001], q_ng_01_02[1001], q_ng_02_03[1001];
struct acc q_ng_03_04[1001], q_ng_04_05[1001], q_ng_05_06[1001], q_ng_06_07[1001];
struct acc q_ng_07_08[1001], q_ng_08_09[1001], q_ng_09_10[1001], q_ng_10[1001];

struct acc s_ng_0[1001], s_ng_0_106[1001], s_ng_106_104[1001], s_ng_104_102[1001], s_ng_102_101[1001];
struct acc s_ng_01_02[1001], s_ng_02_04[1001], s_ng_04_06[1001], s_ng_06_10[1001], s_ng_10[1001];

struct acc qs_ng_0_106[1001], qs_ng_106_104[1001], qs_ng_104_102[1001], qs_ng_102_101[1001];
struct acc qs_ng_01_02[1001], qs_ng_02_04[1001], qs_ng_04_06[1001], qs_ng_06_10[1001], qs_ng_10[1001];

struct acc hs_0_106[1001], hs_106_104[1001], hs_104_102[1001], hs_102_101[1001];
struct acc hs_01_02[1001], hs_02_04[1001], hs_04_06[1001], hs_06_10[1001], hs_10[1001];

struct acc s_ng_0_fix[1001], s_ng_0_106_fix[1001], s_ng_106_104_fix[1001], s_ng_104_102_fix[1001], s_ng_102_101_fix[1001];
struct acc s_ng_01_02_fix[1001], s_ng_02_04_fix[1001], s_ng_04_06_fix[1001], s_ng_06_10_fix[1001], s_ng_10_fix[1001];

struct acc s_ng_0_lost[1001], s_ng_0_106_lost[1001], s_ng_106_104_lost[1001], s_ng_104_102_lost[1001], s_ng_102_101_lost[1001];
struct acc s_ng_01_02_lost[1001], s_ng_02_04_lost[1001], s_ng_04_06_lost[1001], s_ng_06_10_lost[1001], s_ng_10_lost[1001];

struct acc hs_ng_00_02[1001], hs_ng_02_04[1001], hs_ng_04_06[1001], hs_ng_06_08[1001], hs_ng_08_10[1001];

struct acc hs_ng_00_02_fix[1001], hs_ng_02_04_fix[1001], hs_ng_04_06_fix[1001], hs_ng_06_08_fix[1001], hs_ng_08_10_fix[1001];

struct acc hs_ng_00_02_lost[1001], hs_ng_02_04_lost[1001], hs_ng_04_06_lost[1001], hs_ng_06_08_lost[1001], hs_ng_08_10_lost[1001];

FILE *fptr, *fgen, *frep, *fdat, *fpop, *fout, *fped, *fpedf, *fpedv, *fshqNP, *fshq0, *fshqLG, *fdata, *fmap, *fphe;

/* ***************************************************** */

main()
{
	fgen = fopen ("genfile.dat","w");
//	fout = fopen ("outfile.dat","w");
//	fped = fopen ("pedfile.dat","w");
//	fshq0 = fopen ("shqG0file.dat","w");
//	fshqLG = fopen ("shqGLfile.dat","w");
//	fshqNP = fopen ("shqNPfile.dat","w");
//	fdata = fopen ("data.ped","w");
//	fmap = fopen ("data.map","w");
//	fphe = fopen ("qt.phe","w");

	getinputs();
	recombination_masks();
	natural_population();

	if (tracelevel!=0) 	fptr = fopen ("dfilename.dat","w");
//	if (fec_via<2)
//	{
//		fpedf = fopen ("fec_pedfile.dat","w");
//		fpedv = fopen ("via_pedfile.dat","w");
//	}	

	for (rep=1; rep<=replicates; rep++)
	{
		// counter for pedfile and data-map
		n = 0;
		ss = 0;

		// counter for SNP
		for (k=0; k<NCRO; k++) 
		for (l=0; l<NLOCI; l++)
		SNP[k][l] = 0;

		frep = fopen ("repfile.dat","a");
		fprintf (frep,"replicate %d\n", rep);
		fclose(frep);

		if (tracelevel!=0) fprintf (fptr,"\n***** Replicate %d *****\n", rep);

		natural_data();

//		if (rep == 1)
//		{
//			for (k=0; k<NCRO; k++) 
//			for (l=0; l<NLOCI; l++)
//			fprintf(fshqNP,"k=%d  l=%d  sNP=%f  hNP=%f  atNP=%f  hatNP=%f  qNP=%f  inigen=%d\n", k, l, sNP[k][l], hsNP[k][l], atNP[k][l], hatNP[k][l], qNP[k][l], initialgen[k][l]);
//		}

		sample();
		NIND[0] = NMAX;

		for (gen=0; gen<=generations; gen++)
		{
			// possibility of change mating type with time 
			if ( ((type==0)&&(CHT==99)) || ((type==0)&&(gen<CHTg)) || ((CHT==0)&&(gen>=CHTg)))		currenttype = 0;
			else if ( ((type==1)&&(CHT==99)) || ((type==1)&&(gen<CHTg)) || ((CHT==1)&&(gen>=CHTg)))		currenttype = 1;
			else if ( ((type==2)&&(CHT==99)) || ((type==2)&&(gen<CHTg)) || ((CHT==2)&&(gen>=CHTg)))		currenttype = 2;
			else if ( ((type==3)&&(CHT==99)) || ((type==3)&&(gen<CHTg)) || ((CHT==3)&&(gen>=CHTg)))		currenttype = 3;
			else if ( ((type==4)&&(CHT==99)) || ((type==4)&&(gen<CHTg)) || ((CHT==4)&&(gen>=CHTg)))		currenttype = 4;
			else if ( ((type==5)&&(CHT==99)) || ((type==5)&&(gen<CHTg)) || ((CHT==5)&&(gen>=CHTg)))		currenttype = 5;

			if (catastrophe == 1)		naturalcatastrophe();

			if ((NIND[gen] == 0) || (NIND[gen] == 1))
			{
//				fprintf(fgen,"%d ", gen);
				accum(&TIME_EXT, gen);
				break;
			}

			if (NIND[gen] == 1)	break;

			if (tracelevel!=0) fprintf (fptr,"\n***** Generation %d *****\n\n", gen);
			if (tracelevel!=0) fprintf (fptr,"NIND = %d  MAXFEC = %f\n", NIND[gen], MAXFEC);

			poisson_tables();
			accum(&mean_N[gen], (double)NIND[gen]);

			mutation_neutral ();
			mutation_selected ();
//			if (tracelevel!=0) dumpoffspringaftermutation();
			neutral_genes ();
			selected_genes ();

//			fprintf(fshqLG,"gen=%d\n", gen);
//			for (k=0; k<NCRO; k++) 
//			for (l=0; l<NLOCI; l++)
//			if ((q[k][l] > 0.0) && (at[k][l] != 0.0))	fprintf(fshqLG,"k=%d  l=%d  SNP%d s=%f  h=%f  a=%f  ah=%f  q=%f\n", k, l, (k*30)+l, s[k][l], hs[k][l], at[k][l], hat[k][l], q[k][l]);

//			if (gen == 0)
//			{
//				for (k=0; k<NCRO; k++) 
//				for (l=0; l<NLOCI; l++)
//				if (q[k][l] > 0.0)	fprintf(fshq0,"k=%d  l=%d  SNP%d s=%f  h=%f  a=%f  ah=%f  q=%f  inigen=%d\n", k, l, (k*30)+l, s[k][l], hs[k][l], at[k][l], hat[k][l], q[k][l], initialgenNP[k][l]);
//			}

			genotypic_values();
			
			if (gen%tgen == 0) 	if (NIND[gen] != 0)
			{
				reps_alive_gen[gen/tgen] ++;
				distribution_qsh();
			}

			phenotypeB();

			/* check that there are parents available for fecunditysearch */
			if ( ((pps<2)&&(selfing==0)) || ((pps==0)&&(selfing==1)) )
			{
//				fprintf(fgen,"%d ", gen);
				accum(&TIME_EXT, gen);
				break;	
			}

			if (tracelevel!=0) 			dumpphenotypes();
			if (currenttype != 2)			disorder_parents();
			if ((currenttype == 1)&&(ECrm==1))	avoid_fullsib_parents_EC();
			if ((currenttype == 3)&&(gen!=0))	sort_parents_PFS();

			coancestry_matrix();
			dumpparents();

//			if (rep == 1)
//			{
//				if (zygotes==0)			plink_files();
//				if ((zygotes==1)&&(gen==0))	pedfile_zygt0();
//			}

			fecundities();

			if ((currenttype==0)||(currenttype==3))	mating_RC();
			else if (currenttype==1)		mating_EC();
			else if (currenttype==2)		mating_CM();
			else if (currenttype==4)		mating_RCpol();
			else					mating_RCW();

			calculations();

//			if (tracelevel!=0) dumpoffspring();
		}
  	}
	printout();
	distribution_out();

	writeseed();
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();
	getintandskip("NMAX (max 1000):",&NMAX,1,2000);
	getrealandskip("MAXFEC per couple:",&MAXFEC,1.0,100.0);
	getintandskip("Random (0), equal contributions (1), circular (2), PFS (3), polygamy (4), RCW (5): ",&type,0,5);
	getintandskip("Avoidance of full-sib matings (EC)(0 no, 1 yes): ",&ECrm,0,1);
	getintandskip("Selfing in case of polygamy (0 no, 1 yes): ",&selfing,0,1);
	getrealandskip("Proportion of partial full-sib mating (0: Random): ",&PFSM,0.0,1.0);
	getrealandskip("Length of genome in Morgans (99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("NCRO (min 1, max 2000):",&NCRO,1,2000);
	getintandskip("NLOCI (first is neutral)(min 2, max 30):",&NLOCI,2,30);
	TOTLOCI = NCRO * NLOCI;

	getrealandskip("Lambda_s :",&Lambda_s,0.0,(double)infinity);
	getrealandskip("Lambda_L (s=1,h=0.02):",&Lambda_L,0.0,(double)infinity);
	getrealandskip("Beta_s :",&beta_s,0.0,(double)infinity);
	getrealandskip("Average |s| :",&ave_s,0.0,1.0);
	alpha_s = beta_s / ave_s;
	getintandskip("dom (constant 0, variable 1):",&dom,0,1);
	getrealandskip("Average h:",&ave_hs,0.0,(double)infinity);
	k_s = alpha_s * (pow((2.0*ave_hs),((-1.0)/beta_s))-1.0);

	getrealandskip("Absolute effect of lethals for QT in heter. Normal (aL,aL/2)(haL=1) :",&ave_aL,0.0,(double)infinity);
	getrealandskip("Beta_a :",&beta_a,0.0,(double)infinity);
	getrealandskip("Average |a| :",&ave_a,0.0,(double)infinity);
	alpha_a = beta_a / ave_a;
	getrealandskip("Average h (mod 2):",&ave_ha,0.0,(double)infinity);
	k_a = alpha_a * (pow((2.0*ave_ha),((-1.0)/beta_a))-1.0);
	getrealandskip("rho (99: a=s):",&rho,-1.0,99.0);
	getrealandskip("Proportion of positives for the trait :",&P_a,0.0,1.0);
	getrealandskip("Stabilizing selection (Vs) :",&Vs,0.0,(double)infinity);
	getrealandskip("Optimal:",&OPT,0.0,(double)infinity);
	getrealandskip("VE :",&VE,0.0,(double)infinity);

	getrealandskip("Relaxation factor (0: no relaxation; 1: full relaxaton) :",&relaxedsel,0.0,1.0);
	getintandskip("Neutral model :",&neutral,0,1);
	getintandskip("Scaling (0: no; 1: yes) :",&scaling,0,1);
	getintandskip("fec_via (0: equal; 1: 1/3-2/3; 2: viability; 3: fecundity) :",&fec_via,0,3);
	getintandskip("Fecundity as a maternal component (0: no; 1: yes) :",&maternalf,0,1);
	getintandskip("Pedfile information (0: survivors; 1: zygotes) :",&zygotes,0,1);

	getintandskip("Environmental source of mortality (0: no; 1: yes) :",&envmort,0,1);
	getintandskip("Natural catastrophes (P=0.058) (0: no; 1: yes) :",&catastrophe,0,1);

	getintandskip("Number of generations :",&generations,1,1000);
	getintandskip("Change system to: RC (0), EC (1), CM (2), PFS (3), without change (99)) : ",&CHT,0,99);
	getintandskip("Change system at generation : ",&CHTg,0,1000);
	getintandskip("Number of replicates :",&replicates,1,infinity);

	tgen = 5;
}

/* **************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ***************************************************** */

natural_population ()
{
	int g0, g1, dinitialgen, dN;
	
	double dadds, dadda, ds, da, dhs, dha;

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("popfile","r");

	fscanf(fpop,"%d", &dN);
	NINDNP = dN;

	for (i=0; i<dN; i++)
	for (k=0; k<NCRO; k++)
	{
		fscanf(fpop,"%d%d", &g0, &g1);
		gmNP[i][k][0] = g0;
		gmNP[i][k][1] = g1;
	}

	fclose(fpop);

//	COMMENT for out
//	printf("\n");
//	for (i=0; i<NINDNP; i++)
//	{
//		if ((gmNP[i][0][0] & RM[5])==RM[5])	printf("1 ");
//		else						printf("0 ");
//		if ((gmNP[i][0][1] & RM[5])==RM[5])	printf("1 ");
//		else						printf("0 ");
//	}

	/* ***** take effects of genes ***** */

	fdat=fopen("datafile","r");

	fscanf(fdat,"%lf%lf", &dadds, &dadda);
	addedfixed_sNP = dadds;
	addedfixed_sNP = addedfixed_sNP + ((1.0 - addedfixed_sNP) * relaxedsel);

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%lf%lf%lf%lf%d", &ds, &da, &dhs, &dha, &dinitialgen);
		sNP[k][l] = ds;
		atNP[k][l] = da;
		if (sNP[k][l] == (-1))	sNP[k][l] = sNP[k][l];
		else			sNP[k][l] = sNP[k][l] - (sNP[k][l] * relaxedsel);
		hsNP[k][l] = dhs;
		hatNP[k][l] = dha;
		
		initialgenNP[k][l] = dinitialgen;
	}

	/* ***** estimate LEQ in the natural population ***** */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gmNP[i][k][0] & RM[l])==RM[l])&&((gmNP[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gmNP[i][k][0] & RM[l])!=RM[l])&&((gmNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));
		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

		if (fec_via == 0)
		{
			if (l%2==0)	LEQ_NPf += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
			else		LEQ_NPv += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
		}
		else if (fec_via == 1)
		{
			if (l%3==0)	LEQ_NPf += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
			else		LEQ_NPv += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
		}
		else if (fec_via == 2)
		{
			LEQ_NPv += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
		}
		else
		{
			LEQ_NPf += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
		}


		TLEQ_NP += -(sNP[k][l] * qNP[k][l]);
		if (qNP[k][l] > 0.0)	NSEGLOCNP ++;
	} 

//	COMMENT for out
//	printf("\n LEQ_NP = %f", LEQ_NP);
//	for (k=0; k<NCRO; k++)	for (l=1; l<NLOCI; l++)
//	printf("\n k=%d l=%d   sNP=%f  hNP=%f  qNP=%f  inigen=%d  aNP=%f  haNP=%f", k, l, sNP[k][l], hsNP[k][l], qNP[k][l], initialgenNP[k][l], atNP[k][l], hatNP[k][l]);

	fclose(fdat);
}

/* ***************************************************** */

natural_data()
{
	addedfixed_sf = sqrt(addedfixed_sNP);
	addedfixed_sv = sqrt(addedfixed_sNP);
	addedfixed_at = 0.0;

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		s[k][l] = sNP[k][l];
		at[k][l] = atNP[k][l];
		hs[k][l] = hsNP[k][l];
		hat[k][l] = hatNP[k][l];
		
		initialgen[k][l] = initialgenNP[k][l];
//COMMENT
//		if (tracelevel!=0)    fprintf(fptr,"\n k=%d l=%d s=%f  h=%f  inigen=%d", k, l, s[k][l], hs[k][l], initialgen[k][l]);
	}
}

/* ***************************************************** */

sample ()
{
	int g;
	
	for (i=0; i<NMAX; i++)
	{
		ran_i = (int)(uniform() * NINDNP);

		for (k=0; k<NCRO; k++)
		{
			g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;
		}
	}
	for (i=0; i<NMAX; i++)
	for (k=0; k<NCRO; k++)
	{
		gm[i][k][0]=gmNP[i][k][0];
		gm[i][k][1]=gmNP[i][k][1];
	}
}

/* ***************************************************** */

poisson_tables ()
{   
    /* SELECTED LOCI WITH POISSON (2NL) NEW MUTATIONS */

    if ( (exp(-2.0*(double)NIND[gen]*Lambda_s) != 0.0)&&(2.0*(double)NIND[gen]*Lambda_s < normalthreshold) )
    generatepoissontable(2.0*(double)NIND[gen]*Lambda_s, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

    if ( (exp(-2.0*(double)NIND[gen]*Lambda_L) != 0.0)&&(2.0*(double)NIND[gen]*Lambda_L < normalthreshold) )
    generatepoissontable(2.0*(double)NIND[gen]*Lambda_L, &lastinmutantspoissontableL, mutantspoissontableL, maxmpt-1);

  /* NUMERO DE RECOMBINACIONES POISSON CON MEDIA L*/

    if ( (exp(-L) != 0.0) && (L < normalthreshold) )
    generatepoissontable(L, &lastinrecombinantpoissontable, recombinantpoissontable, maxmpt-1);
}

/* ***************************************************** */

mutation_neutral ()
{
    /* NEUTRAL GENES: (POISSON) 2N(Lambda_a) NEW MUTATIONS (TWO DIRECTIONS) */

    muts = mutationnumber();

    if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants = %d\n", muts);

    for (m=0; m<muts; m++)
    {
		ran_i = (int)(uniform()*NIND[gen]);
		ran_k = (int)(uniform()*NCRO);
		ran_h = (int)(uniform()*2.0);
		if ( (gm[ran_i][ran_k][ran_h] & RM[0])==RM[0] )
		    gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] & (~RM[0]));
		else	gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[0]);
    }
}

/* ***************************************************** */

void mutation_selected ()
{
	/* SELECTED GENES: (POISSON) 2N(Lambda_s + Lambda_L) NEW MUTATIONS */

	muts = mutationnumber();
	mutsL = mutationnumberL();

	if ( (muts + mutsL) == 0 ) goto label;

	if (tracelevel!=0)    fprintf(fptr,"\nNew selected mutants  muts=%d  mutsL=%d\n", muts, mutsL);

	countNoSS = 0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		if (q[k][l]==0.0)
		{
			countNoSS += 1;
			NoSS_k[countNoSS-1] = k;
			NoSS_l[countNoSS-1] = l;
		}
	}
   
//	if (tracelevel!=0)    fprintf(fptr,"\n countNoSS=%d\n", countNoSS);

	mutants_ocurred = 0;

	if (countNoSS != 0)
	{
		disorder_NoSS (NoSS_k,NoSS_l);

		for (m=0; m<countNoSS; m++)
		{
			if (mutants_ocurred==(muts+mutsL))    goto label;

			ran_i = (int)(uniform()*NIND[gen]);
   			ran_h = (int)(uniform()*2.0);
    			gm[ran_i][NoSS_k[m]][ran_h]=(gm[ran_i][NoSS_k[m]][ran_h] | RM[NoSS_l[m]]);
		   	mutants_ocurred += 1;
			initialgen[NoSS_k[m]][NoSS_l[m]] = gen;

			/* ****** Lethal mutations ****** */

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
				hs[NoSS_k[m]][NoSS_l[m]] = 0.02;
			}
			else
			{
	     			/* ****** values of s, hs ****** */
				if(rho==0.0)
		      		{
		    			za = gengam (alpha_a, beta_a);
		    			zs = gengam (alpha_s, beta_s);
					if (zs > 1.0)	zs = 1.0;
				}
		      		else if(rho==99.0)
		       		{
					za = gengam (alpha_a, beta_a);
					zs = za;
					if (zs > 1.0)	zs=1.0;
				}
				else
				{
					bvgam(alpha_s, alpha_a, beta_s, beta_a, rho, &zs, &za);
				}

				if (zs > 1.0)		zs = 1.0;
				s[NoSS_k[m]][NoSS_l[m]] = (-zs) - ((-zs) * relaxedsel);

				if (dom==0)	hs[NoSS_k[m]][NoSS_l[m]] = ave_hs;
				else		hs[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
			}

	    		if (uniform()<P_a)	at[NoSS_k[m]][NoSS_l[m]] = za;
	    		else			at[NoSS_k[m]][NoSS_l[m]] = (-za);

			if(dom==0)		hat[NoSS_k[m]][NoSS_l[m]] = ave_ha;
			else
			{
			       	if (at[NoSS_k[m]][NoSS_l[m]] < 0.0)	hat[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(k_a*at[NoSS_k[m]][NoSS_l[m]]);
	    			else					hat[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(-k_a*at[NoSS_k[m]][NoSS_l[m]]);
			}
			if(rho==99.0)	hat[NoSS_k[m]][NoSS_l[m]] = hs[NoSS_k[m]][NoSS_l[m]];
	    		if(zs==1.0)	hat[NoSS_k[m]][NoSS_l[m]] = 1.0;

			if (tracelevel!=0)    fprintf(fptr,"alpha_s=%f  beta_s=%f\n", alpha_s, beta_s);
			if (tracelevel!=0)    fprintf(fptr,"s=%f  hs=%f  a=%f  ha=%f\n", s[NoSS_k[m]][NoSS_l[m]], hs[NoSS_k[m]][NoSS_l[m]], at[NoSS_k[m]][NoSS_l[m]], hat[NoSS_k[m]][NoSS_l[m]]);
			if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  k=%d  l=%d   ran_h=%d\n", ran_i, NoSS_k[m],  NoSS_l[m], ran_h);
		}
	}

	for (m=mutants_ocurred; m<muts; m++)
	{
		ran_i = (int)(uniform()*NIND[gen]);
		ran_k = (int)(uniform()*NCRO);
		do {ran_l = (int)(uniform()*NLOCI);}   while (ran_l==0);
		ran_h = (int)(uniform()*2.0);

//		fprintf(frep,"recurrent mutation\n");
		
		gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[ran_l]);

//		if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  ran_k=%d  ran_l=%d  ran_h=%d\n", ran_i, ran_k, ran_l, ran_h);
	}
	label: /* end of mutations */;
}

/* ***************************************************** */

int mutationnumber()
{
	int r;
	if ((2.0*(double)NIND[gen]*Lambda_s < normalthreshold) && (exp(-2.0*(double)NIND[gen]*Lambda_s) != 0.0) )
	{
		r = poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal(2.0*(double)NIND[gen]*Lambda_s, sqrt(2.0*(double)NIND[gen]*Lambda_s)) );
	return(r);
}

/* ***************************************************** */

int mutationnumberL()
{
	int r;
	if ((2.0*(double)NIND[gen]*Lambda_L < normalthreshold) && (exp(-2.0*(double)NIND[gen]*Lambda_L) != 0.0) )
	{
		r = poisson(lastinmutantspoissontableL, mutantspoissontableL);
	}
	else r = (int)( normal(2.0*(double)NIND[gen]*Lambda_L, sqrt(2.0*(double)NIND[gen]*Lambda_L)) );
	return(r);
}

/* ***************************************************** */

/*from naturalvBOT.c*/
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
   	printf("bvgam: alpha1 %lf alpha2 %lf beta1 %lf beta2 %lf rho %lf\n", alpha1,alpha2,beta1,beta2,rho);*/
	/*Check that rho is within bounds (see Figs 1 and 2 of Schmeiser and Lal (1982)*/

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

/* ***************************************************** */

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

/* ***************************************************** */

dumpoffspringaftermutation()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nOffspring after mutation (gm0 gm1)\n");	
	for (i=0; i<NIND[gen]; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}

/* ***************************************************** */

neutral_genes ()
{
	H = 0.0; 
	double dhwe=0.0, Hobs=0.0;

	for (k=0; k<NCRO; k++)
	{
	    AA=0.0; Aa=0.0; aa=0.0; 

	    for (i=0; i<NIND[gen]; i++)
	    {
	        if (((gm[i][k][0] & RM[0])==RM[0])&&((gm[i][k][1] & RM[0])==RM[0]))	aa+=1.0;
	        else    if (((gm[i][k][0] & RM[0])!=RM[0])&&((gm[i][k][1] & RM[0])!=RM[0])) AA+=1.0;
	        else		Aa+=1.0;
	    }

	    q[k][0] = (aa/(double)NIND[gen])+(Aa/(2.0*(double)NIND[gen]));
	    H += 2.0 * q[k][0] * (1.0 - q[k][0]);
	    Hobs += (Aa/(double)NIND[gen]);
	    if (q[k][0] > 0.0)   SNP[k][0] = 1;
	    if (q[k][0] > 0.0) 	 dhwe += 1.0 - ((Aa/(double)NIND[gen]) / (2.0 * q[k][0] * (1.0 - q[k][0])));

//	    if (tracelevel!=0)    fprintf(fptr,"%d\t0\t%1.0f\t%1.0f\t%1.0f\t%f\t%d\n",k,AA,Aa,aa,q[k][0],initialgen[k][0]);	
	}
	accum(&Hw[gen], H/(double)NCRO);
	accum(&Ho[gen], Hobs/(double)NCRO);
	accum(&DHWE[gen], (dhwe/(double)NCRO));
}

/* ***************************************************** */

selected_genes()
{
	double Hobs=0.0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND[gen]; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND[gen])+(Aa/(2.0*(double)NIND[gen]));
		if (q[k][l] > 0.0)   SNP[k][l] = 1;

	    	Hobs += (Aa/(double)NIND[gen]);
	}

	accum(&Ho_del[gen], Hobs/(double)(NCRO*(NLOCI-1)));
}

/* ***************************************************** */

genotypic_values ()
{
//COMMENT
//	if (tracelevel!=0)	fprintf(fptr,"\n Genotypic values \n");

	for (i=0; i<NIND[gen]; i++)
	{
		leqLf[i] = 0.0;
		leqNLf[i] = 0.0;
		leqLv[i] = 0.0;
		leqNLv[i] = 0.0;
		leqT[i] = 0.0;

		genval_sf[i] = addedfixed_sf;
		genval_sv[i] = addedfixed_sv;
		genval_at[i] = 0.0;

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if (initialgen[k][l] != (-99))
		{
//			if (tracelevel!=0)
//			{
//				if ( ((gm[i][k][0] & RM[l])==RM[l]) || ((gm[i][k][1] & RM[l])==RM[l]) )	fprintf(fptr,"gen=%d i=%d k=%d l=%d s=%f h=%f q=%f ini=%d\n", gen, i, k, l, s[k][l], hs[k][l], q[k][l], initialgen[k][l]);
//			}

			// Trait
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_at[i] += at[k][l];
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else	genval_at[i] += (at[k][l]*hat[k][l]);

			if (fec_via == 0)
			{

				if (l%2==0)
				{
					// fecundity
			    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sf[i] *= (1.0 + s[k][l]);
					else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
					else	genval_sf[i] *= (1.0 + (s[k][l]*hs[k][l]));
				}
				else
				{
					// viability
			    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sv[i] *= (1.0 + s[k][l]);
					else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
					else	genval_sv[i] *= (1.0 + (s[k][l]*hs[k][l]));
				}
			}
			else if (fec_via == 1)
			{
				if (l%3==0)
				{
					// fecundity
			    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sf[i] *= (1.0 + s[k][l]);
					else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
					else	genval_sf[i] *= (1.0 + (s[k][l]*hs[k][l]));
				}
				else
				{
					// viability
			    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sv[i] *= (1.0 + s[k][l]);
					else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
					else	genval_sv[i] *= (1.0 + (s[k][l]*hs[k][l]));
				}
			}
			else if (fec_via == 2)
			{
				genval_sf[i] *= 1.0;

				// viability
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sv[i] *= (1.0 + s[k][l]);
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
				else	genval_sv[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
			else
			{
				genval_sv[i] *= 1.0;

				// fecundity
			    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sf[i] *= (1.0 + s[k][l]);
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
				else	genval_sf[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}

			if (s[k][l]==(-1.0))
			{
				leqT[i] += -(s[k][l] * q[k][l]);
				if (fec_via == 0)
				{
					if (l%2==0)	leqLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
					else		leqLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
				else if (fec_via == 1)
				{
					if (l%3==0)	leqLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
					else		leqLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
				else if (fec_via == 2)
				{
					leqLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
				else
				{
					leqLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
			}
			else 
			{
				leqT[i] += -(s[k][l] * q[k][l]);
				if (fec_via == 0)
				{
					if (l%2==0)	leqNLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
					else		leqNLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
				else if (fec_via == 1)
				{
					if (l%3==0)	leqNLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
					else		leqNLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
				else if (fec_via == 2)
				{
					leqNLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
				else
				{
					leqNLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				}
			}
		}

//COMMENT
//		if (tracelevel!=0)	fprintf(fptr," %d    genval_sf = %f    genval_sv = %f\n", i, genval_sf[i], genval_sv[i]);
	}
}

/* ***************************************************** */

distribution_qsh()
{
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	if (initialgen[k][l] != (-99)) 
	{
		if (q[k][l] == 0.0)					accum (&q_ng_00[gen/tgen],  1.0);
		else if ( (q[k][l] > 0.0) && (q[k][l] < 0.1) )	accum (&q_ng_00_01[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.1) && (q[k][l] < 0.2) )	accum (&q_ng_01_02[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.2) && (q[k][l] < 0.3) )	accum (&q_ng_02_03[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.3) && (q[k][l] < 0.4) )	accum (&q_ng_03_04[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.4) && (q[k][l] < 0.5) )	accum (&q_ng_04_05[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.5) && (q[k][l] < 0.6) )	accum (&q_ng_05_06[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.6) && (q[k][l] < 0.7) )	accum (&q_ng_06_07[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.7) && (q[k][l] < 0.8) )	accum (&q_ng_07_08[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.8) && (q[k][l] < 0.9) )	accum (&q_ng_08_09[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.9) && (q[k][l] < 1.0) )	accum (&q_ng_09_10[gen/tgen],  1.0);
		else if (q[k][l] == 1.0)	accum (&q_ng_10[gen/tgen],  1.0);

		if (q[k][l] == 0.0)
		{
			if (s[k][l] == 0.0)						accum (&s_ng_0_lost[gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10_lost[gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))					accum (&s_ng_10_lost[gen/tgen],  1.0);
		}
		else if (q[k][l] == 1.0)
		{
			if (s[k][l] == 0.0)						accum (&s_ng_0_fix[gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10_fix[gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))					accum (&s_ng_10_fix[gen/tgen],  1.0);
		}
		else
		{
			if (s[k][l] == 0.0)						accum (&s_ng_0[gen/tgen],  1.0);	
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		{accum (&s_ng_0_106[gen/tgen],  1.0);   accum (&qs_ng_0_106[gen/tgen],  q[k][l]);   accum (&hs_0_106[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	{accum (&s_ng_106_104[gen/tgen],  1.0); accum (&qs_ng_106_104[gen/tgen],  q[k][l]); accum (&hs_106_104[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	{accum (&s_ng_104_102[gen/tgen],  1.0); accum (&qs_ng_104_102[gen/tgen],  q[k][l]); accum (&hs_104_102[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		{accum (&s_ng_102_101[gen/tgen],  1.0); accum (&qs_ng_102_101[gen/tgen],  q[k][l]); accum (&hs_102_101[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		{accum (&s_ng_01_02[gen/tgen],  1.0);   accum (&qs_ng_01_02[gen/tgen],  q[k][l]);   accum (&hs_01_02[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		{accum (&s_ng_02_04[gen/tgen],  1.0);   accum (&qs_ng_02_04[gen/tgen],  q[k][l]);   accum (&hs_02_04[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		{accum (&s_ng_04_06[gen/tgen],  1.0);   accum (&qs_ng_04_06[gen/tgen],  q[k][l]);   accum (&hs_04_06[gen/tgen],  hs[k][l]);}
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		{accum (&s_ng_06_10[gen/tgen],  1.0);   accum (&qs_ng_06_10[gen/tgen],  q[k][l]);   accum (&hs_06_10[gen/tgen],  hs[k][l]);}
			else if (s[k][l] == (-1.0))					{accum (&s_ng_10[gen/tgen],  1.0);      accum (&qs_ng_10[gen/tgen],  q[k][l]);      accum (&hs_10[gen/tgen],  hs[k][l]);}
//			if (tracelevel!=0)	fprintf(fptr,"k=%d l=%d s=%f hs=%f q=%f ini=%d\n", k, l, s[k][l], hs[k][l], q[k][l], initialgen[k][l]);
		}

		if (q[k][l] == 0.0)
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10_lost[gen/tgen],  1.0);
		}
		else if (q[k][l] == 1.0)
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10_fix[gen/tgen],  1.0);
		}
		else
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10[gen/tgen],  1.0);
		}
//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"k=%d l=%d s=%f hs=%f q=%f ini=%d\n", k, l, s[k][l], hs[k][l], q[k][l], initialgen[k][l]);
	}
}

/* **************************************************************** */

phenotypeB()
{
	int ii, it;
	double gsum_sf=0.0, gsum_sv=0.0, gsum_at=0.0, gsum2_at=0.0, sum_leqNLf=0.0, sum_leqLf=0.0, sum_leqNLv=0.0, sum_leqLv=0.0, sum_leqT=0.0;
	double sq_pm_a;

	for (i=0; i<NIND[gen]; i++)
	{
		pm_at[i] = genval_at[i] + normal(0.0, sqrt(VE));
	}

	// scaling to generation 0 fitness

	if (scaling == 1)	if (gen == 0)
	{
		Wfo=0.0; Wvo=0.0;
		for (i=0; i<NIND[gen]; i++)
		{
			if(Vs==0.0)
			{
				Wfo += genval_sf[i]/(double)NIND[gen];
				Wvo += genval_sv[i]/(double)NIND[gen];
			}
			else
			{
				sq_pm_a = pow ((pm_at[i]-OPT), 2.0);
				Wfo += genval_sf[i]* exp(-0.5*(sq_pm_a)/Vs)/(double)NIND[gen];
				Wvo += genval_sv[i]* exp(-0.5*(sq_pm_a)/Vs)/(double)NIND[gen];
			}
		}
	}

	for (i=0; i<NIND[gen]; i++)
	{
		if(Vs==0.0) 
		{
			if (scaling == 1)
			{
				/* environmental source of mortality A = 0.3 */
				if (envmort==1)
				{
					pm_sf[i] = (genval_sf[i] / Wfo) - 0.3;
					pm_sv[i] = (genval_sv[i] / Wvo) - 0.3;
					if (tracelevel!=0)   fprintf(fptr,"wfo=%f wvo=%f i%d genval_sf=%f (scaled %f)  genval_sv=%f (scaled %f)\n", Wfo, Wvo, i, genval_sf[i],(genval_sf[i] / Wfo), genval_sv[i],(genval_sv[i] / Wvo));
				}
				else
				{
					pm_sf[i] = genval_sf[i] / Wfo;
					pm_sv[i] = genval_sv[i] / Wvo;
					if (tracelevel!=0)   fprintf(fptr,"wfo=%f wvo=%f i%d genval_sf=%f (scaled %f)  genval_sv=%f (scaled %f)\n", Wfo, Wvo, i, genval_sf[i],(genval_sf[i] / Wfo), genval_sv[i],(genval_sv[i] / Wvo));
				}

				if (pm_sf[i] < 0.0)	pm_sf[i] = 0.0;
				if (pm_sf[i] > 1.0)	pm_sf[i] = 1.0;
				if (pm_sv[i] < 0.0)	pm_sv[i] = 0.0;
				if (pm_sv[i] > 1.0)	pm_sv[i] = 1.0;
			}
			else
			{
				pm_sf[i] = genval_sf[i];
				pm_sv[i] = genval_sv[i];
			}
		}
		else
		{
			sq_pm_a = pow ((pm_at[i]-OPT), 2.0);
			if (scaling == 1)
			{
				/* environmental source of mortality A = 0.3 */
				if (envmort==1)
				{
					pm_sf[i] = (genval_sf[i] * exp(-0.5*(sq_pm_a)/Vs)/ Wfo) - 0.3;
					pm_sv[i] = (genval_sv[i] * exp(-0.5*(sq_pm_a)/Vs)/ Wvo) - 0.3;
				}
				else
				{
					pm_sf[i] = genval_sf[i] * exp(-0.5*(sq_pm_a)/Vs)/ Wfo;
					pm_sv[i] = genval_sv[i] * exp(-0.5*(sq_pm_a)/Vs)/ Wvo;
				}

				if (pm_sf[i] < 0.0)	pm_sf[i] = 0.0;
				if (pm_sf[i] > 1.0)	pm_sf[i] = 1.0;
				if (pm_sv[i] < 0.0)	pm_sv[i] = 0.0;
				if (pm_sv[i] > 1.0)	pm_sv[i] = 1.0;
			}
			else
			{
				pm_sf[i] = genval_sf[i] * exp(-0.5*(sq_pm_a)/Vs);
				pm_sv[i] = genval_sv[i] * exp(-0.5*(sq_pm_a)/Vs);
			}
		}

		sum_leqNLf += leqNLf[i];
		sum_leqLf += leqLf[i];
		sum_leqNLv += leqNLv[i];
		sum_leqLv += leqLv[i];
		sum_leqT += leqT[i];

		gsum_sf += pm_sf[i];
		gsum_sv += pm_sv[i];
		gsum_at += pm_at[i];
		gsum2_at += (pm_at[i]*pm_at[i]);
	}

	accum (&LEQ_NLf[gen], sum_leqNLf/(double)NIND[gen]);
	accum (&LEQ_Lf[gen], sum_leqLf/(double)NIND[gen]);
	accum (&LEQ_NLv[gen], sum_leqNLv/(double)NIND[gen]);
	accum (&LEQ_Lv[gen], sum_leqLv/(double)NIND[gen]);
	accum (&LEQ_T[gen], sum_leqT/(double)NIND[gen]);

	accum (&gmean_sf[gen], gsum_sf/(double)NIND[gen]);
	accum (&gmean_sv[gen], gsum_sv/(double)NIND[gen]);
	accum (&gmean_at[gen], gsum_at/(double)NIND[gen]);
	accum (&gvar_at[gen], (gsum2_at - (gsum_at*gsum_at / (double)NIND[gen])) / ((double)NIND[gen] - 1.0));

	if (tracelevel!=0)   fprintf(fptr,"\nLEQ_NLf = %f  LEQ_Lf = %f  LEQ_NLv = %f  LEQ_Lv = %f  LEQ_T = %f  mean fec = %f  mean via = %f  mean a = %f\n",
				sum_leqNLf/(double)NIND[gen], sum_leqLf/(double)NIND[gen], 
				sum_leqNLv/(double)NIND[gen], sum_leqLv/(double)NIND[gen], sum_leqT/(double)NIND[gen],
				gsum_sf/(double)NIND[gen], gsum_sv/(double)NIND[gen], gsum_at/(double)NIND[gen]);

	// Neutral case
	if (neutral == 1)
	for (i=0; i<NIND[gen]; i++)
	{
		pm_sf[i] = 1.0;
		pm_sv[i] = 1.0;
	}

	/* check that there are parents available for fecunditysearch */
	pps = NIND[gen];
	if ((currenttype==4)||(currenttype==5))
	{
		pps = 0;
		for (i=0; i<NIND[gen]; i++)
		{
			if (pm_sf[i]>0) pps += 1;
		}
	}
}

/* ***************************************************** */

naturalcatastrophe()
{
	int alive[NN];
	double pcat;
	int j = 0;
	dead = 0;

	/* probability from a beta distribution */
	pcat = genbet (0.5, 8);

	for (i=0; i<NIND[gen]; i++)	alive[i] = 1;
	if (tracelevel!=0)	fprintf (fptr,"\nNatural catastrophes\n");
	if (tracelevel!=0)	for (i=0; i<NIND[gen]; i++)	fprintf (fptr,"i %d  alive %d  gm[i][0][0]=%d\n", i, alive[i],gm[i][0][0]);

	for (i=0; i<NIND[gen]; i++)	if (uniform() < pcat)	{ alive[i] = 0; dead += 1;}
	if (tracelevel!=0)	fprintf (fptr,"\nAfter catastrophe (P = %f  dead = %d  N = %d)\n", pcat, dead, NIND[gen]-dead);
	if (tracelevel!=0)	for (i=0; i<NIND[gen]; i++)	fprintf (fptr,"i %d  alive %d  gm[i][0][0]=%d\n", i, alive[i],gm[i][0][0]);

	for (i=0; i<(NIND[gen]-dead); i++)
	{
		if (alive[i] == 0)
		{
			ran_i = NIND[gen] - dead + j;
			if (alive[ran_i] == 0)	do {j += 1; ran_i += 1;}	while((alive[ran_i] == 0)&&(ran_i<=NIND[gen]));
	
			for (k=0; k<NCRO; k++)
			{
				a=gm[i][k][0]; gm[i][k][0]=gm[ran_i][k][0]; gm[ran_i][k][0]=a;
				a=gm[i][k][1]; gm[i][k][1]=gm[ran_i][k][1]; gm[ran_i][k][1]=a;
			}

			a=father[i]; father[i]=father[ran_i]; father[ran_i]=a;
			a=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=a;
			if(zygotes==1)	a=mark_zyg[i]; mark_zyg[i]=mark_zyg[ran_i]; mark_zyg[ran_i]=a;
			a=alive[i]; alive[i]=alive[ran_i]; alive[ran_i]=a;

			j += 1;
		}
	}

	accum (&mean_pcat[gen], pcat);
	NIND[gen] = NIND[gen] - dead;

	if (tracelevel!=0)	fprintf (fptr,"\nAfter ordering\n");
	if (tracelevel!=0)	for (i=0; i<(NIND[gen]+dead); i++)	fprintf (fptr,"i %d  alive %d  gm[i][0][0]=%d\n", i, alive[i],gm[i][0][0]);
}

/* ***************************************************** */

dumpphenotypes()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nFecundity, viability and trait phenotype of parents\n");
	for (i=0; i<NIND[gen]; i++)		fprintf(fptr,"%d  %f %f %f\n", i, pm_sf[i], pm_sv[i], pm_at[i]);
}

/* ***************************************************** */

disorder_parents()
{
	/***** disorder parents before mating (except for circular mating - type=2) *****/

	for (i=0; i<NIND[gen]; i++)
	{
		ran_i=(int)(uniform()*NIND[gen]);
		for (k=0; k<NCRO; k++)
		{
			a=gm[i][k][0]; gm[i][k][0]=gm[ran_i][k][0]; gm[ran_i][k][0]=a;
			a=gm[i][k][1]; gm[i][k][1]=gm[ran_i][k][1]; gm[ran_i][k][1]=a;
		}
		w=pm_sf[i]; pm_sf[i]=pm_sf[ran_i]; pm_sf[ran_i]=w;
		w=pm_sv[i]; pm_sv[i]=pm_sv[ran_i]; pm_sv[ran_i]=w;
		w=pm_at[i]; pm_at[i]=pm_at[ran_i]; pm_at[ran_i]=w;

		a=father[i]; father[i]=father[ran_i]; father[ran_i]=a;
		a=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=a;
		if(zygotes==1)	a=mark_zyg[i]; mark_zyg[i]=mark_zyg[ran_i]; mark_zyg[ran_i]=a;
	}
}

/* ***************************************************** */

avoid_fullsib_parents_EC()
{
	/***** Avoid full-sib couples with EC (type=1) *****/

	for (j=1; j<=3; j++)
	for (i=1; i<NIND[gen]; i+=2)
	{
		if (father[i] == father[i-1])
		{
			ran_i=(int)(uniform()*NIND[gen]);
			for (k=0; k<NCRO; k++)
			{
				a=gm[i][k][0]; gm[i][k][0]=gm[ran_i][k][0]; gm[ran_i][k][0]=a;
				a=gm[i][k][1]; gm[i][k][1]=gm[ran_i][k][1]; gm[ran_i][k][1]=a;
			}
			w=pm_sf[i]; pm_sf[i]=pm_sf[ran_i]; pm_sf[ran_i]=w;
			w=pm_sv[i]; pm_sv[i]=pm_sv[ran_i]; pm_sv[ran_i]=w;
			w=pm_at[i]; pm_at[i]=pm_at[ran_i]; pm_at[ran_i]=w;

			a=father[i]; father[i]=father[ran_i]; father[ran_i]=a;
			a=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=a;
			if(zygotes==1)	a=mark_zyg[i]; mark_zyg[i]=mark_zyg[ran_i]; mark_zyg[ran_i]=a;
		}
	}
}

/* ***************************************************** */

sort_parents_PFS()
{
	int nb[NN], brot=0, tcouplebrot=0, realFS, mark_brot[NN], bro1, bro2;	
	for (i=0; i<NIND[gen]; i++)	{mark_brot[i]=0; mark[i]=0;}

	/**** couples of brothers ****/ 

	for (i=0; i<NIND[gen]; i++)
	for (j=0; j<NIND[gen]; j++)
	{
		if ((i!=j)&&(father[i]==father[j])&&(mother[i]==mother[j]))
		{
			if ((mark_brot[i]==0)&&(mark_brot[j]==0))
			{
				brot+=1;
				mark_brot[i]=brot;
				mark_brot[j]=brot;
			}
			else if ((mark_brot[i]!=0)&&(mark_brot[j]==0))	mark_brot[j]=mark_brot[i];
		}
	}

	for (b=1; b<=brot; b++)
	{
		nb[b]=0;
		for (i=0; i<NIND[gen]; i++)	if(mark_brot[i]==b)	nb[b] += 1;
		tcouplebrot += nb[b]*0.5;
	}

	/**** real FS ****/

	FS = (int)((PFSM * NIND[gen])/2);
	if (tcouplebrot<FS)	realFS=tcouplebrot;
	else			realFS=FS;

	accum(&mean_FS[gen], realFS);
	if (tracelevel!=0)	fprintf(fptr,"FS = %d couples (%f*%dindividuals)\nnº of brother couples = %d\nrealFS=%d\n\n",FS,PFSM,NIND[gen],tcouplebrot,realFS);

	/**** sort brothers ****/

	for (b=0; b<(realFS*2); b+=2)
	{
		for (i=0; i<NIND[gen]; i++)	if ((mark_brot[i]!=0)&&(mark[i]==0)&&((nb[mark_brot[i]]-2)>=0))	{bro1 = i; mark[i]=1; break;}
		for (i=0; i<NIND[gen]; i++)	if ((mark_brot[i]==mark_brot[bro1])&&(mark[i]==0))		{bro2 = i; mark[i]=1; break;}
		nb[mark_brot[bro1]] = nb[mark_brot[bro1]]-2;

		if (tracelevel!=0)	fprintf(fptr,"bro1=%d fa=%d  mo=%d mark_brot=%d		bro2=%d fa=%d  mo=%d mark_brot=%d\n", bro1, father[bro1], mother[bro1], mark_brot[bro1], bro2, father[bro2], mother[bro2], mark_brot[bro2]);
			
		for (k=0; k<NCRO; k++)
		{
		   	a = gm[b][k][0]; gm[b][k][0] = gm[bro1][k][0]; gm[bro1][k][0] = a;
		  	a = gm[b][k][1]; gm[b][k][1] = gm[bro1][k][1]; gm[bro1][k][1] = a;

		 	a = gm[b+1][k][0]; gm[b+1][k][0] = gm[bro2][k][0]; gm[bro2][k][0] = a;
		 	a = gm[b+1][k][1]; gm[b+1][k][1] = gm[bro2][k][1]; gm[bro2][k][1] = a;
		}

		a = father[b]; father[b] = father[bro1]; father[bro1] = a;
		a = mother[b]; mother[b] = mother[bro1]; mother[bro1] = a;
		a = tagPS[b]; tagPS[b] = tagPS[bro1]; tagPS[bro1] = a;
		a = tagFS[b]; tagFS[b] = tagFS[bro1]; tagFS[bro1] = a;
		a = tagHS[b]; tagHS[b] = tagHS[bro1]; tagHS[bro1] = a;
		a = mark[b]; mark[b] = mark[bro1]; mark[bro1] = a;
		a = mark_brot[b]; mark_brot[b] = mark_brot[bro1]; mark_brot[bro1] = a;
		if (zygotes==1)	a = mark_zyg[b]; mark_zyg[b] = mark_zyg[bro1]; mark_zyg[bro1] = a;

		a = father[b+1]; father[b+1] = father[bro2]; father[bro2] = a;
		a = mother[b+1]; mother[b+1] = mother[bro2]; mother[bro2] = a;
		a = tagPS[b+1]; tagPS[b+1] = tagPS[bro2]; tagPS[bro2] = a;
		a = tagFS[b+1]; tagFS[b+1] = tagFS[bro2]; tagFS[bro2] = a;
		a = tagHS[b+1]; tagHS[b+1] = tagHS[bro2]; tagHS[bro2] = a;
		a = mark[b+1]; mark[b+1] = mark[bro2]; mark[bro2] = a;
		a = mark_brot[b+1]; mark_brot[b+1] = mark_brot[bro2]; mark_brot[bro2] = a;
		if (zygotes==1)	a = mark_zyg[b+1]; mark_zyg[b+1] = mark_zyg[bro2]; mark_zyg[bro2] = a;
	}

	if (tracelevel!=0)	fprintf (fptr,"\n");
	if (tracelevel!=0)	for (i=0; i<NIND[gen]; i++)	fprintf (fptr,"Randomised born FS i %d  mark_brot %d  fa %d  mo %d\n", i, mark_brot[i], father[i], mother[i]);
}

/* ***************************************************** */

coancestry_matrix()
{
	if (gen==0)
	{
		for (i=0; i<NIND[gen]; i++)
		for (j=0; j<NIND[gen]; j++)
		{
			if (i == j)	mat[i][j] = 0.5;
			else		mat[i][j] = 0.0;
		}
	}
	else
	{
		for (i=0; i<NIND[gen]; i++)
		for (j=0; j<NIND[gen]; j++)
		{
			if (i == j)	mat[i][j] = 0.5 * (1.0 + smat[father[i]][mother[i]]);
			else		mat[i][j] = 0.25 * (smat[father[i]][father[j]] + smat[father[i]][mother[j]] + smat[mother[i]][father[j]] + smat[mother[i]][mother[j]]);
		}
	}

	for (i=0; i<NIND[gen]; i++)
	for (j=0; j<NIND[gen]; j++)
	{
		smat[i][j] = mat[i][j];
		accum(&fp[gen], mat[i][j]); 
		accum(&Fp[gen], (2.0*mat[i][i]-1.0)); 
	}

//	if (tracelevel!=0)    
//	for (i=0; i<NIND[gen]; i++)
//	{
//		fprintf(fptr, "gen=%d   i=%d   fat=%d   mot=%d   Fp=%f\n", gen, i, father[i], mother[i], 2.0*mat[i][i]-1.0);
//	}

}

/* ***************************************************** */

dumpparents()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nParents\n");
	for (i=0; i<NIND[gen]; i++)   fprintf(fptr,"i=%d  father=%d  mother=%d  pm_sf=%f pm_sv=%f pm_at=%f  Fp=%f\n", i, father[i], mother[i], pm_sf[i], pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0);
}

/* ***************************************************** */

plink_files()
{
	if (gen==0)
	{
		fprintf(fped, "ind,fath,moth,fit,QT,Fp,gen\n");
		if (fec_via<2)
		{
			fprintf(fpedf, "ind,fath,moth,fit,QT,Fp,gen\n");
			fprintf(fpedv, "ind,fath,moth,fit,QT,Fp,gen\n");
		}

		// data.map

/*		ss = 0;
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			ss ++;
			if (SNP[k][l] == 1)	fprintf(fmap,"1 SNP%d 0 %d\n", ss-1, (k*30)+l);
		}
*/
	}

	for (i=0; i<NIND[gen]; i++)
	{
		n ++;

		// pedfile

		if (gen == 0)		fprintf(fped, "%d,%d,%d,%f,%f,%f,%d\n", n, 0, 0, pm_sf[i]*pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
		else if (gen == 1)	fprintf(fped, "%d,%d,%d,%f,%f,%f,%d\n", n, father[i]+1, mother[i]+1, pm_sf[i]*pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
		else			fprintf(fped, "%d,%d,%d,%f,%f,%f,%d\n", n, father[i]+1+nn[gen-2], mother[i]+1+nn[gen-2], pm_sf[i]*pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
		nn[gen] = n;

		// pedfile - fecundity - viability

		if (fec_via<2)
		{
			if (gen == 0)
			{
				fprintf(fpedf, "%d,%d,%d,%f,%f,%f,%d\n", n, 0, 0, pm_sf[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
				fprintf(fpedv, "%d,%d,%d,%f,%f,%f,%d\n", n, 0, 0, pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
			}
			else if (gen == 1)
			{
				fprintf(fpedf, "%d,%d,%d,%f,%f,%f,%d\n",n,father[i]+1, mother[i]+1, pm_sf[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
				fprintf(fpedv, "%d,%d,%d,%f,%f,%f,%d\n",n,father[i]+1, mother[i]+1, pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
			}
			else
			{
				fprintf(fpedf, "%d,%d,%d,%f,%f,%f,%d\n",n,father[i]+1+nn[gen-2], mother[i]+1+nn[gen-2], pm_sf[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
				fprintf(fpedv, "%d,%d,%d,%f,%f,%f,%d\n",n,father[i]+1+nn[gen-2], mother[i]+1+nn[gen-2], pm_sv[i], pm_at[i], 2.0*mat[i][i]-1.0, gen);
			}
			nn[gen] = n;
		}

		// qt.phe

//		fprintf(fphe,"1 IND%d %f\n", n, pm_at[i] + normal(0.0, sqrt(VE)));

		// data.ped
/*
		fprintf(fdata,"1 IND%d 0 0 1 0   ", n);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (SNP[k][l] == 1)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(fdata,"2 2  ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(fdata,"1 1  ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(fdata,"2 1  ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(fdata,"1 2  ");
		}
		fprintf(fdata,"\n");
*/
	}

	return(0);
}

/* ***************************************************** */

pedfile_zygt0()
{
	fprintf(fped, "ind,fath,moth,fit,QT,gen\n");
	if (fec_via<2)
	{
		fprintf(fpedf, "ind,fath,moth,fit,QT,gen\n");
		fprintf(fpedv, "ind,fath,moth,fit,QT,gen\n");
	}
	
	for (i=0; i<NMAX; i++)
	{
		n ++;
		mark_zyg[i] = n;
		fprintf(fped, "%d,%d,%d,%f,%f,%d\n", n, 0, 0, pm_sf[i]*pm_sv[i], pm_at[i], gen);
			
		if (fec_via<2)
		{
			fprintf(fpedf, "%d,%d,%d,%f,%f,%d\n", n, 0, 0, pm_sf[i], pm_at[i], gen);
			fprintf(fpedv, "%d,%d,%d,%f,%f,%d\n", n, 0, 0, pm_sv[i], pm_at[i], gen);
		}
	}
}

/* ***************************************************** */

pedfile_zyg()
{
	n ++;

	fprintf(fped, "%d,%d,%d,%f,%f,%d\n", n, mark_zyg[p1], mark_zyg[p2], pm_sfb[j]*pm_svb[j], pm_atb[j], gen+1);
	if (fec_via<2)
	{
		fprintf(fpedf, "%d,%d,%d,%f,%f,%d\n", n, mark_zyg[p1], mark_zyg[p2], pm_sfb[j], pm_atb[j], gen+1);
		fprintf(fpedv, "%d,%d,%d,%f,%f,%d\n", n, mark_zyg[p1], mark_zyg[p2], pm_svb[j], pm_atb[j], gen+1);
	}

	mark_zygb[j] = n;
}

/* ***************************************************** */

fecundities()
{
	double cum=0.0;

	/********* keep genotypes of parents ***********/ 

	for (i=0; i<NIND[gen]; i++)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];
	}

	/********* Accumulate psf fecundity values ***********/

	for (i=0; i<NIND[gen]; i++)
	{
		cum += pm_sf[i];
		cum_pm_sf[i] = cum;
	}
	for (i=0; i<NIND[gen]; i++)	cum_pm_sf[i] = cum_pm_sf[i] / cum;

//	COMMENT
	if (tracelevel!=0)	for (i=0; i<NIND[gen]; i++)	fprintf (fptr," %d    cum_pm_sf = %f\n", i, cum_pm_sf[i]);
}

/* ***************************************************** */

int fecunditysearch(double *array, int size)
{
	int i, p;
	double r;

	r = uniform();
	if ( (r >= 0.0) && (r <= array[0]) )	p = 0; 
	else for (i=1; i<size; i++)
		if ( (r > array[i-1]) && (r <= array[i]) )	p = i; 

	return(p);
}

/* ***************************************************** */

mating_RC()
{
	int rnd, a;
	double K;
	for (i=0; i<NIND[gen]; i++)	family[i] = 0;
	for (i=0; i<NMAX; i++)	failed[i] = 1;
	
	/******************************************/

	if (tracelevel!=0)	fprintf(fptr,"\nMatings MAXFEC %f\n", MAXFEC);

	/******************* Round of random contributions mating ***********************/

	if (tracelevel!=0)	fprintf (fptr,"\nRound of random contributions mating\n\n");

	born = 0;
	dead = 0;
	if (NIND[gen]%2 != 0)	NIND[gen] = NIND[gen] - 1;

	for (i=0; i<NIND[gen]; i+=2)
	{
		p1 = i;
		p2 = p1 + 1;

		if (maternalf==0)	K = MAXFEC * sqrt(pm_sf[p1] * pm_sf[p2]);
		else			K = MAXFEC * pm_sf[p1];

		maxfec[i] = (int)(K*2);
//		maxfec[i] = (int)(K*2 + 0.5);
//		if (tracelevel!=0)	fprintf(fptr,"\nK  %f    macfec %d\n", K, maxfec[i]);

		for (j=0; j<maxfec[i]; j++)
		{
			if (tracelevel!=0)
			{
				if (zygotes==0)	fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d    p2 = %d\n", i, j, maxfec[i], p1, p2);
				else		fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d (zyg%d)    p2 = %d (zyg%d)\n", i, j, maxfec[i], p1, mark_zyg[p1], p2, mark_zyg[p2]);
			}
			gametes_offspring();
//			if ((zygotes==1)&&(rep == 1)) 	pedfile_zyg();
			bornoffspring();
		}
	}

	// First assignation of progeny to population

	if (tracelevel!=0)	fprintf (fptr,"First assignation of progeny to population\n\n");

	for (b=1; b<=born; b++)	mark[b] = 0;
	assign = 0;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	{
		if(zygotes==0)	fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d\n", b, fathero[b], mothero[b], mark[b]);
		else		fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d  zyg_n %d\n", b, fathero[b], mothero[b], mark[b], mark_zygo[b]);
	}
	if (tracelevel!=0)	fprintf (fptr,"\n");

	// Randomise born progeny

	for (b=1; b<=born; b++)
	{
		rnd = (int)(uniform()*born) + 1;

		for (k=0; k<NCRO; k++)
		{
		   	a = gmo[b][k][0]; gmo[b][k][0] = gmo[rnd][k][0]; gmo[rnd][k][0] = a;
		   	a = gmo[b][k][1]; gmo[b][k][1] = gmo[rnd][k][1]; gmo[rnd][k][1] = a;
		}
	   	a = fathero[b]; fathero[b] = fathero[rnd]; fathero[rnd] = a;
	   	a = mothero[b]; mothero[b] = mothero[rnd]; mothero[rnd] = a;
	   	a = tagPSo[b]; tagPSo[b] = tagPSo[rnd]; tagPSo[rnd] = a;
	   	a = tagFSo[b]; tagFSo[b] = tagFSo[rnd]; tagFSo[rnd] = a;
	   	a = tagHSo[b]; tagHSo[b] = tagHSo[rnd]; tagHSo[rnd] = a;
	   	a = mark[b]; mark[b] = mark[rnd]; mark[rnd] = a;
	   	if(zygotes==1)	a = mark_zygo[b]; mark_zygo[b] = mark_zygo[rnd]; mark_zygo[rnd] = a;
	}
	if (tracelevel!=0)	for (b=1; b<=born; b++)
	{
		if(zygotes==0)	fprintf (fptr,"Randomised born b %d  mark_b %d  fa %d  mo %d\n", b, mark[b], fathero[b], mothero[b]);
		else 		fprintf (fptr,"Randomised born b %d  mark_b %d  zyg_n %d  fa %d  mo %d\n", b, mark[b], mark_zygo[b], fathero[b], mothero[b]);
	}

	// Assignation

	if (assign == NMAX) goto endofmating_RC;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"\nAvailable born before final assignment b %d  fa %d  mo %d  mark %d", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n\n");

	for (i=0; i<NMAX; i++)  if (failed[i] == 1)
	{
		if (assign == NMAX) goto endofmating_RC;

		for (b=1; b<=born; b++)  if (mark[b] == 0)
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
			failed[i] = 0;

			if (tracelevel!=0)
			{
				if(zygotes==0)	fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
				else		fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  zyg_n %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], mark_zygo[b], fathero[b], mothero[b], failed[i]);
			}

			break;
		}
	}

	compl_fail=0;
	parc_fail=0;
	for (i=0; i<NIND[gen]; i+=2)
	{
		if (family[i] == 0)		compl_fail += 1;
		else if (family[i] == 1)	parc_fail += 1;
	}
	accum (&cfail[gen], (double)compl_fail);
	accum (&pfail[gen], (double)parc_fail);
	if (tracelevel!=0)	fprintf (fptr,"\nFailed families (0): %d\nParcial failed families (1): %d\n", compl_fail, parc_fail);

	/******************* End of rounds of mating ***********************/

	endofmating_RC: /******/;

	NIND[gen+1] = assign;
}

/* ***************************************************** */

mating_EC()
{
	double K;
	for (i=0; i<NIND[gen]; i++)	family[i] = 0;
	for (i=0; i<NMAX; i++)		failed[i] = 1;
	
	/******************************************/

	if (tracelevel!=0)	fprintf(fptr,"\nMatings MAXFEC %f\n", MAXFEC);

	/******************* Round of equal contributions mating ***********************/

	if (tracelevel!=0)	fprintf (fptr,"\nRound of equal contributions mating\n\n");

	born = 0;
	dead = 0;

	if (NIND[gen]%2 != 0)	NIND[gen] = NIND[gen] - 1;

	for (i=0; i<NIND[gen]; i+=2)
	{
		p1 = i;
		p2 = p1 + 1;

		if (maternalf==0)	K = MAXFEC * sqrt(pm_sf[p1] * pm_sf[p2]);
		else			K = MAXFEC * pm_sf[p1];

		maxfec[i] = (int)(K*2);
//		maxfec[i] = (int)(K*2 + 0.5);
//		if (tracelevel!=0)	fprintf(fptr,"\nK  %f    macfec %d\n", K, maxfec[i]);

		for (j=0; j<maxfec[i]; j++)
		{
			if (tracelevel!=0)
			{
				if (zygotes==0)	fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d    p2 = %d\n", i, j, maxfec[i], p1, p2);
				else		fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d (zyg%d)    p2 = %d (zyg%d)\n", i, j, maxfec[i], p1, mark_zyg[p1], p2, mark_zyg[p2]);
			}
			gametes_offspring();
//			if ((zygotes==1)&&(rep == 1)) 	pedfile_zyg();
			bornoffspring();
		}
	}

	// First assignation of progeny to population

	if (tracelevel!=0)	fprintf (fptr,"\nFirst assignation of progeny to population\n\n");

	for (b=1; b<=born; b++)	mark[b] = 0;
	assign = 0;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	{
		if(zygotes==0)	fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d\n", b, fathero[b], mothero[b], mark[b]);
		else		fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d  zyg_n %d\n", b, fathero[b], mothero[b], mark[b], mark_zygo[b]);
	}
	if (tracelevel!=0)	fprintf (fptr,"\n");

	for (i=0; i<NMAX; i++)
	for (b=1; b<=born; b++)
	{
		if (assign == NMAX) goto endofmating_EC;

		if (i%2 == 0)		C = (i - (NIND[gen] * (i / NIND[gen])));
		else			C = (i - 1 - (NIND[gen] * (i / NIND[gen])));

		if ( (mark[b] == 0) && (fathero[b] == C) )
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
			failed[i] = 0;
			break;
		}
	}

	// Final assignation of progeny to population

	if (assign == NMAX) goto endofmating_EC;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"\nAvailable born before final assignment b %d  fa %d  mo %d  mark %d", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n");

	for (i=0; i<NMAX; i++)  if (failed[i] == 1)
	{
		if (assign == NMAX) goto endofmating_EC;

		for (b=1; b<=born; b++)  if (mark[b] == 0)
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
			failed[i] = 0;

			if (tracelevel!=0)
			{
				if(zygotes==0)	fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
				else		fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  zyg_n %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], mark_zygo[b], fathero[b], mothero[b], failed[i]);
			}

			break;
		}
	}

	compl_fail=0;
	parc_fail=0;
	for (i=0; i<NIND[gen]; i+=2)
	{
		if (family[i] == 0)		compl_fail += 1;
		else if (family[i] == 1)	parc_fail += 1;
	}
	accum (&cfail[gen], (double)compl_fail);
	accum (&pfail[gen], (double)parc_fail);

	/******************* End of rounds of mating ***********************/

	endofmating_EC: /******/;

	// Offspring renumbering

	if (tracelevel!=0)	fprintf (fptr,"\nProgeny before renumbering\n\n");
	if (tracelevel!=0)	for (i=0; i<NMAX; i++)	if (failed[i] == 0)
	fprintf (fptr,"i=%d  fa=%d  mo=%d  tPS=%d  tFS=%d  tHS=%d\n", i, father[i], mother[i], tagPS[i], tagFS[i], tagHS[i]);

	count = 0;
	for (i=0; i<NMAX; i++)  if (failed[i] == 0)
	{
		count ++;

		for (k=0; k<NCRO; k++)
		{
		   	gm[count-1][k][0] = gm[i][k][0];
		   	gm[count-1][k][1] = gm[i][k][1];
		}
		father[count-1] = father[i];
		mother[count-1] = mother[i];
		tagPS[count-1] = tagPS[i];
		tagFS[count-1] = tagFS[i];
		tagHS[count-1] = tagHS[i];
		if(zygotes==1)	mark_zyg[count-1] = mark_zyg[i];
	}

	if (tracelevel!=0)	fprintf (fptr,"\nProgeny after renumbering\n\n");
	if (tracelevel!=0)	for (i=0; i<count; i++)
	fprintf (fptr,"i=%d  fa=%d  mo=%d  tPS=%d  tFS=%d  tHS=%d\n", i, father[i], mother[i], tagPS[i], tagFS[i], tagHS[i]);

	if (tracelevel!=0)
	{
		fprintf (fptr,"\nNMAX = %d  born = %d  assign = %d\n\n", NMAX, born, assign);
		for (i=0; i<assign; i++)
		fprintf (fptr,"i %d  fa %d  mo %d  tagPS %d  tagFS %d  tagHS %d\n", i, father[i], mother[i], tagPS[i], tagFS[i], tagHS[i]);
		fprintf (fptr,"\nFailed families (0): %d\nParcial failed families (1): %d\n", compl_fail, parc_fail);
	}

	NIND[gen+1] = assign;
}

/* ***************************************************** */

mating_CM()
{
	double K;
	for (i=0; i<NIND[gen]; i++)	family[i] = 0;
	for (i=0; i<NMAX; i++)	failed[i] = 1;
	
	/******************************************/

	if (tracelevel!=0)	fprintf(fptr,"\nMatings MAXFEC %f\n", MAXFEC);

	/******************* Round of circular mating ***********************/

	if (tracelevel!=0)	fprintf (fptr,"\nRound of circular mating\n\n");

	born = 0;
	dead = 0;

	if (NIND[gen]%2 != 0)	NIND[gen] = NIND[gen] - 1;

	for (i=0; i<NIND[gen]; i++)
	{
		p1 = i;
		if (p1 == NIND[gen]-1)	p2 = 0;
		else p2 = p1 + 1;

		if (maternalf==0)	K = MAXFEC * sqrt(pm_sf[p1] * pm_sf[p2]);
		else			K = MAXFEC * pm_sf[p1];

		maxfec[i] = (int)(K);
//		maxfec[i] = (int)(K + 0.5);
		if (tracelevel!=0)	fprintf (fptr, "pm_sf=%f pm_sf=%f K=%f maxfec=%d\n", pm_sf[p1], pm_sf[p2], K, maxfec[i]);

		for (j=0; j<maxfec[i]; j++)
		{
			if (tracelevel!=0)
			{
				if (zygotes==0)	fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d    p2 = %d\n", i, j, maxfec[i], p1, p2);
				else		fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d (zyg%d)    p2 = %d (zyg%d)\n", i, j, maxfec[i], p1, mark_zyg[p1], p2, mark_zyg[p2]);
			}
			gametes_offspring();
//			if ((zygotes==1)&&(rep == 1)) 	pedfile_zyg();
			bornoffspring();
		}
	}

	// First assignation of progeny to population

	if (tracelevel!=0)	fprintf (fptr,"\nFirst assignation of progeny to population\n\n");

	for (b=1; b<=born; b++)	mark[b] = 0;
	assign = 0;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	{
		if(zygotes==0)	fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d\n", b, fathero[b], mothero[b], mark[b]);
		else		fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d  zyg_n %d\n", b, fathero[b], mothero[b], mark[b], mark_zygo[b]);
	}
	if (tracelevel!=0)	fprintf (fptr,"\n");

	for (i=0; i<NMAX; i++)
	for (b=1; b<=born; b++)
	{
		if (assign == NMAX) goto endofmating_CM;

		C = (i - (NIND[gen] * (i / NIND[gen])));

		if ( (mark[b] == 0) && (fathero[b] == C) )
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
			failed[i] = 0;
			if (tracelevel!=0)	fprintf (fptr,"After first assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
			break;
		}
	}

	// Second assignation of progeny to population

	if (tracelevel!=0)	fprintf (fptr,"\nSecond assignation of progeny to population\n\n");

	if (assign == NMAX) goto endofmating_CM;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"Available born before second assignment b %d  fa %d  mo %d  mark %d\n", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n");

	for (i=0; i<NMAX; i++)  if (failed[i] == 1)
	{
		if (assign == NMAX) goto endofmating_CM;

		C = (i - (NIND[gen] * (i / NIND[gen])));

		for (b=1; b<=born; b++)  if (mark[b] == 0)
		{
			if (C == NIND[gen]-1)
			{
				if ( (fathero[b]==C-1) || (fathero[b]==1) )
				{
					assign ++;
					for (k=0; k<NCRO; k++)
					{
					   	gm[i][k][0] = gmo[b][k][0];
					   	gm[i][k][1] = gmo[b][k][1];
					}
					father[i] = fathero[b];
					mother[i] = mothero[b];
					tagPS[i] = tagPSo[b];
					tagFS[i] = tagFSo[b];
					tagHS[i] = tagHSo[b];
					family[father[i]] += 1;
					family[mother[i]] += 1;
					mark[b] = 1;
					if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
					failed[i] = 0;
					if (tracelevel!=0)	fprintf (fptr,"After second assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
					break;
				}
			}
			else
			{
				if ( (fathero[b]==C-1) || (fathero[b]==C+1) )
				{
					assign ++;
					for (k=0; k<NCRO; k++)
					{
					   	gm[i][k][0] = gmo[b][k][0];
					   	gm[i][k][1] = gmo[b][k][1];
					}
					father[i] = fathero[b];
					mother[i] = mothero[b];
					tagPS[i] = tagPSo[b];
					tagFS[i] = tagFSo[b];
					tagHS[i] = tagHSo[b];
					family[father[i]] += 1;
					family[mother[i]] += 1;
					mark[b] = 1;
					if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
					failed[i] = 0;
					if (tracelevel!=0)	fprintf (fptr,"After second assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
					break;
				}
			}
		}
	}

	// Final assignation of progeny to population

	if (assign == NMAX) goto endofmating_CM;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"\nAvailable born before final assignment b %d  fa %d  mo %d  mark %d", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n");

	for (i=0; i<NMAX; i++)  if (failed[i] == 1)
	{
		if (assign == NMAX) goto endofmating_CM;

		for (b=1; b<=born; b++)  if (mark[b] == 0)
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
			failed[i] = 0;

			if (tracelevel!=0)
			{
				if(zygotes==0)	fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
				else		fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  zyg_n %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], mark_zygo[b], fathero[b], mothero[b], failed[i]);
			}

			break;
		}
	}

	compl_fail=0;
	parc_fail=0;
	for (i=0; i<NIND[gen]; i+=2)
	{
		if (family[i] == 0)		compl_fail += 1;
		else if (family[i] == 1)	parc_fail += 1;
	}
	accum (&cfail[gen], (double)compl_fail);
	accum (&pfail[gen], (double)parc_fail);
	
	/******************* End of rounds of mating ***********************/

	endofmating_CM: /******/;

	// Offspring renumbering

	if (tracelevel!=0)	fprintf (fptr,"\nProgeny before renumbering\n\n");
	if (tracelevel!=0)	for (i=0; i<NMAX; i++)	if (failed[i] == 0)
	fprintf (fptr,"i=%d  fa=%d  mo=%d  tPS=%d  tFS=%d  tHS=%d\n", i, father[i], mother[i], tagPS[i], tagFS[i], tagHS[i]);

	count = 0;
	for (i=0; i<NMAX; i++)  if (failed[i] == 0)
	{
		count ++;

		for (k=0; k<NCRO; k++)
		{
		   	gm[count-1][k][0] = gm[i][k][0];
		   	gm[count-1][k][1] = gm[i][k][1];
		}
		father[count-1] = father[i];
		mother[count-1] = mother[i];
		tagPS[count-1] = tagPS[i];
		tagFS[count-1] = tagFS[i];
		tagHS[count-1] = tagHS[i];
		if(zygotes==1)	mark_zyg[count-1] = mark_zyg[i];
	}

	if (tracelevel!=0)	fprintf (fptr,"\nProgeny after renumbering\n\n");
	if (tracelevel!=0)	for (i=0; i<count; i++)
	fprintf (fptr,"i=%d  fa=%d  mo=%d  tPS=%d  tFS=%d  tHS=%d\n", i, father[i], mother[i], tagPS[i], tagFS[i], tagHS[i]);

	if (tracelevel!=0)
	{
		fprintf (fptr,"\nNMAX = %d  born = %d  assign = %d\n\n", NMAX, born, assign);
		for (i=0; i<assign; i++)
		fprintf (fptr,"i %d  fa %d  mo %d  tagPS %d  tagFS %d  tagHS %d\n", i, father[i], mother[i], tagPS[i], tagFS[i], tagHS[i]);
		fprintf (fptr,"\nFailed families (0): %d\nParcial failed families (1): %d\n", compl_fail, parc_fail);
	}

	NIND[gen+1] = assign;
}

/* ***************************************************** */

mating_RCpol()
{
	int rnd, a, tries;
	for (i=0; i<NIND[gen]; i++)	family[i] = 0;
	
	/******************* Round of random contributions mating ***********************/

	if (tracelevel!=0)	fprintf (fptr,"\nRound of random contributions mating\n\n");
	born = 0;
	dead = 0;
	tries = 0;
	j = 0;

	for (i=0; i<NMAX; i++)
	{
		generahijo: /* ***** */;

		p1 = fecunditysearch(cum_pm_sf, NIND[gen]);
		if (selfing==0) 
		{
			if (maternalf==0)	do { p2 = fecunditysearch(cum_pm_sf, NIND[gen]); }	while (p2==p1);
			else			do { p2 = (int)(uniform()*NIND[gen]); }			while (p2==p1);
		}
		else
		{
			if (maternalf==0)	p2 = fecunditysearch(cum_pm_sf, NIND[gen]);
			else 			p2 = (int)(uniform()*NIND[gen]);
		}

		if (tracelevel!=0)
		{
			if (zygotes==0)	fprintf(fptr,"Progeny %d  p1 = %d    p2 = %d\n", i, p1, p2);
			else		fprintf(fptr,"Progeny %d  p1 = %d (zyg%d)   p2 = %d (zyg%d)\n", i, p1, mark_zyg[p1], p2, mark_zyg[p2]);
		}

		gametes_offspring();
//		if ((zygotes==1)&&(rep == 1))	pedfile_zyg();
		
		/********* Born or dead individuals *********/

		if (uniform() <= pm_svb[j])
		{
			born ++;
			for (k=0; k<NCRO; k++)
			{
			   	gmo[i][k][0] = gmb[j][k][0];
			   	gmo[i][k][1] = gmb[j][k][1];
			}
			fathero[i] = p1;
			mothero[i] = p2;
			if(zygotes==1)	mark_zygo[i] = mark_zygb[j];

			if (p1 == p2)	tagPSo[born] = 1;
			if (gen > 0)
			{	
				if ( ((father[p1] == father[p2]) && (mother[p1] == mother[p2])) ||
					((father[p1] == mother[p2]) && (mother[p1] == father[p2])) )	tagFSo[born] = 1;
				else if ((father[p1] == father[p2]) || (father[p1] == mother[p2]) ||
					(mother[p1] == father[p2]) || (mother[p1] == mother[p2]))	tagHSo[born] = 1;
			}
			if (tracelevel!=0)
			{
				fprintf (fptr,"tagPSo=%d tagFSo=%d tagHSo=%d\n", tagPSo[born], tagFSo[born], tagHSo[born]);
				if(zygotes==0)	fprintf (fptr,"Born %d\n", born);
				else		fprintf (fptr,"Born %d  zyg%d\n", born, mark_zygb[i]);
			}
		}
		else
		{
			dead ++;
			if (tracelevel!=0)
			{
				if(zygotes==0)	fprintf (fptr,"Dead\n");
				else		fprintf (fptr,"Dead  zyg%d \n", mark_zygb[i]);
			}

			tries += 1;
			if (tries<1000)	goto generahijo;
			else 		goto bornprogeny;
		}
	}

	/******* born progeny *******/
	bornprogeny: /* ***** */;
	if (tracelevel!=0)	fprintf (fptr,"\n");

	for (i=0; i<born; i++)
	{
		for (k=0; k<NCRO; k++)
		{
		   	gm[i][k][0] = gmo[i][k][0];
		   	gm[i][k][1] = gmo[i][k][1];
		}

		father[i] = fathero[i];
		mother[i] = mothero[i];
		tagPS[i] = tagPSo[i];
		tagFS[i] = tagFSo[i];
		tagHS[i] = tagHSo[i];
		family[father[i]] += 1;
		family[mother[i]] += 1;
		if(zygotes==1)	mark_zyg[i] = mark_zygo[i];

		if (tracelevel!=0)	fprintf (fptr,"Born Progeny %d  p1 = %d    p2 = %d\n", i, father[i], mother[i]);
	}
	NIND[gen+1] = born;
	assign = born;
}

/* ***************************************************** */

mating_RCW()
{
	int rnd, a;
	double K;
	for (i=0; i<NIND[gen]; i++)	family[i] = 0;
	for (i=0; i<NMAX; i++)	failed[i] = 1;
	
	/******************************************/

	if (tracelevel!=0)	fprintf(fptr,"\nMatings MAXFEC %f\n", MAXFEC);

	/******************* Round of random contributions mating ***********************/

	if (tracelevel!=0)	fprintf (fptr,"\nRound of random contributions mating\n\n");

	born = 0;
	dead = 0;

	for (i=0; i<NIND[gen]; i+=2)
	{
		if (NIND[gen]==2)	{p1 = 0; p2 = 1;}
		else
		{
			p1 = fecunditysearch(cum_pm_sf, NIND[gen]);
			if (selfing==0)	do { p2 = fecunditysearch(cum_pm_sf, NIND[gen]); }	while (p2==p1);
			else	p2 = fecunditysearch(cum_pm_sf, NIND[gen]);
		}

		if (maternalf==0)	K = MAXFEC * sqrt(pm_sf[p1] * pm_sf[p2]);
		else			K = MAXFEC * pm_sf[p1];

		maxfec[i] = (int)(K*2);
//		maxfec[i] = (int)(K*2 + 0.5);
//		if (tracelevel!=0)	fprintf(fptr,"\nK  %f    macfec %d\n", K, maxfec[i]);

		for (j=0; j<maxfec[i]; j++)
		{
			if (tracelevel!=0)
			{
				if (zygotes==0)	fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d    p2 = %d\n", i, j, maxfec[i], p1, p2);
				else		fprintf(fptr,"Progeny %d - %d  maxfec %d  p1 = %d (zyg%d)    p2 = %d (zyg%d)\n", i, j, maxfec[i], p1, mark_zyg[p1], p2, mark_zyg[p2]);
			}
			gametes_offspring();
//			if ((zygotes==1)&&(rep == 1)) 	pedfile_zyg();
			bornoffspring();
		}
	}

	// First assignation of progeny to population

	if (tracelevel!=0)	fprintf (fptr,"First assignation of progeny to population\n\n");

	for (b=1; b<=born; b++)	mark[b] = 0;
	assign = 0;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	{
		if(zygotes==0)	fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d\n", b, fathero[b], mothero[b], mark[b]);
		else		fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d  zyg_n %d\n", b, fathero[b], mothero[b], mark[b], mark_zygo[b]);
	}
	if (tracelevel!=0)	fprintf (fptr,"\n");

	// Randomise born progeny

	for (b=1; b<=born; b++)
	{
		rnd = (int)(uniform()*born) + 1;

		for (k=0; k<NCRO; k++)
		{
		   	a = gmo[b][k][0]; gmo[b][k][0] = gmo[rnd][k][0]; gmo[rnd][k][0] = a;
		   	a = gmo[b][k][1]; gmo[b][k][1] = gmo[rnd][k][1]; gmo[rnd][k][1] = a;
		}
	   	a = fathero[b]; fathero[b] = fathero[rnd]; fathero[rnd] = a;
	   	a = mothero[b]; mothero[b] = mothero[rnd]; mothero[rnd] = a;
	   	a = tagPSo[b]; tagPSo[b] = tagPSo[rnd]; tagPSo[rnd] = a;
	   	a = tagFSo[b]; tagFSo[b] = tagFSo[rnd]; tagFSo[rnd] = a;
	   	a = tagHSo[b]; tagHSo[b] = tagHSo[rnd]; tagHSo[rnd] = a;
	   	a = mark[b]; mark[b] = mark[rnd]; mark[rnd] = a;
	   	if(zygotes==1)	a = mark_zygo[b]; mark_zygo[b] = mark_zygo[rnd]; mark_zygo[rnd] = a;
	}
	if (tracelevel!=0)	for (b=1; b<=born; b++)
	{
		if(zygotes==0)	fprintf (fptr,"Randomised born b %d  mark_b %d  fa %d  mo %d\n", b, mark[b], fathero[b], mothero[b]);
		else 		fprintf (fptr,"Randomised born b %d  mark_b %d  zyg_n %d  fa %d  mo %d\n", b, mark[b], mark_zygo[b], fathero[b], mothero[b]);
	}

	// Assignation

	if (assign == NMAX) goto endofmating_RCw;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"\nAvailable born before final assignment b %d  fa %d  mo %d  mark %d", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n\n");

	for (i=0; i<NMAX; i++)  if (failed[i] == 1)
	{
		if (assign == NMAX) goto endofmating_RCw;

		for (b=1; b<=born; b++)  if (mark[b] == 0)
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			if(zygotes==1)	mark_zyg[i] = mark_zygo[b];
			failed[i] = 0;
			if (tracelevel!=0)
			{
				if(zygotes==0)	fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
				else		fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  zyg_n %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], mark_zygo[b], fathero[b], mothero[b], failed[i]);
			}
			break;
		}
	}

	/******************* End of rounds of mating ***********************/

	endofmating_RCw: /******/;

	NIND[gen+1] = assign;
}

/* ***************************************************** */

printtriedoffspring()
{
	if (tracelevel!=0)
	{
		if (currenttype==2)	fprintf (fptr,"Progeny %d - %d  maxfec %d  p1 = %d    p2 = %d\n",
						i, j, maxfec[i], p1, p2);
		else	fprintf (fptr,"try %d  Tried parents  p1 = %d    p2 = %d   (p1f=%f p1v=%f p1at=%f p2f=%f p2v=%f p2at=%f)\n",
						j, p1, p2, pm_sf[p1], pm_sv[p1], pm_at[p1], pm_sf[p2], pm_sv[p2], pm_at[p2]);
	}
}

/* ***************************************************** */

gametes_offspring()
{
	int EE[MM], FF[MM], marker;
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;
	double rnd, sq_pm_a;

	/* ******************* Gametes with free recombination ******************* */

	    if(L==99.0)
	    {	    /* ******************* Free recombination ******************* */

		for (k=0; k<NCRO; k++)
		{
		   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			FF[k] = ~EE[k];
		   	gmb[j][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
		}
// COMMENT
//		if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n",
//		j, EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1],
//		sm[p1][2][0], sm[p1][2][1], gmb[j][0][0], gmb[j][1][0], gmb[j][2][0]);

		for (k=0; k<NCRO; k++)
		{
		   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
		   	FF[k] = ~EE[k];
		   	gmb[j][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
		}
// COMMENT
//		if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", j,
//		EE[0], EE[1], EE[2], sm[p2][0][0], sm[p2][0][1], sm[p2][1][0], sm[p2][1][1], 
//		sm[p2][2][0], sm[p2][2][1], gmb[j][0][1], gmb[j][1][1], gmb[j][2][1]);
	    }

	    else
	    {	    /* ************** Restricted recombination ***************** */

			/* ****** Chromosome from father ****** */

			numberrecs = recombinationnumber();

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)
				{
					EE[k] = ~EE[k];
				}
				FF[k] = ~EE[k];
				gmb[j][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}

// COMMENT
/*			if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n",
			j, EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1],
			sm[p1][2][0], sm[p1][2][1], gmb[j][0][0], gmb[j][1][0], gmb[j][2][0]);

			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
			/* ****** Chromosome from mother ****** */

			numberrecs = recombinationnumber();
			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}

			marker = 1;


			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)
				{
					EE[k] = ~EE[k];
				}
				FF[k] = ~EE[k];
				gmb[j][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
			}
// COMMENT
/*			if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", j,
			EE[0], EE[1], EE[2], sm[p2][0][0], sm[p2][0][1], sm[p2][1][0], sm[p2][1][1], 
			sm[p2][2][0], sm[p2][2][1], gmb[j][0][1], gmb[j][1][1], gmb[j][2][1]);

			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
	    }

	/* ****************** Genotypic value of offspring ****************** */

	pm_sfb[j] = addedfixed_sf;
	pm_svb[j] = addedfixed_sv;
	pm_atb[j] = addedfixed_at;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		// Trait
		if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_atb[j] += at[k][l];
		else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
		else	pm_atb[j] += (at[k][l]*hat[k][l]);

		// fecundity
		if (fec_via == 0)
		{
			if ( (l%2 == 0) && (initialgen[k][l] != (-99)) ) 
			{
   				if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_sfb[j] *= (1.0 + s[k][l]);
				else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_sfb[j] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		else if (fec_via == 1)
		{
			if ( (l%3 == 0) && (initialgen[k][l] != (-99)) ) 
			{
   				if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_sfb[j] *= (1.0 + s[k][l]);
				else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_sfb[j] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		// viability
		if (fec_via == 0)
		{
			if ( (l%2 != 0) && (initialgen[k][l] != (-99)) ) 
			if (initialgen[k][l] != (-99))
			{
				if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_svb[j] *= (1.0 + s[k][l]);
				else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_svb[j] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		else if (fec_via == 1)
		{
			if ( (l%3 != 0) && (initialgen[k][l] != (-99)) ) 
			if (initialgen[k][l] != (-99))
			{
				if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_svb[j] *= (1.0 + s[k][l]);
				else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_svb[j] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		// only viability
		else if (fec_via == 2)
		{
			if (initialgen[k][l] != (-99)) 
			{
				if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_svb[j] *= (1.0 + s[k][l]);
				else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_svb[j] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		// only fecundity
		else
		{
			if (initialgen[k][l] != (-99)) 
			{
				if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_sfb[j] *= (1.0 + s[k][l]);
				else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_sfb[j] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
	}

	// scaling by generation 0 fitness

	if (scaling == 1)
	{
		/* environmental source of mortality A = 0.3 */
		if (envmort==1)
		{
			if(Vs==0.0)
			{
				pm_sfb[j] = (pm_sfb[j] / Wfo) - 0.3;
				pm_svb[j] = (pm_svb[j] / Wvo) - 0.3;
			}
			else
			{
				sq_pm_a = pow ((pm_atb[j]-OPT), 2.0);
				pm_sfb[j] = (pm_sfb[j] * exp(-0.5*(sq_pm_a)/Vs) / Wfo) - 0.3;
				pm_svb[j] = (pm_svb[j] * exp(-0.5*(sq_pm_a)/Vs) / Wvo) - 0.3;
			}
		}
		else
		{
			if(Vs==0.0)
			{
				pm_sfb[j] = pm_sfb[j] / Wfo;
				pm_svb[j] = pm_svb[j] / Wvo;
			}
			else
			{
				sq_pm_a = pow ((pm_atb[j]-OPT), 2.0);
				pm_sfb[j] = pm_sfb[j] * exp(-0.5*(sq_pm_a)/Vs) / Wfo;
				pm_svb[j] = pm_svb[j] * exp(-0.5*(sq_pm_a)/Vs) / Wvo;
			}
		}
		if (pm_sfb[j] < 0.0)	pm_sfb[j] = 0.0;
		if (pm_sfb[j] > 1.0)	pm_sfb[j] = 1.0;
		if (pm_svb[j] < 0.0)	pm_svb[j] = 0.0;
		if (pm_svb[j] > 1.0)	pm_svb[j] = 1.0;
	}
	else
	{
		if(Vs==0.0)
		{
			pm_sfb[j] = pm_sfb[j];
			pm_svb[j] = pm_svb[j];
		}
		else
		{
			sq_pm_a = pow ((pm_atb[j]-OPT), 2.0);
			pm_sfb[j] = pm_sfb[j] * exp(-0.5*(sq_pm_a)/Vs);
			pm_svb[j] = pm_svb[j] * exp(-0.5*(sq_pm_a)/Vs);
		}
	}

	if (tracelevel!=0)	fprintf (fptr,"Offspring   fec = %f   via = %f   trait = %f\n", pm_sfb[j], pm_svb[j], pm_atb[j]);

	if (tracelevel!=0)	fprintf (fptr,"Offspring   born = %d   dead = %d\n", born, dead);

	accum (&gmean_sf_zygotes[gen], pm_sfb[j]);
	accum (&gmean_sv_zygotes[gen], pm_svb[j]);
	accum (&gmean_a_zygotes[gen], pm_atb[j]);

	if (tracelevel!=0)	fprintf (fptr,"gen = %d   mean_pm_sfb = %f   mean_pm_svb = %f   mean_pm_a = %f\n", gen, accmean(&gmean_sf_zygotes[gen]), accmean(&gmean_sv_zygotes[gen]), accmean(&gmean_a_zygotes[gen]));

	if (neutral == 1) pm_svb[j] = 1.0;
}

/* ***************************************************** */

int recombinationnumber ()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecombinantpoissontable, recombinantpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}

/* ***************************************************** */

bornoffspring()
{
	/********* Born or dead individuals *********/

	if (uniform() <= pm_svb[j])
	{
		born ++;
		for (k=0; k<NCRO; k++)
		{
		   	gmo[born][k][0] = gmb[j][k][0];
		   	gmo[born][k][1] = gmb[j][k][1];
		}
		fathero[born] = p1;
		mothero[born] = p2;
		if(zygotes==1)	mark_zygo[born] = mark_zygb[j];

		if (p1 == p2)	tagPSo[born] = 1;
		if (gen > 0)
		{	
			if ( ((father[p1] == father[p2]) && (mother[p1] == mother[p2])) ||
				((father[p1] == mother[p2]) && (mother[p1] == father[p2])) )	tagFSo[born] = 1;
			else if ((father[p1] == father[p2]) || (father[p1] == mother[p2]) ||
				(mother[p1] == father[p2]) || (mother[p1] == mother[p2]))	tagHSo[born] = 1;
		}
		if (tracelevel!=0)
		{
			fprintf (fptr,"tagPSo=%d tagFSo=%d tagHSo=%d\n", tagPSo[born], tagFSo[born], tagHSo[born]);
			if(zygotes==0)	fprintf (fptr,"Born %d\n", born);
			else		fprintf (fptr,"Born %d  zyg%d\n", born, mark_zygb[j]);
		}
	}
	else
	{
		dead ++;
		if (tracelevel!=0)
		{
			if(zygotes==0)	fprintf (fptr,"Dead\n");
			else		fprintf (fptr,"Dead  zyg%d\n", mark_zygb[j]);

		}
	}
}

/* ***************************************************** */

calculations()
{
	double family_sum=0.0, family_sum2=0.0;

	numPSborn=0.0;
	numFSborn=0.0;
	numHSborn=0.0;

	for (i=0; i<assign; i++)
	{
		numPSborn += (double)tagPS[i] / assign;
		numFSborn += (double)tagFS[i] / assign;
		numHSborn += (double)tagHS[i] / assign;
	}

	/************ P_SELF *************/

	if (assign > 0)		accum(&P_PSBORN[gen], numPSborn);
	if (assign > 0)		accum(&P_FSBORN[gen], numFSborn);
	if (assign > 0)		accum(&P_HSBORN[gen], numHSborn);

	/************ SK2, Ne(Sk2) and FVk *************/

	for (i=0; i<assign; i++)
	{
		family_sum += family[i];
		family_sum2 += (family[i] * family[i]);
	}
	sk2 = ((family_sum2-(family_sum*family_sum/assign))/assign);
	if (assign != 0)	accum(&SK2[gen], sk2);

	ne = (4.0 * assign) / (2.0 + sk2);
	if (assign != 0)	accum(&NeVk[gen], ne);

	if (currenttype==3)	alpha = PFSM / (4.0 - (3.0 * PFSM));
	else			alpha = 0.0;

	ne_alpha = (4.0 * assign) / ((2.0 * (1.0 - alpha)) + (sk2 * (1.0 + (3 * alpha))));
	if (assign != 0)	accum(&Nevk_alpha[gen], ne_alpha); 

	fvk = 1.0 - pow(1.0 - (1/(2.0*ne)), gen);
	if (assign != 0)	accum(&Fvk[gen], fvk);

	fvk_alpha = 1.0 - pow(1.0 - (1/(2.0*ne_alpha)), gen);
	if (assign != 0)	accum(&Fvk_alpha[gen], fvk_alpha);

	if (tracelevel!=0)  fprintf (fptr,"\nContributions from parents\n");
	if (tracelevel!=0)  for (i=0; i<NIND[gen]; i++)  fprintf (fptr,"family(%d) = %d\n", i, family[i]);
	if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", sk2);

	for (i=0; i<assign; i++) { tagPS[i]=0; tagFS[i]=0; tagHS[i]=0; }
	for (i=1; i<=born; i++) { tagPSo[i]=0; tagFSo[i]=0; tagHSo[i]=0; }

	/************ MEAN PRODUCTIVITY *************/

	accum (&PROD[gen], born);

}

/* ***************************************************** */

dumpoffspring()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Offspring before mutation (gm0 gm1)\n");	
	for (i=0; i<NIND[gen]; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}

/* ***************************************************** */

printout()
{
	fprintf(fgen,"\n\nNMAX=%d   type=%d   PFSM(3)=%4.2f  changed_to=%d  gen_change=%d	L=%4.2f	N.S-LOCI=%d    N.N-LOCI=%d    reps=%d   gens=%d\n    MAXFEC=%f    Lambda_s=%f   Lambda_L=%f\n    ave_s=%f   relax=%f    neutral=%d   scaling=%d   fec_via=%d   \n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f	A(0.3)=%d\n",
		NMAX, type, PFSM, CHT, CHTg, L, TOTLOCI-NCRO, NCRO, replicates, generations, MAXFEC, Lambda_s, Lambda_L, ave_s, relaxedsel, neutral, scaling, fec_via, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP,envmort);
	fprintf(fgen,"    rho = %f    ave_a = %f    ave_aL = %f    beta_a = %f    ave_ha = %f    P_a = %f    Vs = %f    OPT = %f    VE  = %f\n", rho, ave_a, ave_aL, beta_a, ave_ha, P_a, Vs, OPT, VE);

	fprintf(fgen,"LEQ(2dpq) = %f    LEQf = %f    LEQv = %f    TLEQ(sq) = %f\n", LEQ_NP, LEQ_NPf, LEQ_NPv, TLEQ_NP);
	fprintf(fgen,"TIME TO EXTINCTION = %f   var = %f    NSEGLOCNP = %d\n", accmean(&TIME_EXT), variance(&TIME_EXT), NSEGLOCNP);
	
	fprintf(fgen,"\ngen   PROD      Wfzyg      Wvzyg     N          fp        Fp        Fvk      F_PFSM     NeVk    Ne_PFSM      Wf        Wv       QT       var_QT    LEQ_NLf   LEQ_Lf    LEQ_NLv   LEQ_Lv    LEQ_T       SK2      FSborn    HSborn     Hw	     Ho	      DHWE	     Ho(del)	P_catastrophe\n");

	for (gen=0; gen<=generations; gen++)
	{
		if (accmean(&mean_N[gen]) > 0)
			fprintf(fgen, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
			gen, accmean(&PROD[gen]), accmean(&gmean_sf_zygotes[gen]), accmean(&gmean_sv_zygotes[gen]), accmean(&mean_N[gen]), accmean(&fp[gen]), accmean(&Fp[gen]), accmean(&Fvk[gen]), accmean(&Fvk_alpha[gen]), accmean(&NeVk[gen]), accmean(&Nevk_alpha[gen]), accmean(&gmean_sf[gen]), accmean(&gmean_sv[gen]), accmean(&gmean_at[gen]), accmean(&gvar_at[gen]),
			accmean(&LEQ_NLf[gen]), accmean(&LEQ_Lf[gen]), accmean(&LEQ_NLv[gen]), accmean(&LEQ_Lv[gen]), accmean(&LEQ_T[gen]), accmean(&SK2[gen]),
			accmean(&P_FSBORN[gen]), accmean(&P_HSBORN[gen]), accmean(&Hw[gen]), accmean(&Ho[gen]), 1.0 - (accmean(&Ho[gen])/accmean(&Hw[gen])), accmean(&Ho_del[gen]), accmean(&mean_pcat[gen]));
		else
			fprintf(fgen, "%d\n", gen);
	}

	//	outfile summary  Ls - Ll - s - hs - NMAX - MAXFEC - type - PFSM - gen - TimeEXT - N - Mf - Mv - LEQ_NLf - LEQ_Lf - LEQ_NLv - LE_Lv - LE_T - SK2 - selfB - Hw  

	// gen 0
//	fprintf(fout, "%f  %f  %f  %f  %d  %d  %d  %f  0  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", Lambda_s, Lambda_L, ave_s, ave_hs, NMAX, MAXFEC, type, PFSM,
//		accmean(&TIME_EXT), accmean(&mean_N[0]), accmean(&gmean_sf[0]), accmean(&gmean_sv[0]), accmean(&LEQ_NLf[0]),
//		accmean(&LEQ_Lf[0]), accmean(&LEQ_NLv[0]), accmean(&LEQ_Lv[0]), accmean(&LEQ_T[0]), accmean(&SK2[0]), 
//		accmean(&P_FSBORN[0]), accmean(&P_HSBORN[0]), accmean(&Hw[0]));
	// gen 10
//	fprintf(fout, "%f  %f  %f  %f  %d  %d  %d  %f  10  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", Lambda_s, Lambda_L, ave_s, ave_hs, NMAX, MAXFEC, type, PFSM,
//		accmean(&TIME_EXT), accmean(&mean_N[10]), accmean(&gmean_sf[10]), accmean(&gmean_sv[10]), accmean(&LEQ_NLf[10]),
//		accmean(&LEQ_Lf[10]), accmean(&LEQ_NLv[10]), accmean(&LEQ_Lv[10]), accmean(&LEQ_T[10]), accmean(&SK2[10]),
//		accmean(&P_FSBORN[10]), accmean(&P_HSBORN[10]), accmean(&Hw[10]));
	// gen generations
//	fprintf(fout, "%f  %f  %f  %f  %d  %d  %d  %f  %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", Lambda_s, Lambda_L, ave_s, ave_hs, NMAX, MAXFEC, type, PFSM,
//		generations, accmean(&TIME_EXT), accmean(&mean_N[generations]), accmean(&gmean_sf[generations]), accmean(&gmean_sv[generations]), accmean(&LEQ_NLf[generations]),
//		accmean(&LEQ_Lf[generations]), accmean(&LEQ_NLv[generations]), accmean(&LEQ_T[generations]), accmean(&LEQ_Lv[generations]),
//		accmean(&P_FSBORN[generations]), accmean(&P_HSBORN[generations]), accmean(&Hw[generations]));
//	fprintf(fout, "\n\n\n");
}

/* ***************************************************** */

distribution_out()
{
	fprintf (fgen, "\n");
	for (i=0; i<=(generations/tgen); i++) 	fprintf (fgen, "Extinction risk gen %d = %f\n", i*tgen, (replicates-reps_alive_gen[i])*100.0/replicates);

	fprintf (fgen, "\nGen		");
	for (i=0; i<=(generations/tgen); i++) 	fprintf (fgen, "%d         ", i*tgen);
	fprintf (fgen, "\n");

	// distribution of q values in the population

		fprintf (fgen, "q=0.0          ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_00[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.0-0.1      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_00_01[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.1-0.2      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_01_02[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.2-0.3      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_02_03[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.3-0.4      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_03_04[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.4-0.5      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_04_05[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.5-0.6      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_05_06[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.6-0.7      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_06_07[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.7-0.8      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_07_08[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.8-0.9      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_08_09[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.9-1.0      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_09_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=1.0          ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

	// distribution of s values in the population

		fprintf (fgen, "LOST\n");

		fprintf (fgen, "s=0.0         ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_106_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_106_104_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.01 ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_104_102_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.01-0.1  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_102_101_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_01_02_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_02_04_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_04_06_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_06_10_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_10_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "SEGREGATING\n");

		fprintf (fgen, "s=0.0         ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_106[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_106_104[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.001 ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_104_102[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.001-0.1  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_102_101[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_01_02[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_02_04[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_04_06[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_06_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "FIXED\n");

		fprintf (fgen, "s=0.0         ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_106_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_106_104_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.001 ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_104_102_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.001-0.1  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_102_101_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_01_02_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_02_04_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_04_06_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_06_10_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_10_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n\n");

	// distribution of hs values in the population

		fprintf (fgen, "LOST\n");

		fprintf (fgen, "hs=0.0-0.2      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_00_02_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.2-0.4      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_02_04_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.4-0.6      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_04_06_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.6-0.8      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_06_08_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.8-1.0      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_08_10_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "SEGREGATING\n");

		fprintf (fgen, "hs=0.0-0.2      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_00_02[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.2-0.4      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_02_04[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.4-0.6      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_04_06[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.6-0.8      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_06_08[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.8-1.0      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_08_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "FIXED\n");

		fprintf (fgen, "hs=0.0-0.2      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_00_02_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.2-0.4      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_02_04_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.4-0.6      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_04_06_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.6-0.8      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_06_08_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.8-1.0      ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_08_10_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

	// distribution of s (q and h values) in the population

		fprintf (fgen, "\nSEGREGATING - frequencies (q)\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_0_106[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_106_104[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.001 ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_104_102[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.001-0.1  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_102_101[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_01_02[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_02_04[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_04_06[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_06_10[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accmean(&qs_ng_10[i]));
		fprintf (fgen, "\n");


		fprintf (fgen, "\nSEGREGATING - dominance (h)\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_0_106[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_106_104[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.001 ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_104_102[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.001-0.1  ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_102_101[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_01_02[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_02_04[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_04_06[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_06_10[i]));
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=(generations/tgen); i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.6f  ", accmean(&hs_10[i]));
		fprintf (fgen, "\n");

}

/* ***************************************************** */


