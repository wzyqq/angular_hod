#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_expint.h>
#include <unistd.h>
#include "mpi.h"

//#define K_num  237
//#define K_num  10001
#define header_number 200
#define file_char_num 200
#define r_num  500
#define cut_number 15
#define node   16
#define y_node 100000
#define lightspeed 299792.458
#define pi     3.141592653589793
#define G   6.67259E-11
#define little_h     0.7
#define H0    100*little_h
#define delta_c   1.686
//#define z   1.75
//#define Omega_m0   0.3
//#define Omega_lambda0   0.7
#define triangle   200.0
#define theta_max 20
#define theta_observe   16
#define sigma_log_M 0.2
#define alpha 1.0
//#define sigma_log_M 0.6
//#define alpha 1.6
#define m_max   70          //compare to m_max=70000,ka_square 30.0476 and 30.0458 ,but when the largest m is 15,ka_square is about 38,big difference
#define ndim  2
//#define walks  50000
#define walks  1
#define chains 1
//#define r_test 431


char Prefix[file_char_num];
char Linear_matter_power_spectrum[file_char_num];     
char Nonlinear_matter_power_spectrum[file_char_num]; 
int  Start_point;
int  Theta_min;	
char W_mean[file_char_num];
char Inverse_matrix[file_char_num];
double Number_density;			
double Sigma_ng;	
int K_num;
double Omega_m0;
double Omega_lambda0;
char Redshift_file[file_char_num];
double y_low,y_high;
long long Rr_pair[theta_max]={0};
double Norm=0;


double M[m_max]={0},R[m_max]={0},M_add[m_max+1]={0},R_add[m_max+1]={0};
double *Sigma,*Sigma_add, *F, *Dlninvsigma_divide_dm, *Dn_divide_dm, *Dn_divide_dm_eff, *Bias, *Bias_eff, *Dr_divide_dz;
double Xi[node]={0},Wi[node]={0};
double *Y, *Yweight;
double R_corr[r_num]={0};
double *K, *Pkl, *Pk;
double *Theta;
double Rh=0;
int Redshift_number;
double *Redshift;
double *Redshift_distace;
double *Redshift_distribution;
double *Redshift_weight;
double *Dredshift_divide_dr;
double *Corr_mm;
int r_cut[cut_number+1]= {0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 150, 180,  250,  300,  400, r_num};
//int r_cut_value[cut_number]={15, 20, 25, 30, 40, 50, 60, 80, 100, 100, 500, 1000, 5000, 10000, 100000};
int r_cut_value[cut_number]={15, 20, 25, 30, 40, 50, 60, 80, 100, 100, 500, 1000, 3000, 5000, 5000};


struct Linear
{		
	double kl_r;
	double pkl_r;
};
		
struct Nonlinear
{		
	double k_r;
	double pk_r;
};



double Rho_m(double redshift)
{

	double E,Omega_m,Omega_lambda,H;
	double rho_m;
	E = sqrt(Omega_lambda0+Omega_m0*pow(1+redshift,3));
//	WRITE(*,*)Omega_m0,Omega_lambda0
//	printf(*,*)'E=',E
    Omega_m = Omega_m0*pow((1+redshift),3)/pow(E,2);
//	 printf(*,*)'Omega_m=',Omega_m
    Omega_lambda = Omega_lambda0/pow(E,2);
//	printf(*,*)'Omega_lambda=',Omega_lambda
    H = H0*E*1000/(3.08568E22);
//	printf(*,*)'H=',H
    rho_m = 3.0*H*H/(8.0*pi*G);
 // printf(*,*)'Rho_m=',Rho_m
    return  rho_m;
}



void Sigma_f(double *sigma_s, double *R, int length)
{
    int i,j,mum,num;
	double xi_res[node],wi_res[node];
	gsl_interp_accel *acc0[Redshift_number];
	gsl_spline *spline0[Redshift_number];
	
	for(i=0; i<Redshift_number; i++)
	{	
		acc0[i]= gsl_interp_accel_alloc ();
		spline0[i]= gsl_spline_alloc (gsl_interp_cspline, K_num);
		gsl_spline_init (spline0[i], &K[i*K_num], &Pkl[i*K_num], K_num);

		for(j=0; j<length; j++)
		{
			sigma_s[i*length+j] = 0;
			for(mum=0; mum<100; mum++)
			{	
				for(num=0; num<node; num++)
				{
					xi_res[num] = Xi[num] + 2*pi*mum;
					wi_res[num] = Wi[num];
					sigma_s[i*length+j] = sigma_s[i*length+j] + wi_res[num]*xi_res[num]*xi_res[num]/R[j]/R[j]/R[j]*\
					pow(3.0*(sin(xi_res[num])-xi_res[num]*cos(xi_res[num]))/pow(xi_res[num],3),2)*gsl_spline_eval(spline0[i], xi_res[num]/R[j], acc0[i]);
				}
			} 
			sigma_s[i*length+j] = sqrt(sigma_s[i*length+j]/2.0/pi/pi);
		}
		gsl_spline_free (spline0[i]);
		gsl_interp_accel_free (acc0[i]);
	}
}

void Multiplicity(double *f)
{
	int i,j;
	double A[Redshift_number],a_temp[Redshift_number],beta[Redshift_number],b[Redshift_number],c[Redshift_number];
	for(i=0; i<Redshift_number; i++)
	{
//		A = (0.1*log10(triangle)-0.05)*pow(1+redshift,-0.14);
		A[i] = 0.186*pow(1+Redshift[i],-0.14);
//		printf(*,*)"A=",A;
//		a_temp = (1.43+pow(log10(triangle)-2.3,1.5))*pow(1+redshift,-0.06);
		a_temp[i] = 1.47*pow(1+Redshift[i],-0.06);
// 		printf(*,*)"a_temp=",a_temp;
    	beta[i] = pow(10,-pow(0.75/log10(triangle/75.0),1.2));
//		printf(*,*)"beta=",beta;
//		b = (1.0+pow(log10(triangle)-1.6,-1.5))*pow(1+redshift,-beta);
		b[i] = 2.57*pow(1+Redshift[i],-beta[i]);
//		printf(*,*)"b=",b;
//		c = 1.2+pow(log10(triangle)-2.35,1.6);
		c[i] = 1.19;
//		printf(*,*)"c=",c;
//		printf("red%lf\n", Redshift[i]);
//		printf("%lf\t%lf\t%lf\t%lf\t%lf\n", A[i], a_temp[i], beta[i], b[i], c[i]);
    	for(j=0; j<m_max; j++)
    	{
    		f[i*m_max+j] = A[i]*(pow(Sigma[i*m_max+j]/b[i],-a_temp[i])+1)*exp(-c[i]/pow(Sigma[i*m_max+j],2));
	    }
	}
}

void bias_f(double *bias)
{
	int i,j;
	double nu[Redshift_number*m_max];
	double y, A, a_temp, B, b_temp, C, c_temp;
    y = log10(triangle);
    A = 1 + 0.24*y*exp(-pow(4.0/y,4));
    a_temp = 0.44*y - 0.88;
    B = 0.183;
    b_temp = 1.5;
    C = 0.019 + 0.107*y +0.19*exp(-pow(4.0/y,4));
    c_temp = 2.4;
    for(i=0; i<Redshift_number; i++)
	{
		for(j=0; j<m_max; j++)
		{
			nu[i*m_max+j] = delta_c/Sigma[i*m_max+j];      
			bias[i*m_max+j] = 1 - A*pow(nu[i*m_max+j],a_temp)/(pow(nu[i*m_max+j],a_temp)+pow(delta_c,a_temp)) + B*pow(nu[i*m_max+j],b_temp) + C*pow(nu[i*m_max+j],c_temp);
    	}
    }
}

void profile(double *prof, double *rvir)
{
	int i,j,s;
	double m_nl=3.9E12;
    double c[Redshift_number*m_max], rho_h, rs[Redshift_number*m_max], rhos[Redshift_number*m_max];
    double *si_more, *ci_more, *si, *ci;
   	si_more = (double *)malloc(Redshift_number*m_max*K_num*sizeof(double));
   	ci_more = (double *)malloc(Redshift_number*m_max*K_num*sizeof(double));
   	si = (double *)malloc(Redshift_number*m_max*K_num*sizeof(double));
   	ci = (double *)malloc(Redshift_number*m_max*K_num*sizeof(double));
    rho_h = triangle*Rh;
    for(s=0; s<m_max; s++)
    {
		rvir[s] = pow(3.0*M[s]/4.0/pi/rho_h,1.0/3.0);
    }
	for(i=0; i<Redshift_number; i++)
	{
		for(s=0; s<m_max; s++)
		{
//			printf("red%lf\n",Redshift[i]);
//			printf("m_nl%le\n", m_nl[i]);
//			printf("m%le\n", M[s]);
			c[i*m_max+s] = 11.0*pow(M[s]/m_nl,-0.13)/(1+Redshift[i]);
			rs[i*m_max+s] = rvir[s]/c[i*m_max+s];
			rhos[i*m_max+s] = Rh*triangle/3.0*pow(c[i*m_max+s],3)/(log(1+c[i*m_max+s])-c[i*m_max+s]/(1+c[i*m_max+s]));
//			printf("rho_h = %lf\trvir= %lf\tc = %lf\trs = %lf\trhos = %lf\n ", rho_h, rvir[s],  c[i*m_max+s], rs[i*m_max+s], rhos[i*m_max+s]);
    		for(j=0; j<K_num; j++)
    		{
				si_more[i*m_max*K_num+s*K_num+j] = gsl_sf_Si((1+c[i*m_max+s])*K[i*K_num+j]*rs[i*m_max+s]);
				ci_more[i*m_max*K_num+s*K_num+j] = gsl_sf_Ci((1+c[i*m_max+s])*K[i*K_num+j]*rs[i*m_max+s]);               
//				printf(*,*)'si_more=',si_more
//				printf(*,*)'ci_more=',ci_more
				si[i*m_max*K_num+s*K_num+j] = gsl_sf_Si(K[i*K_num+j]*rs[i*m_max+s]);
				ci[i*m_max*K_num+s*K_num+j] = gsl_sf_Ci(K[i*K_num+j]*rs[i*m_max+s]);
//				printf(*,*)'si=',si
//				printf(*,*)'ci=',ci
				prof[i*m_max*K_num+s*K_num+j] = 4.0*pi*rhos[i*m_max+s]*pow(rs[i*m_max+s],3)/M[s]*
				( sin(K[i*K_num+j]*rs[i*m_max+s]) * (si_more[i*m_max*K_num+s*K_num+j]-si[i*m_max*K_num+s*K_num+j]) - 
				sin(c[i*m_max+s]*K[i*K_num+j]*rs[i*m_max+s])/(1+c[i*m_max+s])/K[i*K_num+j]/rs[i*m_max+s] + cos(K[i*K_num+j]*rs[i*m_max+s]) 
				* (ci_more[i*m_max*K_num+s*K_num+j]-ci[i*m_max*K_num+s*K_num+j]) ); 
//				printf("prof[%d] = %.10lf\n", j, (*prof)[j]);
			}
		}
	}
	free(si_more);
	free(ci_more);
	free(si);
	free(ci);
}

void exclusion(double *exc, double *rvir)
{
	int i,j;
	int loc;
	for(i=0; i<m_max; i++)
	{
		loc = (int)100.0*(log10(2.0*rvir[i])+2.005);
		for(j=0; j<r_num; j++)
		{
			if(j<=loc)
			{
				exc[i*r_num+j] = 0;
			}
			else 
			{
				exc[i*r_num+j] = 1;
			}
		}
	}
}


double get_w(double *param, double *w_analysis, double *halo_mass, double *beff ,double *fsat)
{
//	double ng=0;
//	FILE *fcorr;
//	FILE *fwrong;
	int i, j, mum, num, pa, pe, loop, s;
	double b[Redshift_number];
	double ng[Redshift_number],nc[Redshift_number],ns[Redshift_number];
	double w_true[theta_max]={0};
	double temp1[m_max]={0};
	double N_cen_average[m_max]={0},N_sat_average[m_max]={0},N_average[m_max]={0};
	double rvir[m_max]={0};
	double *u;
	double *p1;
	double *pcs1, *pss1, *p;
	double exc[m_max*r_num],fexc[Redshift_number*r_num];
	double corr_1h[Redshift_number*r_num];
	double corr_2hnew[Redshift_number*r_num];
	double log_M_min, M0, M1_primer;
	double xi_res[node]={0},wi_res[node]={0};
	double corrrrr[Redshift_number*r_num];
	double n_all=0;
	double *r12;
	double *temp3;
	int *location;
	double IC=0;
	long long sum=0;
	gsl_interp_accel *acc2[Redshift_number];
	gsl_spline *spline2[Redshift_number];

//	printf("good\n");


	p1 = (double *)malloc(Redshift_number*K_num*sizeof(double));
	pcs1 = (double *)malloc(Redshift_number*K_num*sizeof(double));
	pss1 = (double *)malloc(Redshift_number*K_num*sizeof(double));
	p = (double *)malloc(Redshift_number*K_num*sizeof(double));
	u = (double *)malloc(Redshift_number*m_max*K_num*sizeof(double));

	for(i=0; i<Redshift_number; i++)
	{
		for(j=0; j<K_num; j++)
		{
			p1[i*K_num+j] = 0;
			pcs1[i*K_num+j] = 0;
			pss1[i*K_num+j] = 0;
			p[i*K_num+j] = 0;
		}
		for(j=0; j<r_num; j++)
		{
			corr_1h[i*r_num+j]=0;
			corr_2hnew[i*r_num+j]=0;
			corrrrr[i*r_num+j]=0;
			fexc[i*r_num+j]=0;
		}
	}

	
//	FILE *ft;
//	double corrrrr[r_test]={0};

	log_M_min = *(param+0);
	M1_primer = pow( 10.0,*(param+1) );
	M0 		  = pow( 10.0,0.76*(*(param+1))+2.3 );
//	printf("log_M_min=%lf\tsigma_log_M=%lf\tM0=%le\tM1_primer=%le\talpha=%lf\n", log_M_min, sigma_log_M, M0, M1_primer, alpha);
	profile(u, rvir);
	exclusion(exc, rvir);
	for(i=0;i<m_max;i++)
	{
		temp1[i] = (log10(M[i]) - log_M_min)/sigma_log_M;
//		printf(*,*)'temp1=',temp1
		N_cen_average[i] = (gsl_sf_erf(temp1[i])+1.0)/2.0;
//		printf(*,*)'N_cen_average=',N_cen_average   
		if((M[i] - M0)<=0)
		{
			N_sat_average[i]=0;
		}
		else        
		{
			N_sat_average[i] = N_cen_average[i]*pow((M[i] - M0)/M1_primer,alpha); 		
		}
//		printf(*,*)'N_sat_average=',N_sat_average
		N_average[i] = N_cen_average[i] + N_sat_average[i];
//		printf(*,*)'N_average=',N_average

//		printf("N_cen_average[%d] = %le\tN_sat_average[%d] = %le\tN_average[%d] = %le\n", i, N_cen_average[i], i, N_sat_average[i], i, N_average[i]);
//		printf("rvir[%d] = %lf\n", i, rvir[i]);
//		for(j=0; j<r_num; j++)
//		{
//			printf("exc[%d][%d]=%lf\n", i, j ,exc[i*r_num+j]);
//		}
	}

	for(i=0; i<Redshift_number; i++)
	{
//		for(j=0; j<m_max; j++)
//		{
//			for(s=0; s<K_num; s++)
//			{
//				printf("u[%d][%d][%d]=%lf\n", i, j, s, u[i*m_max*K_num+j*K_num+s]);
//			}
//		}


		nc[i]=0;
		ns[i]=0;
		ng[i]=0;
		for(j=0; j<(m_max-1); j++)
		{
			nc[i] = nc[i] + (Dn_divide_dm[i*m_max+j]*N_cen_average[j] + Dn_divide_dm[i*m_max+j+1]*N_cen_average[j+1])*(M[j+1]-M[j]);
			ns[i] = ns[i] + (Dn_divide_dm[i*m_max+j]*N_sat_average[j] + Dn_divide_dm[i*m_max+j+1]*N_sat_average[j+1])*(M[j+1]-M[j]);
    		ng[i] = ng[i] + (Dn_divide_dm[i*m_max+j]*N_average[j] + Dn_divide_dm[i*m_max+j+1]*N_average[j+1])*(M[j+1]-M[j]);
		}
		nc[i] = nc[i]/2.0;
		ns[i] = ns[i]/2.0;
		ng[i] = ng[i]/2.0;
//		printf("nc = %lf\n",nc[i]);
//		printf("ns = %lf\n",ns[i]);
//		printf("ng = %lf\n",ng[i]);		
		n_all = n_all + Redshift_distace[i]*Redshift_distace[i]*Dr_divide_dz[i]*ng[i]*Redshift_weight[i];    //Redshift_weight is normalized and Redshift_weight = z_distribution.
	}

	n_all = n_all/Norm;
//	printf("n_all=%lf\n", n_all);

	for(i=0; i<Redshift_number; i++)
	{
		b[i] = 0;
		for(j=0; j<(m_max-1); j++)
		{
    		b[i] = b[i] + (Dn_divide_dm[i*m_max+j]*N_average[j]*Bias[i*m_max+j] + Dn_divide_dm[i*m_max+j+1]*N_average[j+1]*Bias[i*m_max+j+1])
    		*(M[j+1]-M[j])/2.0;
		}
		*beff = *beff + Redshift_distace[i]*Redshift_distace[i]*Dr_divide_dz[i]*b[i]*Redshift_weight[i];    //Redshift_weight is normalized and Redshift_weight = z_distribution.
	}
	*beff = *beff/Norm/n_all;

	for(i=0; i<(m_max-1); i++)
	{	
		*halo_mass = *halo_mass + (M[i+1]-M[i])
		*(Dn_divide_dm_eff[i+1]*N_average[i+1]*M[i+1]+Dn_divide_dm_eff[i]*N_average[i]*M[i])/n_all/2.0;
//		*beff      = *beff      + (M[i+1]-M[i])
//		*(Dn_divide_dm_eff[i+1]*N_average[i+1]*Bias_eff[i+1]+Dn_divide_dm_eff[i]*N_average[i]*Bias_eff[i])/n_all/2.0;
		*fsat      = *fsat      + (M[i+1]-M[i])
		*(Dn_divide_dm_eff[i+1]*N_sat_average[i+1]+Dn_divide_dm_eff[i]*N_sat_average[i])/n_all/2.0;
	}
//	printf("halo_mass=%lf\n", log10(*halo_mass));
//	printf("beff=%lf\n", *beff);
//	printf("fsat=%lf\n", *fsat);



	for(s=0; s<Redshift_number; s++)
	{
		for(j=0; j<K_num; j++)
		{
			for(i=0; i<(m_max-1); i++)
    		{
    			pcs1[s*K_num+j] = pcs1[s*K_num+j] + (N_cen_average[i]*N_sat_average[i]*Dn_divide_dm[s*m_max+i]*u[s*m_max*K_num+i*K_num+j]+\
    				N_cen_average[i+1]*N_sat_average[i+1]*Dn_divide_dm[s*m_max+i+1]*u[s*m_max*K_num+(i+1)*K_num+j])*(M[i+1]-M[i]);
    			pss1[s*K_num+j] = pss1[s*K_num+j] + (N_sat_average[i]*N_sat_average[i]*Dn_divide_dm[s*m_max+i]*u[s*m_max*K_num+i*K_num+j]*u[s*m_max*K_num+i*K_num+j]+\
    	    		N_sat_average[i+1]*N_sat_average[i+1]*Dn_divide_dm[s*m_max+i+1]*u[s*m_max*K_num+(i+1)*K_num+j]*u[s*m_max*K_num+(i+1)*K_num+j])*(M[i+1]-M[i]); 
			}
			pcs1[s*K_num+j] = pcs1[s*K_num+j]/ng[s]/ng[s];
			pss1[s*K_num+j] = pss1[s*K_num+j]/ng[s]/ng[s]/2.0;
			p1[s*K_num+j] =  pcs1[s*K_num+j] + pss1[s*K_num+j];
//			printf("pcs1[%d] = %.8lf\tpss1[%d] = %.8lf\tp1[%d] = %.8lf\n",j, pcs1[j], j, pss1[j], j, p1[j]);
		}
	}
	
	for(s=0; s<Redshift_number; s++)
	{	
		acc2[s]= gsl_interp_accel_alloc ();
		spline2[s]= gsl_spline_alloc (gsl_interp_cspline, K_num);
		gsl_spline_init (spline2[s], &K[s*K_num], &p1[s*K_num], K_num);
		for(mum=0; mum<cut_number; mum++)
		{	
			for(j=r_cut[mum]; j<r_cut[mum+1]; j++)
			{
				for(i=0; i<r_cut_value[mum]; i++)
				{	
					for(num=0; num<node; num++)
					{
						xi_res[num] = Xi[num] + 2*pi*i;
						wi_res[num] = Wi[num];
						corr_1h[s*r_num+j] = corr_1h[s*r_num+j] + wi_res[num]*xi_res[num]/R_corr[j]/R_corr[j]/R_corr[j]*gsl_spline_eval(spline2[s], xi_res[num]/R_corr[j], acc2[s])*sin(xi_res[num]);
					}
				}
				corr_1h[s*r_num+j] = corr_1h[s*r_num+j]/2.0/pi/pi;
			}
		}
		gsl_spline_free (spline2[s]);
		gsl_interp_accel_free (acc2[s]);


		for(j=0; j<r_num; j++)
		{
			for(i=0; i<(m_max-1); i++)
			{
				fexc[s*r_num+j] = fexc[s*r_num+j] + (Dn_divide_dm[s*m_max+i]*N_average[i]*Bias[s*m_max+i]*exc[i*r_num+j] + 
					Dn_divide_dm[s*m_max+i+1]*N_average[i+1]*Bias[s*m_max+i+1]*exc[(i+1)*r_num+j])*(M[i+1]-M[i]);
			}
			fexc[s*r_num+j] = fexc[s*r_num+j]/ng[s]/2.0;
			corr_2hnew[s*r_num+j] = corr_2hnew[s*r_num+j] + fexc[s*r_num+j]*fexc[s*r_num+j]*Corr_mm[s*r_num+j];
			corrrrr[s*r_num+j] = corr_1h[s*r_num+j] + corr_2hnew[s*r_num+j];
//			printf("fexc[%d] = %.8lf\n", j, fexc[j]);
//			printf("%d\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", j, corr_2h[j], corr_2hold1[j], corr_2hold2[j], corr_2hnew[j], corr_1h[j]+corr_2hnew[j]);
//			printf("%d\t%.8lf\t%.8lf\t%.8lf\n", j, corr_1h[s*r_num+j], corr_2hnew[s*r_num+j], corr_1h[s*r_num+j]+corr_2hnew[s*r_num+j]);
		}
	}
//
/*	for(s=0; s<Redshift_number; s++)
	{
		fcorr = fopen("xi_camb","r");
		for(i=0; i<r_num; i++)
		{
			fscanf(fcorr, "%lf", &corrrrr[s*r_num+i]);
			printf("%lf\n", corrrrr[s*r_num+i]);
		}
	}
	
	fwrong = fopen("redshift_info", "r");	
	for(i=0; i<Redshift_number; i++)
	{
		fscanf(fwrong, "%lf%lf%lf%lf", &Redshift[i], &Redshift_weight[i], &Redshift_distace[i], &Redshift_distribution[i]);		
		Dr_divide_dz[i] = 100.0/299792.458;	
		printf("%lf\t%lf\t%lf\t%lf\n", Redshift[i], Redshift_weight[i], Redshift_distace[i], Redshift_distribution[i]);
	}
*/
	r12 = (double *)malloc(theta_max*Redshift_number*y_node*sizeof(double));
	location = (int *)malloc(theta_max*Redshift_number*y_node*sizeof(int));
	temp3 = (double *)malloc(theta_max*Redshift_number*sizeof(double));

	for(j=0; j<Theta_min; j++)
	{
		w_analysis[j] = 0;
	}
	for(j=0; j<theta_max*Redshift_number; j++)
	{
		temp3[j] = 0;
	}

	for(s=0; s<theta_max; s++)
	{
		for(i=0; i<Redshift_number; i++)
		{
			for(j=0; j<y_node; j++)
			{
				r12[s*Redshift_number*y_node+i*y_node+j] = sqrt(Y[j]*Y[j]+Redshift_distace[i]*Redshift_distace[i]*Theta[s]*Theta[s]);
				location[s*Redshift_number*y_node+i*y_node+j] = (int)(100*log10(r12[s*Redshift_number*y_node+i*y_node+j])+200.5);
				if(location[s*Redshift_number*y_node+i*y_node+j]<0)
				{	
					location[s*Redshift_number*y_node+i*y_node+j] = 0;
				}
				else if(location[s*Redshift_number*y_node+i*y_node+j]>=r_num)
				{
					location[s*Redshift_number*y_node+i*y_node+j] = r_num-1;
				}
				temp3[s*Redshift_number+i] = temp3[s*Redshift_number+i] + corrrrr[i*r_num+location[s*Redshift_number*y_node+i*y_node+j]]*Yweight[j];
			}
//			if(s==0)
//			{
//				printf("temp3[0][%d]=%lf\n", i, temp3[0][i]);
//			}
			w_true[s] = w_true[s] + temp3[s*Redshift_number+i]*Redshift_distribution[i]*Redshift_distribution[i]*Dr_divide_dz[i]*Redshift_weight[i];
		}
		printf("%lf\t", w_true[s]);
	}

	for(i=0; i<theta_max; i++)
	{
		sum = sum + Rr_pair[i];
		IC = IC + Rr_pair[i]*w_true[i];
	}
//	printf("pre_IC=%lf\n", IC);	
//	printf("sum=%lld\n", sum);	
	IC = IC*1.0/sum;
	printf("IC=%lf\n", IC);
	for(i=0; i<Theta_min; i++)
	{
		w_analysis[i] = w_true[i+Start_point] - IC;
	}

	free(p1);
	free(pcs1);
	free(pss1);
	free(p);
	free(u);
	free(r12);
	free(location);
	free(temp3);

	return n_all;
}





void getMCMC(double *init, int *argc_address, char ***argv_address)
{
	FILE *foutput;
	FILE *fmean, *finvconvar;
	int num,i,j,loop;
	double *ng;
	double *halo_mass;
	double *beff;
	double *fsat;	
	double mag[ndim]={0.01,0.01};
	double random0,random1[ndim]={0}, random2[ndim]={0}, rand_delta[ndim]={0};
	double *temp;
	double *delta_w;
	double random_xi;
	double *observe,*w_observe;
	double *covariance_inv;
	double judge;
	double *pos;
	double *w_analysis;
	double *ka_square;
	double *result1, *result2, *result3, *result4, *result5, *result6, *result7;
	double *result_log_M_min;
	double *result_M1_primer;
	double *result_ng;
	double *result_halo_mass;
	double *result_beff;
	double *result_fsat;	
	double *result_ka;
	int root =0;
	int myid;
		
	temp = 				(double *)malloc(Theta_min*sizeof(double));
	observe = 			(double *)malloc(theta_observe*sizeof(double));
	w_observe = 		(double *)malloc(Theta_min*sizeof(double));
	w_analysis =		(double *)malloc(Theta_min*sizeof(double));
	covariance_inv =	(double *)malloc(Theta_min*Theta_min*sizeof(double));
	pos = 				(double *)malloc(walks*ndim*sizeof(double *));
	halo_mass = 		(double *)malloc(walks*sizeof(double));
	beff = 				(double *)malloc(walks*sizeof(double));
	fsat =				(double *)malloc(walks*sizeof(double));
	ng = 				(double *)malloc(walks*sizeof(double));
	delta_w = 			(double *)malloc(walks*sizeof(double));
	ka_square = 		(double *)malloc(walks*sizeof(double));
	result1 = 			(double *)malloc(walks*sizeof(double));
	result2 = 			(double *)malloc(walks*sizeof(double));
	result3 = 			(double *)malloc(walks*sizeof(double));
	result4 = 			(double *)malloc(walks*sizeof(double));
	result5 = 			(double *)malloc(walks*sizeof(double));
	result6 = 			(double *)malloc(walks*sizeof(double));
	result7 = 			(double *)malloc(walks*sizeof(double));	
	result_log_M_min =	(double *)malloc(walks*chains*sizeof(double));
	result_M1_primer =	(double *)malloc(walks*chains*sizeof(double));
	result_ng = 		(double *)malloc(walks*chains*sizeof(double));
	result_halo_mass =	(double *)malloc(walks*chains*sizeof(double));
	result_beff =		(double *)malloc(walks*chains*sizeof(double));
	result_fsat = 		(double *)malloc(walks*chains*sizeof(double));	
	result_ka = 		(double *)malloc(walks*chains*sizeof(double));

	for(i=0; i<walks; i++)
	{
		ng[i]       =0;
		halo_mass[i]=0;
		beff[i]     =0;
		fsat[i]     =0;
		delta_w[i]  =0;
		ka_square[i]=0;
		result1[i]  =0;
		result2[i]  =0;
		result3[i]  =0;
		result4[i]  =0;
		result5[i]  =0;
		result6[i]  =0;
		result7[i]  =0;		
	}

	for(i=0; i<walks*chains; i++)
	{
		result_log_M_min[i]=0;
		result_M1_primer[i]=0;
		result_ng[i]=0;
		result_halo_mass[i] =0;
		result_beff[i] =0;
		result_fsat[i] =0;			
		result_ka[i]=0;
	}

	for(i=0; i<Theta_min; i++)
	{
		w_analysis[i]=0;
		temp[i]=0;
//		printf("temp[%d]=%lf\n",i, temp[i]);
	}

	fmean = fopen(W_mean, "r");
	finvconvar = fopen(Inverse_matrix, "r");
	
	for(i=0; i<theta_observe; i++)
	{	
		fscanf(fmean, "%lf", &observe[i]);
//		printf("observe[%d]=%lf\n", i, observe[i]);
	}
	fclose(fmean);

	for(i=0; i<Theta_min; i++)
	{	
		w_observe[i] = observe[Start_point+i];
//		printf("w_observe[%d] = %lf\n", i, w_observe[i]);
	}


	for(i=0; i<Theta_min; i++)
	{
		for(j=0; j<Theta_min; j++)
		{
			fscanf(finvconvar, "%lf", &covariance_inv[i*Theta_min+j]);
//			printf("%lf\t", covariance_inv[i*Theta_min+j]);
		}
//		printf("\n");
	}
	fclose(finvconvar);

	MPI_Init(argc_address, argv_address);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	sleep(myid);
	srand((unsigned) time(NULL));

	for(i=0; i<ndim; i++)
	{
		random0 = rand()/(RAND_MAX+1.0);
		pos[0*ndim+i] = init[i];
//		printf("%lf\t", pos[0][i]);
	}

//	pos[0*ndim+0] = pos[0*ndim+0] + random0*1.0; 
//	pos[0*ndim+1] = pos[0*ndim+1] + random0*1.0;


	for(loop=0; loop<walks; loop++)
	{
		if(loop!=0)
		{
			for(i=0; i<ndim; i++)
			{	
				random1[i] = rand()/(RAND_MAX+1.0);
				random2[i] = rand()/(RAND_MAX+1.0);
				rand_delta[i] = sqrt(-2.0*pi*log(random1[i]))*cos(2.0*pi*random2[i])*mag[i];
				pos[loop*ndim+i] = pos[(loop-1)*ndim+i] + rand_delta[i];
//				printf("early%lf\t%lf\t%lf\n", rand_delta[i], pos[loop-1][i], pos[loop][i]);
			}
		}

		for(i=0; i<Theta_min; i++)
		{
			w_analysis[i]=0;
		}
		ng[loop] = get_w(&pos[loop*ndim], w_analysis, &halo_mass[loop], &beff[loop], &fsat[loop]);

//		printf("ng=%lf\n", ng[loop]);

		for(i=0; i<Theta_min; i++)
		{	
			printf("%.10lf\t", w_analysis[i]);
//			printf("\n");
			for(j=0; j<Theta_min; j++)
			{	
//				printf("temp[%d]=%lf\t", i, temp[i]);			
				temp[i] = temp[i]+(w_observe[j]-w_analysis[j])*covariance_inv[j*Theta_min+i];	
//				printf("covariance_inv[%d][%d]=%lf\n", j, i,covariance_inv[j][i]);
			}
			delta_w[loop] = delta_w[loop]+temp[i]*(w_observe[i]-w_analysis[i]);
//			printf("temp[%d]=%lf\t", i, temp[i]);
//			printf("w_observe[%d]=%lf\tw_analysis[%d]=%lf\n", i, w_observe[i], i, w_analysis[i]);
//			printf("delta_w[%d]=%lf\n", loop, delta_w[loop]);
		}
		ka_square[loop] = delta_w[loop] + (Number_density-ng[loop])*(Number_density-ng[loop])/(Sigma_ng*Sigma_ng);
//		printf("delta_w[%d]=%lf\n", loop, delta_w[loop]);
//		printf("ka_square[%d]=%lf\n", loop, ka_square[loop]);
			
		for(i=0 ;i<Theta_min; i++)
		{
			temp[i] = 0;
			w_analysis[i] = 0; 
		}
		for(i=0; i<ndim; i++)
		{
//			printf("prepos%lf\n",pos[loop][i]);
		}
//		printf("preng%lf\n", ng[loop]);
//		printf("preka%lf\n", ka_square[loop]);
		if(loop!=0)
		{
			judge = exp((ka_square[loop]-ka_square[loop-1])*(-0.5));
//			printf("judge = %le\n", judge);
			if((pos[loop*ndim+0]<10)||(pos[loop*ndim+1]>16)||(pos[loop*ndim+1]<10)||(pos[loop*ndim+1]>16))
			{
//				printf("fuck\n");
				for(i=0; i<ndim; i++)
				{
					pos[loop*ndim+i] = pos[(loop-1)*ndim+i];
				}
				ka_square[loop] = ka_square[loop-1];
				ng[loop]        = ng[loop-1];
				halo_mass[loop] = halo_mass[loop-1];
				beff[loop]      = beff[loop-1];
				fsat[loop]      = fsat[loop-1];
			}
			else if(judge<1)
			{
				random_xi = rand()/(RAND_MAX+1.0);
//				printf("random_xi = %lf\n", random_xi);
				if(random_xi>judge)
				{	
//					printf("fuckagain\n");
					for(i=0; i<ndim; i++)
					{
						pos[loop*ndim+i] = pos[(loop-1)*ndim+i];
					}
					ka_square[loop] = ka_square[loop-1];
					ng[loop]        = ng[loop-1];
					halo_mass[loop] = halo_mass[loop-1];
					beff[loop]      = beff[loop-1];
					fsat[loop]      = fsat[loop-1];					
				}
			}
		}
		for(j=0; j<ndim; j++)
		{
//			printf("pos[%d][%d] = %lf\t", loop, j, pos[loop*ndim+j]);
		}
//		printf("ng[%d] = %lf\tka_square[%d] = %lf\n", loop, ng[loop], loop, ka_square[loop]);

	}

	for(loop=0; loop<walks; loop++)
	{	
		for(j=0; j<ndim; j++)
		{
//			printf("pos[%d][%d] = %.8lf\t", loop, j, pos[loop*ndim+j]);
		}
//		printf("ng[%d] = %lf\tka_square[%d] = %lf\n", loop, ng[loop], loop, ka_square[loop]);
		result1[loop] = pos[loop*ndim+0];
		result2[loop] = pos[loop*ndim+1];
		result3[loop] = ng[loop];
		result4[loop] = halo_mass[loop];
		result5[loop] = beff[loop];
		result6[loop] = fsat[loop];
		result7[loop] = ka_square[loop];
	}


	MPI_Gather(result1, walks, MPI_DOUBLE, result_log_M_min,   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Gather(result2, walks, MPI_DOUBLE, result_M1_primer,   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result3, walks, MPI_DOUBLE, result_ng,          walks,   MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Gather(result4, walks, MPI_DOUBLE, result_halo_mass,   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result5, walks, MPI_DOUBLE, result_beff,        walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result6, walks, MPI_DOUBLE, result_fsat,        walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result7, walks, MPI_DOUBLE, result_ka,          walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 

	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0)
	{	
		foutput = fopen("foutput", "w");
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_log_M_min[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_M1_primer[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_ng[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", log10(result_halo_mass[i]));
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_beff[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_fsat[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_ka[i]);
		}
		fclose(foutput);
	}
	MPI_Finalize();

	free(temp);
	free(observe);
	free(w_observe);
	free(w_analysis);
	free(covariance_inv);
	free(pos);
	free(ng);
	free(delta_w);
	free(ka_square);
	free(result1);
	free(result2);
	free(result3);
	free(result4);
	free(result5);
	free(result6);
	free(result7);
	free(result_log_M_min);
	free(result_M1_primer);
	free(result_ng);
	free(result_halo_mass);
	free(result_beff);
	free(result_fsat);	
	free(result_ka);
}

int main(int argc, char **argv)
{
	double omegam;
	double rho;
	double *e_temp;
	char midfix[file_char_num];
	char **pre_file1,**pre_file2;
	struct Linear *li;
	struct Nonlinear *nonli;
	double init[ndim]={ 1.26514731e+01,   1.34501766e+01};
	int i, j, s, mum, num, pa, pe, loop;
	double start,end,cost;
	double xi_res[node]={0},wi_res[node]={0};	
	char **pl_file,**pnl_file;
	char lheader[header_number];
	char nlheader[header_number];
	char fallrr[file_char_num];

//	FILE *fv;
//	FILE *fcorr_test;
	FILE *parafile;
	FILE *redshiftfile;
	FILE *frr;
//	fl = fopen("spread_l_noinp_matterpower.dat", "r");
//	fnl = fopen("spread_noinp_matterpower.dat", "r");
	parafile = fopen("infile16","r");
	fscanf(parafile, "%s ", Prefix);
	fscanf(parafile, "%d ", &Redshift_number);	
	fscanf(parafile, "%s ", Linear_matter_power_spectrum);	
	fscanf(parafile, "%s",  Nonlinear_matter_power_spectrum);	
	fscanf(parafile, "%d",  &Start_point);
	fscanf(parafile, "%d",  &Theta_min);	
	fscanf(parafile, "%s",  W_mean);	
	fscanf(parafile, "%s",  Inverse_matrix);	
	fscanf(parafile, "%lf", &Number_density);	
	fscanf(parafile, "%lf", &Sigma_ng);
	fscanf(parafile, "%d",  &K_num);	
	fscanf(parafile, "%lf", &Omega_m0);	
	fscanf(parafile, "%lf", &Omega_lambda0);
	fscanf(parafile, "%s",  Redshift_file);	
	fscanf(parafile, "%lf", &y_low);
	fscanf(parafile, "%lf", &y_high);
	fscanf(parafile, "%s", 	fallrr);

	fclose(parafile);

	FILE *fl[Redshift_number],*fnl[Redshift_number];

//	Theta = (double *)malloc(Theta_min*sizeof(double));
	Theta = (double *)malloc(theta_max*sizeof(double));
	K = (double *)malloc(Redshift_number*K_num*sizeof(double));
	Pkl = (double *)malloc(Redshift_number*K_num*sizeof(double));
	Pk = (double *)malloc(Redshift_number*K_num*sizeof(double));
	li = (struct Linear *)malloc(Redshift_number*K_num*sizeof(struct Linear));
	nonli = (struct Nonlinear *)malloc(Redshift_number*K_num*sizeof(struct Nonlinear));
	Redshift = (double *)malloc(Redshift_number*sizeof(double));
	Redshift_distace = (double *)malloc(Redshift_number*sizeof(double));
	Redshift_distribution = (double *)malloc(Redshift_number*sizeof(double));	
	Redshift_weight = (double *)malloc(Redshift_number*sizeof(double));
	Dredshift_divide_dr = (double *)malloc(Redshift_number*sizeof(double));
	Sigma = (double *)malloc(Redshift_number*m_max*sizeof(double));
	Sigma_add = (double *)malloc(Redshift_number*(m_max+1)*sizeof(double));
	F = (double *)malloc(Redshift_number*K_num*sizeof(double));
	Dlninvsigma_divide_dm = (double *)malloc(Redshift_number*m_max*sizeof(double));
	Dn_divide_dm = (double *)malloc(Redshift_number*m_max*sizeof(double));
	Dn_divide_dm_eff = (double *)malloc(m_max*sizeof(double));
	Bias = (double *)malloc(Redshift_number*m_max*sizeof(double));
	Bias_eff = (double *)malloc(m_max*sizeof(double));	
	e_temp = (double *)malloc(Redshift_number*sizeof(double));
	Dr_divide_dz = (double *)malloc(Redshift_number*sizeof(double));	
	pre_file1 = (char **)malloc(Redshift_number*sizeof(char *));
	pre_file2 = (char **)malloc(Redshift_number*sizeof(char *));	
	pl_file = (char **)malloc(Redshift_number*sizeof(char *));
	pnl_file = (char **)malloc(Redshift_number*sizeof(char *));
	Corr_mm = (double *)malloc(Redshift_number*r_num*sizeof(double));
	gsl_interp_accel *acc1[Redshift_number];
	gsl_spline *spline1[Redshift_number];

	redshiftfile = fopen(Redshift_file,"r");
	for(i=0; i<Redshift_number; i++)
	{
		fscanf(redshiftfile, "%lf%lf%lf%lf", &Redshift[i], &Redshift_weight[i], &Redshift_distace[i], &Redshift_distribution[i]);		
		pre_file1[i] = (char *)malloc(file_char_num*sizeof(char));
		pre_file2[i] = (char *)malloc(file_char_num*sizeof(char));	
		pl_file[i]   = (char *)malloc(file_char_num*sizeof(char));
		pnl_file[i]  = (char *)malloc(file_char_num*sizeof(char));		
		strcpy(pre_file1[i],Prefix);
		strcpy(pre_file2[i],Prefix);		
//		printf("%s\n", pre_file[i]);	
//		printf("%lf\t%lf\t%lf\t%lf\n", Redshift[i], Redshift_weight[i], Redshift_distace[i], Redshift_distribution[i]);
	}
	fclose(redshiftfile);

/*
	printf("linear_matter_power_spectrum=%s\n",  linear_matter_power_spectrum);	
	printf("nonlinear_matter_power_spectrum=%s\n",  nonlinear_matter_power_spectrum);	
	printf("start_point=%d\n",  start_point);
	printf("Theta_min=%d\n",  Theta_min);	
	printf("w_mean=%s\n",  w_mean);	
	printf("inverse_matrix=%s\n",  inverse_matrix);	
	printf("number_density=%lf\n", number_density);	
	printf("sigma_ng=%lf\n", sigma_ng);	
	printf("z=%lf\n", z);
	printf("K_num=%d\n",  K_num);	
	printf("Omega_m0=%lf\n", Omega_m0);	
	printf("Omega_lambda0=%lf\n", Omega_lambda0);	
*/
	
	for(i=0; i<Redshift_number; i++)
	{
		sprintf(midfix, "%.6lf", Redshift[i]);
		strcpy(pl_file[i],strcat( strcat(pre_file1[i],Linear_matter_power_spectrum ),midfix) );
		strcpy(pnl_file[i],strcat( strcat(pre_file2[i],Nonlinear_matter_power_spectrum ),midfix) );		
//		printf("%s\n", pl_file[i]);		
//		printf("%s\n", pnl_file[i]);				
		fl[i]  = fopen( pl_file[i]  ,"r"); 
		if(fl[i] ==NULL)printf("cao");
		fnl[i] = fopen( pnl_file[i] ,"r");
		if(fnl[i] ==NULL)printf("cao");

		fgets(lheader,  header_number, fl[i]);
//		printf("%s\n", lheader);
		fgets(nlheader, header_number ,fnl[i]);
//		printf("%s\n", nlheader);

		for(j=0; j<K_num; j++)
		{
			fscanf(fl[i], "%lf%lf", &li[i*K_num+j].kl_r, &li[i*K_num+j].pkl_r);
			K[i*K_num+j] = li[i*K_num+j].kl_r;
			Pkl[i*K_num+j] = li[i*K_num+j].pkl_r;
			fscanf(fnl[i], "%lf%lf", &nonli[i*K_num+j].k_r, &nonli[i*K_num+j].pk_r);
			Pk[i*K_num+j] = nonli[i*K_num+j].pk_r;
//			printf("%14.4le\t%14.4le\t%14.4le\n", K[i*K_num+j], Pkl[i*K_num+j], Pk[i*K_num+j]);
		}
		fclose(fl[i]);
		fclose(fnl[i]);		
	}
/*
	for(i=0; i<Redshift_number; i++)
	{	
		fcorr_test = fopen("elucid_matterpower.dat", "r");
		for(j=0; j<K_num; j++)
		{
			fscanf(fcorr_test, "%lf%lf", &li[i*K_num+j].kl_r, &li[i*K_num+j].pkl_r);
			K[i*K_num+j] = li[i*K_num+j].kl_r;
			Pkl[i*K_num+j] = li[i*K_num+j].pkl_r;
			//fscanf(fcorr_test, "%lf%lf", &nonli[i*K_num+j].k_r, &nonli[i*K_num+j].pk_r);
			Pk[i*K_num+j] = Pkl[i*K_num+j];
			printf("%14.4le\t%14.4le\t%14.4le\n", K[i*K_num+j], Pkl[i*K_num+j], Pk[i*K_num+j]);
		}
	}
*/
	Y = (double *)malloc(y_node*sizeof(double));
	Yweight = (double *)malloc(y_node*sizeof(double));
	gsl_integration_glfixed_table *table_y = gsl_integration_glfixed_table_alloc(y_node);


	for(num=0; num<y_node; num++)
	{
		gsl_integration_glfixed_point(y_low,y_high,num,&Y[num],&Yweight[num],table_y);		
//		printf("%lf\t%lf\n", Y[num], Yweight[num]);
	}

	for(i=0; i<r_num; i++)
	{
		R_corr[i] = pow(10,(0.01*i -2.0));
	}

	Rh= Rho_m(0)/(1.9891E30/little_h)*pow((3.08568E22/little_h),3)*Omega_m0;       //(M_sun/h)/(Mpc/h)^3

	for(i=0; i<Redshift_number; i++)
	{
		e_temp[i] = sqrt(Omega_lambda0+Omega_m0*pow(1+Redshift[i],3));
		Dr_divide_dz[i] = H0/little_h*e_temp[i]/lightspeed;
	}

	gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(node);

	for(num=0; num<node; num++)
	{
		gsl_integration_glfixed_point(0,2*pi,num,&Xi[num],&Wi[num],table);
//		printf("%lf\t%lf\n", xi[num], wi[num]);
	}


	for(i=0; i<(m_max+1); i++)
	{
		M_add[i] = pow(10,(9.0+0.1*i));
		R_add[i] = pow(M_add[i]/Rh/(4.0/3*pi),(1.0/3));
//		printf("m_add[%d]=%le\tr_add[%d]=%lf\n", i, m_add[i], i, r_add[i]);	
	}
	Sigma_f(Sigma_add,R_add,m_max+1);

//	for(i=0; i<Redshift_number; i++)
//	{
//		for(j=0; j<(m_max+1); j++)
//		{
//			printf("Sigma_add[%d][%d]=%lf\n ", i, j, Sigma_add[i*(m_max+1)+j]);  
//		}
//	}
	
	for(i=0; i<m_max; i++)
	{
		M[i] = pow(10,(9.05+0.1*i));
		R[i] = pow((M[i]/Rh/(4.0/3*pi)),(1.0/3));
//		printf("m[%d]=%le\tr[%d]=%lf\n", i, m[i], i, r[i]);	
	}
	Sigma_f(Sigma,R,m_max);

//	for(i=0; i<Redshift_number; i++)
//	{
//		for(j=0; j<m_max; j++)
//		{
//			printf("Sigma[%d][%d]=%lf\n ", i, j, Sigma[i*m_max+j]);  
//		}
//	}
	
	Multiplicity(F);
	bias_f(Bias);
	for(i=0; i<Redshift_number; i++)
	{
		for(j=0; j<m_max; j++)
		{
    		Dlninvsigma_divide_dm[i*m_max+j] =(log(1.0/Sigma_add[i*(m_max+1)+j+1])-log(1.0/Sigma_add[i*(m_max+1)+j]))/(M_add[j+1]-M_add[j]);
    		Dn_divide_dm[i*m_max+j] = F[i*m_max+j]*Rh/M[j]*Dlninvsigma_divide_dm[i*m_max+j];
//			printf("F[%d][%d]=%lf\n",i, j, F[i*m_max+j]);
//			printf("Dlninvsigma_divide_dm[%d][%d] = %le\t Dn_divide_dm[%d][%d] = %le\n", i, j, Dlninvsigma_divide_dm[i*m_max+j], i, j, Dn_divide_dm[i*m_max+j]);
//			printf("Dn_divide_dm[%d][%d] = %le\n", i, j, Dn_divide_dm[i*m_max+j]);
//			printf("Bias[%d][%d] = %lf\n", i, j, Bias[i*m_max+j]);
		}
	}

	Norm = 0;
	for(i=0; i<Redshift_number; i++)
	{
		Norm = Norm + Redshift_distace[i]*Redshift_distace[i]*Dr_divide_dz[i]*Redshift_weight[i];	
	}
//	printf("Norm=%lf\n", Norm);

	for(j=0; j<m_max; j++)
	{
		Dn_divide_dm_eff[j]=0;
		Bias_eff[j]=0;
		for(i=0; i<Redshift_number; i++)
		{
			Dn_divide_dm_eff[j] = Dn_divide_dm_eff[j] + 
			Redshift_distace[i]*Redshift_distace[i]*Dr_divide_dz[i]*Dn_divide_dm[i*m_max+j]*Redshift_weight[i];    //Redshift_weight is normalized and Redshift_weight = z_distribution.
			Bias_eff[j] = Bias_eff[j]+Redshift_distace[i]*Redshift_distace[i]*Dr_divide_dz[i]*Bias[i*m_max+j]*Redshift_weight[i];    //Redshift_weight is normalized and Redshift_weight = z_distribution.
		}
		Dn_divide_dm_eff[j] = Dn_divide_dm_eff[j]/Norm;
		Bias_eff[j] = Bias_eff[j]/Norm;
//		printf("Dn_divide_dm_eff[%d]=%le\tBias_eff[%d]=%lf\n", j, Dn_divide_dm_eff[j], j, Bias_eff[j]);
	}

/*	for(i=0; i<Theta_min; i++)
	{
		Theta[i] = pow(10.0,0.2*i+0.2*Start_point-3.0)*pi/180.0;
		printf("%.10lf\n", Theta[i]);
	}
*/

	frr = fopen(fallrr, "r");
	for(i=0; i<theta_max; i++)
	{
		fscanf(frr, "%lld", &Rr_pair[i]); 
//		printf("%lld\n", rr_pair[i]);
	}
	fclose(frr);

	for(i=0; i<theta_max; i++)
	{
		Theta[i] = pow(10.0,0.2*i-3.0)*pi/180.0;
//		printf("%.10lf\n", Theta[i]);
	}

	for(i=0; i<Redshift_number*r_num; i++)
	{
		Corr_mm[i] = 0;
	}

	for(s=0; s<Redshift_number; s++)
	{	
		acc1[s]= gsl_interp_accel_alloc ();
		spline1[s]= gsl_spline_alloc (gsl_interp_cspline, K_num);
		gsl_spline_init (spline1[s], &K[s*K_num], &Pk[s*K_num], K_num);

		for(mum=0; mum<cut_number; mum++)
		{	
			for(j=r_cut[mum]; j<r_cut[mum+1]; j++)
			{
				for(i=0; i<r_cut_value[mum]; i++)
				{	
					for(num=0; num<node; num++)
					{
						xi_res[num] = Xi[num] + 2*pi*i;
						wi_res[num] = Wi[num];
						Corr_mm[s*r_num+j] = Corr_mm[s*r_num+j] + wi_res[num]*xi_res[num]/R_corr[j]/R_corr[j]/R_corr[j]*gsl_spline_eval(spline1[s], xi_res[num]/R_corr[j], acc1[s])*sin(xi_res[num]);
					}
				}
				Corr_mm[s*r_num+j] = Corr_mm[s*r_num+j]/2.0/pi/pi;
//				printf("%d\t%d\t%lf\n", s, j, Corr_mm[s*r_num+j]);
			}
		}
		gsl_spline_free (spline1[s]);
		gsl_interp_accel_free (acc1[s]);
	}
//---------------------------------------------------start MCMC----------------------------------------------------------

	start =clock();

	getMCMC(init, &argc, &argv);

	end =clock();
	cost = (end -start)/CLOCKS_PER_SEC ;
	printf("time=%lf\n",cost);

	return 0;
}
