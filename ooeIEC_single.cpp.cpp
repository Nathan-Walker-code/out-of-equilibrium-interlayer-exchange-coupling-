#include<iostream>
#include<cmath>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<complex>
#include<iomanip>
#include <nag.h>
#include<chrono>

using namespace std;
using namespace Eigen;
using namespace chrono;


/* details of mathematics in- "Bias controlled Interlayer Exchange Coupling", N. Walker, A. Durie, A. Umerski, arXiv:2604.05705*/


#ifdef __cplusplus
extern "C"
{
#endif
static void NAG_CALL f(const double x[], Integer nx, double fv[],
                       Integer *iflag, Nag_Comm *comm);
#ifdef __cplusplus
}
#endif

double x,y,theta;

double right_fermi=0.0;                                     /*right fermi level- VARIABLE*/
   
double voltage=-0.04;                                     /*bias - VARIABLE*/
 
double Temp=316;                                              /*temperature in kelvin- VARIABLE*/        
double KbT=Temp*(0.00008617/13.6057);                             /*Boltz.*temperature*/

const int N=2;                                             /*Number of orbitals- VARIABLE*/

double grid_size=55;                                        /*grid density for STT term- VARIABLE*/

double delta=0.0000000001;                                   /*imaginary part of energy- VARIABLE*/


/*NOTE: Need cleavage to be in a non-magnetic layer (FM spacer!!)  and have a plane to the left of it that has the required hopping*/
/*NOTE: SETUP FOR SPACER 2!!*/


const int left_insulator_thickness=4;                     /*thickness of insulating layer- VARIABLE*/



const int FM_l_num=5;                                    /*left FM thickness- VARIABLE*/
const int S_num=2;                                         /*spacer thickness NEEDS TO BE EVEN- FIXED IN THIS CODE!!!*/
const int FM_r_num=5;                                              /*right FM thickness- VARIABLE*/

const int num_spacer_eqm=5;                          /*Number of spacer thicknesses for eqm (even starting from 2)- VARIABLE*/
int num_spacer;                                     /*number of loops for ooe part of code*/


double arr[num_spacer_eqm];                         /*stores integrand values for eqm (multiple spacers) at specific E,x,y in f*/
double left_eqm[num_spacer_eqm];                        /*stores final left eqm IEC values (multiple spacers)*/
double right_eqm[num_spacer_eqm];                         /*stores final right eqm IEC values (multiple spacers)*/
double holder[num_spacer_eqm];                         /*stores final ooe IEC values (multiple spacers)*/

double W_B=1.0;                           /*band gap width parameter- VARIABLE*/
double W_FM=1.0;                              /*HG width parameter- VARIABLE*/

int big_flag;                               /*determines whether f is being called for eqm or ooe*/
 



double f(complex<double> E)
{

/*This function evalutes the spin current densities for the STT and IEC terms as well as carry out the STT theta integration*/



  int total_thickness=1+left_insulator_thickness+FM_l_num+S_num+FM_r_num+1;                                                                                                                                            /*total thickness ie. all layers plus two- VARIABLE*/

  int cleavage=left_insulator_thickness+FM_l_num+(0.5*S_num);                                                                                                                                      /*location of cleavage ie. index of plane to left of cleavage- VARIABLE*/

  
  complex<double> w1(2*cos(x)+2*cos(y),0);          /*Structure factor for SC- VARIABLE*/
  complex<double> w2=E;                    

  MatrixXcd I(N,N);  
  complex<double> I1(1,0),I2(0,0),I3(0,0),I4(1,0);                                /*Identity- VARIABLE*/
  I<<I1,I2,I3,I4;                                                 /*VARIABLE*/
  MatrixXcd Z(N,N);
  complex<double> Z1(0,0),Z2(0,0),Z3(0,0),Z4(0,0);                               /*Zero matrix- VARIABLE*/
  Z<<Z1,Z2,Z3,Z4;                                                  /*VARIABLE*/

  ComplexEigenSolver<MatrixXcd> ces;

  
                                  /*declare all the on-sites and hopping elements*/
  
  complex<double> u_l_lead_1(-0.600624+voltage,0),u_l_lead_2(0,0),u_l_lead_3(0,0),u_l_lead_4(0.513965+voltage,0);                                                                                      /*on-site of left lead- VARIABLE*/
  complex<double> u_l_ins_1(-0.22371,0),u_l_ins_2(0,0),u_l_ins_3(0,0),u_l_ins_4(0.401852,0);                                                                                                              /*on-site for left insulator- VARIABLE*/
  complex<double> u_r_lead_1(-0.600624,0),u_r_lead_2(0,0),u_r_lead_3(0,0),u_r_lead_4(0.513965,0);                                                                                                               /*on-site of left lead- VARIABLE*/
  complex<double> u_spacer_1(-0.600624,0),u_spacer_2(0,0),u_spacer_3(0,0),u_spacer_4(0.513965,0);                                                                                          /*on-site of insulator spacer- VARIABLE*/
  complex<double> u_l_FM_up_1(-0.339172,0),u_l_FM_up_2(0,0),u_l_FM_up_3(0,0),u_l_FM_up_4(0.653281,0);                                                                                             /*on-site for left up spin FM- VARIABLE*/
  complex<double> u_l_FM_down_1(-0.22371,0),u_l_FM_down_2(0,0),u_l_FM_down_3(0,0),u_l_FM_down_4(0.401852,0);                                                                                     /*on-site for left down spin FM- VARIABLE*/
  complex<double> u_r_FM_up_1(-0.339172,0),u_r_FM_up_2(0,0),u_r_FM_up_3(0,0),u_r_FM_up_4(0.653281,0);                                                                                                     /*on-site for right up spin FM- VARIABLE*/
  complex<double> u_r_FM_down_1(-0.22371,0),u_r_FM_down_2(0,0),u_r_FM_down_3(0,0),u_r_FM_down_4(0.401852,0);                                                                                             /*on-site for right down spin FM- VARIABLE*/
  

  
  

  complex<double> t_l_lead_1(0.0620793,0),t_l_lead_2(0.0244402,0),t_l_lead_3(0.0244402,0),t_l_lead_4(-0.160936,0);                                                                                    /*hopping of left lead- VARIABLE*/
  complex<double> t_l_ins_1(0.0497084,0),t_l_ins_2(0.0101115*W_B,0),t_l_ins_3(0.0101115*W_B,0),t_l_ins_4(-0.0866541,0);                                                                               /*hopping for insulator- VARIABLE*/
  complex<double> t_r_lead_1(0.0620793,0),t_r_lead_2(0.0244402,0),t_r_lead_3(0.0244402,0),t_r_lead_4(-0.160936,0);                                                                         /*hopping of left lead- VARIABLE*/
  complex<double> t_spacer_1(0.0620793,0),t_spacer_2(0.0244402,0),t_spacer_3(0.0244402,0),t_spacer_4(-0.160936,0);                                                                         /*on-site of spacer- VARIABLE*/
  complex<double> t_l_FM_up_1(-0.00751873,0),t_l_FM_up_2(-0.0409919,0),t_l_FM_up_3(-0.0409919,0),t_l_FM_up_4(-0.181094,0);                                                                               /*hopping for left up spin FM- VARIABLE*/
  complex<double> t_l_FM_down_1(0.0497084,0),t_l_FM_down_2(0.0101115*W_FM,0),t_l_FM_down_3(0.0101115*W_FM,0),t_l_FM_down_4(-0.0866541,0);                                                                          /*hopping for left down spin FM- VARIABLE*/
  complex<double> t_r_FM_up_1(-0.00751873,0),t_r_FM_up_2(-0.0409919,0),t_r_FM_up_3(-0.0409919,0),t_r_FM_up_4(-0.181094,0);                                                                           /*hopping for right up spin FM- VARIABLE*/
  complex<double> t_r_FM_down_1(0.0497084,0),t_r_FM_down_2(0.0101115*W_FM,0),t_r_FM_down_3(0.0101115*W_FM,0),t_r_FM_down_4(-0.0866541,0);                                                                       /*hopping for right down spin FM- VARIABLE*/
 
  

                              /*declare and construct all on-site, hopping and v matrices from elements*/

  MatrixXcd u_l_lead(N,N);
  u_l_lead<<u_l_lead_1,u_l_lead_2,u_l_lead_3,u_l_lead_4;                                          /*VARIABLE*/
  MatrixXcd u_r_lead(N,N);
  u_r_lead<<u_r_lead_1,u_r_lead_2,u_r_lead_3,u_r_lead_4;                                          /*VARIABLE*/
  MatrixXcd u_l_ins(N,N);
  u_l_ins<<u_l_ins_1,u_l_ins_2,u_l_ins_3,u_l_ins_4;                                                  /*VARIABLE*/
  MatrixXcd u_spacer(N,N);
  u_spacer<<u_spacer_1,u_spacer_2,u_spacer_3,u_spacer_4;                              /*VARIABLE*/
  MatrixXcd u_l_FM_up(N,N);
  u_l_FM_up<<u_l_FM_up_1,u_l_FM_up_2,u_l_FM_up_3,u_l_FM_up_4;                                                     /*VARIABLE*/
  MatrixXcd u_l_FM_down(N,N);
  u_l_FM_down<<u_l_FM_down_1,u_l_FM_down_2,u_l_FM_down_3,u_l_FM_down_4;                                             /*VARIABLE*/
  MatrixXcd u_r_FM_up(N,N);
  u_r_FM_up<<u_r_FM_up_1,u_r_FM_up_2,u_r_FM_up_3,u_r_FM_up_4;                                                     /*VARIABLE*/
  MatrixXcd u_r_FM_down(N,N);
  u_r_FM_down<<u_r_FM_down_1,u_r_FM_down_2,u_r_FM_down_3,u_r_FM_down_4;                                            /*VARIABLE*/
  

  MatrixXcd t_l_lead(N,N);
  t_l_lead<<t_l_lead_1,t_l_lead_2,t_l_lead_3,t_l_lead_4;                                       /*VARIABLE*/
  MatrixXcd t_r_lead(N,N);
  t_r_lead<<t_r_lead_1,t_r_lead_2,t_r_lead_3,t_r_lead_4;                                      /*VARIABLE*/
  MatrixXcd t_l_ins(N,N);
  t_l_ins<<t_l_ins_1,t_l_ins_2,t_l_ins_3,t_l_ins_4;                                                          /*VARIABLE*/
  MatrixXcd t_spacer(N,N);
  t_spacer<<t_spacer_1,t_spacer_2,t_spacer_3,t_spacer_4;                              /*VARIABLE*/
  MatrixXcd t_l_FM_up(N,N);
  t_l_FM_up<<t_l_FM_up_1,t_l_FM_up_2,t_l_FM_up_3,t_l_FM_up_4;                                                     /*VARIABLE*/
  MatrixXcd t_l_FM_down(N,N);
  t_l_FM_down<<t_l_FM_down_1,t_l_FM_down_2,t_l_FM_down_3,t_l_FM_down_4;                                        /*VARIABLE*/
  MatrixXcd t_r_FM_up(N,N);
  t_r_FM_up<<t_r_FM_up_1,t_r_FM_up_2,t_r_FM_up_3,t_r_FM_up_4;                                                     /*VARIABLE*/
  MatrixXcd t_r_FM_down(N,N);
  t_r_FM_down<<t_r_FM_down_1,t_r_FM_down_2,t_r_FM_down_3,t_r_FM_down_4;                                           /*VARIABLE*/
 
  

  MatrixXcd v_l_lead(N,N);
  v_l_lead=(w2*I)-u_l_lead-(w1*t_l_lead);
  MatrixXcd v_r_lead(N,N);
  v_r_lead=(w2*I)-u_r_lead-(w1*t_r_lead);
  MatrixXcd v_spacer(N,N);
  v_spacer=(w2*I)-u_spacer-(w1*t_spacer);
  MatrixXcd v_l_FM_up(N,N);
  v_l_FM_up=(w2*I)-u_l_FM_up-(w1*t_l_FM_up);
  MatrixXcd v_l_FM_down(N,N);
  v_l_FM_down=(w2*I)-u_l_FM_down-(w1*t_l_FM_down);
  MatrixXcd v_r_FM_up(N,N);
  v_r_FM_up=(w2*I)-u_r_FM_up-(w1*t_r_FM_up);
  MatrixXcd v_r_FM_down(N,N);
  v_r_FM_down=(w2*I)-u_r_FM_down-(w1*t_r_FM_down);
  


/*Barrier potential drop modelled as linear*/


                               /*Start of shifted left insulator potential*/

  MatrixXcd shifted_v_l_ins[left_insulator_thickness];

  if(left_insulator_thickness==1)
    {
      shifted_v_l_ins[0]=(w2*I)-u_l_ins-(w1*t_l_ins)-(voltage*I);
    }
  else
    {
      complex<double> increment_1(0.5*voltage/(left_insulator_thickness-1),0);

      for(int i=0;i<left_insulator_thickness;i++)
	{
	  complex<double> number_1(i,0);
	  shifted_v_l_ins[i]=(w2*I)-u_l_ins-(w1*t_l_ins)-(number_1*increment_1*I)-(0.5*voltage*I);
	}
    }
  

                                     /*End of shifted left insulator potential*/


                              

  /*store all onsites and hoppings in correct order in arrays for adlayering*/


                                  /*construct the system arrays- VARIABLE setup for specific system type!!!*/

  
 /*spin up v's*/
  MatrixXcd paras_v_up[total_thickness];

  paras_v_up[0]=v_l_lead;
  for(int i=1;i<=left_insulator_thickness;i++)
    {
      paras_v_up[i]=shifted_v_l_ins[left_insulator_thickness-i];
    }
  for(int i=left_insulator_thickness+1;i<=left_insulator_thickness+FM_l_num;i++)
    {
      paras_v_up[i]=v_l_FM_up;
    }
  for(int i=left_insulator_thickness+FM_l_num+1;i<=left_insulator_thickness+FM_l_num+S_num;i++)
    {
      paras_v_up[i]=v_spacer;
    }
  for(int i=left_insulator_thickness+FM_l_num+S_num+1;i<=left_insulator_thickness+FM_l_num+S_num+FM_r_num;i++)
    {
      paras_v_up[i]=v_r_FM_up;
    }
  paras_v_up[total_thickness-1]=v_r_lead;





  
  
  /*spin down v's*/
  MatrixXcd paras_v_down[total_thickness];
  
  paras_v_down[0]=v_l_lead;
  for(int i=1;i<=left_insulator_thickness;i++)
    {
      paras_v_down[i]=shifted_v_l_ins[left_insulator_thickness-i];
    }
  for(int i=left_insulator_thickness+1;i<=left_insulator_thickness+FM_l_num;i++)
    {
      paras_v_down[i]=v_l_FM_down;
    }
  for(int i=left_insulator_thickness+FM_l_num+1;i<=left_insulator_thickness+FM_l_num+S_num;i++)
    {
      paras_v_down[i]=v_spacer;
    }
  for(int i=left_insulator_thickness+FM_l_num+S_num+1;i<=left_insulator_thickness+FM_l_num+S_num+FM_r_num;i++)
    {
      paras_v_down[i]=v_r_FM_down;
    }
  paras_v_down[total_thickness-1]=v_r_lead;
  


  
  
  /*spin up hoppings*/
  MatrixXcd paras_t_up[total_thickness];

  paras_t_up[0]=t_l_lead;
  for(int i=1;i<=left_insulator_thickness;i++)
    {
      paras_t_up[i]=t_l_ins;
    }
  for(int i=left_insulator_thickness+1;i<=left_insulator_thickness+FM_l_num;i++)
    {
      paras_t_up[i]=t_l_FM_up;
    }
  for(int i=left_insulator_thickness+FM_l_num+1;i<=left_insulator_thickness+FM_l_num+S_num;i++)
    {
      paras_t_up[i]=t_spacer;
    }
  for(int i=left_insulator_thickness+FM_l_num+S_num+1;i<=left_insulator_thickness+FM_l_num+S_num+FM_r_num;i++)
    {
      paras_t_up[i]=t_r_FM_up;
    }
  paras_t_up[total_thickness-1]=t_r_lead;





    
  /*spin down hoppings*/
  MatrixXcd paras_t_down[total_thickness];

  paras_t_down[0]=t_l_lead;
  for(int i=1;i<=left_insulator_thickness;i++)
    {
      paras_t_down[i]=t_l_ins;
    }
  for(int i=left_insulator_thickness+1;i<=left_insulator_thickness+FM_l_num;i++)
    {
      paras_t_down[i]=t_l_FM_down;
    }
  for(int i=left_insulator_thickness+FM_l_num+1;i<=left_insulator_thickness+FM_l_num+S_num;i++)
    {
      paras_t_down[i]=t_spacer;
    }
  for(int i=left_insulator_thickness+FM_l_num+S_num+1;i<=left_insulator_thickness+FM_l_num+S_num+FM_r_num;i++)
    {
      paras_t_down[i]=t_r_FM_down;
    }
  paras_t_down[total_thickness-1]=t_r_lead;
  
    
  
                                                         /*end of system arrays*/


  
  
/*SGF's of semi-infnite leads evaluated using mobius transformation method- "Closed-form solutions to surface Green's funtions", A. Umerski, Phys. Rev. B, 1997*/


  


                                                     /*start of left lead SGF*/
                                               /*NOTE- Wont work if LEAD an FM*/
  

  MatrixXcd X_l_lead(2*N,2*N);          /*declares MT matrix using up spin as lead non-FM*/
  X_l_lead<<Z,paras_t_up[0].inverse(),(-1)*paras_t_up[0].adjoint(),paras_v_up[0]*paras_t_up[0].inverse();

  double ev_m_l_lead[2*N];

  int indx_l_lead[2*N];

  complex<double> ev_sorted_l_lead[2*N];

  MatrixXcd L_l_lead(2*N,2*N),O_l_lead(2*N,2*N);

  MatrixXcd L2_l_lead(N,N),O2_l_lead(N,N);

  MatrixXcd g_l_lead(N,N);                                                   /*SGF of left lead*/

  ces.compute(X_l_lead);                                       

  for(int i=0;i<2*N;i++)
    {
      ev_m_l_lead[i]=abs(ces.eigenvalues()[i]);
    }

  for(int i=0;i<2*N;i++)
    {
      indx_l_lead[i]=i;
    }

  int length=sizeof(indx_l_lead)/sizeof(indx_l_lead[0]);

  sort(indx_l_lead, indx_l_lead+length,[&](int i, int j){return ev_m_l_lead[i]<ev_m_l_lead[j];});           /*sort e.values of X by increasing modulus*/

  for(int i=0;i<2*N;i++)
    {
      ev_sorted_l_lead[i]=ces.eigenvalues()[indx_l_lead[i]];
    }


	/*form e.value and e.vector matrices for final formula*/

  L_l_lead<<ev_sorted_l_lead[0],0,0,0,0,ev_sorted_l_lead[1],0,0,0,0,ev_sorted_l_lead[2],0,0,0,0,ev_sorted_l_lead[3];                                                                                                   /*VARIABLE*/

  O_l_lead<<ces.eigenvectors().col(indx_l_lead[0]),ces.eigenvectors().col(indx_l_lead[1]),ces.eigenvectors().col(indx_l_lead[2]),ces.eigenvectors().col(indx_l_lead[3]);                         /*VARIABLE*/

  L2_l_lead=L_l_lead.block<N,N>(N,N);
  O2_l_lead=O_l_lead.block<N,N>(0,N);

  g_l_lead=O2_l_lead*L2_l_lead.inverse()*O2_l_lead.inverse()*paras_t_up[0].inverse();

  MatrixXcd g_l_check(N,N);                                                                      /*SGF check*/

  g_l_check=paras_t_up[0].inverse()*(((-1)*paras_t_up[0].adjoint()*g_l_lead+paras_v_up[0]*paras_t_up[0].inverse()).inverse());

  complex<double> sum_1(0,0);
  complex<double> sum_2(0,0);

  for(int i=0;i<N;i++)
    {
      for(int j=0;j<N;j++)
	{
	  sum_1+=g_l_lead(i,j);
	  sum_2+=g_l_check(i,j);
	}
    }

  double check_1=abs(real(sum_1-sum_2));

  if(check_1>0.00000001)
    {
      cout<<"WARNING! SGF of left lead may not have converged at k_x= "<<x<<", k_y= "<<y<<", E= "<<E<<"."<<endl;
    }

  
                                                       /*end of left lead SGF*/
  
 
                                                 

  
  
  
                                                       /*start of right lead SGF*/
                                                    /*NOTE- wont work if lead an FM*/
 

  MatrixXcd X_r_lead(2*N,2*N);
  X_r_lead<<Z,(paras_t_up[total_thickness-1].adjoint()).inverse(),(-1)*paras_t_up[total_thickness-1],paras_v_up[total_thickness-1]*(paras_t_up[total_thickness-1].adjoint()).inverse();

  double ev_m_r_lead[2*N];

  int indx_r_lead[2*N];

  complex<double> ev_sorted_r_lead[2*N];

  MatrixXcd L_r_lead(2*N,2*N),O_r_lead(2*N,2*N);

  MatrixXcd L2_r_lead(N,N),O2_r_lead(N,N);

  MatrixXcd g_r_lead(N,N);                                                   /*SGF of left lead*/

  ces.compute(X_r_lead);                                       

  for(int i=0;i<2*N;i++)
    {
      ev_m_r_lead[i]=abs(ces.eigenvalues()[i]);
    }

  for(int i=0;i<2*N;i++)
    {
      indx_r_lead[i]=i;
    }

  sort(indx_r_lead, indx_r_lead+length,[&](int i, int j){return ev_m_r_lead[i]<ev_m_r_lead[j];});

  for(int i=0;i<2*N;i++)
    {
      ev_sorted_r_lead[i]=ces.eigenvalues()[indx_r_lead[i]];
    }

  L_r_lead<<ev_sorted_r_lead[0],0,0,0,0,ev_sorted_r_lead[1],0,0,0,0,ev_sorted_r_lead[2],0,0,0,0,ev_sorted_r_lead[3];                                                                              /*VARIABLE*/

  O_r_lead<<ces.eigenvectors().col(indx_r_lead[0]),ces.eigenvectors().col(indx_r_lead[1]),ces.eigenvectors().col(indx_r_lead[2]),ces.eigenvectors().col(indx_r_lead[3]);                           /*VARIABLE*/

  L2_r_lead=L_r_lead.block<N,N>(N,N);
  O2_r_lead=O_r_lead.block<N,N>(0,N);

  g_r_lead=O2_r_lead*L2_r_lead.inverse()*O2_r_lead.inverse()*(paras_t_up[total_thickness-1].adjoint()).inverse();

  MatrixXcd g_r_check(N,N);                                                                      /*SGF check*/

  g_r_check=(paras_t_up[total_thickness-1].adjoint()).inverse()*(((-1)*paras_t_up[total_thickness-1]*g_r_lead+paras_v_up[total_thickness-1]*(paras_t_up[total_thickness-1].adjoint()).inverse()).inverse()); 

  complex<double> sum_3(0,0);
  complex<double> sum_4(0,0);

  for(int i=0;i<N;i++)
    {
      for(int j=0;j<N;j++)
	{
	  sum_3+=g_r_lead(i,j);
	  sum_4+=g_r_check(i,j);
	}
    }

  double check_2=abs(real(sum_3-sum_4));

  if(check_2>0.00000001)
    {
      cout<<"WARNING! SGF of right lead may not have converged at k_x= "<<x<<", k_y= "<<y<<", E= "<<E<<"."<<endl;
    }


                                                       /*end of right lead SGF*/


/*adlayers up to cleavage*/
  
  
                                                       /*start of adlayering*/

  MatrixXcd g_l_up(N,N),g_l_down(N,N),g_r_up(N,N),g_r_down(N,N);

  g_l_up=g_l_lead;
  g_l_down=g_l_lead;
  g_r_up=g_r_lead;
  g_r_down=g_r_lead;

  MatrixXcd t_av(N,N);
  complex<double> two(2,0);

  if(cleavage==0)
    {
    }
  else
    {
      for(int i=1;i<=cleavage;i++)                   /*Left up spin adlayering*/
	{     
	  for(int j=0;j<N;j++)                   /*Geometric average of hoppings*/
	    {
	      for(int k=0;k<N;k++)
		{
		  
		  if(real(paras_t_up[i-1](j,k))>0 && real(paras_t_up[i](j,k))>0)
		    {
		      t_av(j,k)=sqrt(paras_t_up[i-1](j,k)*paras_t_up[i](j,k));
		    }
		  else if(real(paras_t_up[i-1](j,k))<0 && real(paras_t_up[i](j,k))<0)
		    {
		      t_av(j,k)=-sqrt(paras_t_up[i-1](j,k)*paras_t_up[i](j,k));
		    }
		  else
		    {
		      t_av(j,k)=(paras_t_up[i-1](j,k)+paras_t_up[i](j,k))/two;
		    }	  
		}
	    }
	  
	  g_l_up=t_av.inverse()*(((-1)*t_av.adjoint()*g_l_up+paras_v_up[i]*t_av.inverse()).inverse());
	}
    }
  

  if(cleavage==0)
    {
    }
  else
    {
      for(int i=1;i<=cleavage;i++)                /*left down spin adlayering*/
	{
	  for(int j=0;j<N;j++)                   /*Geometric average of hoppings*/
	    {
	      for(int k=0;k<N;k++)
		{
		  
		  if(real(paras_t_down[i-1](j,k))>0 && real(paras_t_down[i](j,k))>0)
		    {
		      t_av(j,k)=sqrt(paras_t_down[i-1](j,k)*paras_t_down[i](j,k));
		    }
		  else if(real(paras_t_down[i-1](j,k))<0 && real(paras_t_down[i](j,k))<0)
		    {
		      t_av(j,k)=-sqrt(paras_t_down[i-1](j,k)*paras_t_down[i](j,k));
		    }
		  else
		    {
		      t_av(j,k)=(paras_t_down[i-1](j,k)+paras_t_down[i](j,k))/two;
		    }	  
		}
	    }
	  
	  g_l_down=t_av.inverse()*(((-1)*t_av.adjoint()*g_l_down+paras_v_down[i]*t_av.inverse()).inverse());
	}
    }
  

  if(cleavage==total_thickness-2)
    {
    }
  else
    {
      for(int i=total_thickness-2;i>=cleavage+1;i=i-1)             /*right up spin adlayering*/
	{
	  for(int j=0;j<N;j++)                   /*Geometric average of hoppings*/
	    {
	      for(int k=0;k<N;k++)
		{
		  
		  if(real(paras_t_up[i+1](j,k))>0 && real(paras_t_up[i](j,k))>0)
		    {
		      t_av(j,k)=sqrt(paras_t_up[i+1](j,k)*paras_t_up[i](j,k));
		    }
		  else if(real(paras_t_up[i+1](j,k))<0 && real(paras_t_up[i](j,k))<0)
		    {
		      t_av(j,k)=-sqrt(paras_t_up[i+1](j,k)*paras_t_up[i](j,k));
		    }
		  else
		    {
		      t_av(j,k)=(paras_t_up[i+1](j,k)+paras_t_up[i](j,k))/two;
		    }	  
		}
	    }
	  
	  g_r_up=(t_av.adjoint()).inverse()*(((-1)*t_av*g_r_up+paras_v_up[i]*(t_av.adjoint()).inverse()).inverse()); 
	}
    }
  

  if(cleavage==total_thickness-2)
    {
    }
  else
    {
      for(int i=total_thickness-2;i>=cleavage+1;i=i-1)                /*right down spin adlayering*/
	{
	  for(int j=0;j<N;j++)                   /*Geometric average of hoppings*/
	    {
	      for(int k=0;k<N;k++)
		{
		  
		  if(real(paras_t_down[i+1](j,k))>0 && real(paras_t_down[i](j,k))>0)
		    {
		      t_av(j,k)=sqrt(paras_t_down[i+1](j,k)*paras_t_down[i](j,k));
		    }
		  else if(real(paras_t_down[i+1](j,k))<0 && real(paras_t_down[i](j,k))<0)
		    {
		      t_av(j,k)=-sqrt(paras_t_down[i+1](j,k)*paras_t_down[i](j,k));
		    }
		  else
		    {
		      t_av(j,k)=(paras_t_down[i+1](j,k)+paras_t_down[i](j,k))/two;
		    }	  
		}
	    }
	  
	  g_r_down=(t_av.adjoint()).inverse()*(((-1)*t_av*g_r_down+paras_v_down[i]*(t_av.adjoint()).inverse()).inverse()); 
	}
    }

                                                         /*end of adlayering*/


  

  if(big_flag==0)               /*does eqm integrand*/
    {
      
      MatrixXcd array_l_up[num_spacer_eqm];             /*find the SGF's for each spacer thickness and store all for efficiency*/
      MatrixXcd array_l_down[num_spacer_eqm];
      MatrixXcd array_r_up[num_spacer_eqm];
      MatrixXcd array_r_down[num_spacer_eqm];
      
      array_l_up[0]=g_l_up;
      array_l_down[0]=g_l_down;
      array_r_up[0]=g_r_up;
      array_r_down[0]=g_r_down;
      
      for(int i=1;i<num_spacer_eqm;i++)
	{
	  array_l_up[i]=t_spacer.inverse()*(((-1)*t_spacer.adjoint()*array_l_up[i-1]+v_spacer*t_spacer.inverse()).inverse());
	}
      
      for(int i=1;i<num_spacer_eqm;i++)
	{
	  array_l_down[i]=t_spacer.inverse()*(((-1)*t_spacer.adjoint()*array_l_down[i-1]+v_spacer*t_spacer.inverse()).inverse());
	}
      
      for(int i=1;i<num_spacer_eqm;i++)
	{
	  array_r_up[i]=(t_spacer.adjoint()).inverse()*(((-1)*t_spacer*array_r_up[i-1]+v_spacer*(t_spacer.adjoint()).inverse()).inverse()); 
	}
      
      for(int i=1;i<num_spacer_eqm;i++)
	{
	  array_r_down[i]=(t_spacer.adjoint()).inverse()*(((-1)*t_spacer*array_r_down[i-1]+v_spacer*(t_spacer.adjoint()).inverse()).inverse()); 
	}
      
      
      
      
      
      
      MatrixXcd R_uu(N,N),R_dd(N,N),R_ud(N,N),R_du(N,N),R(N,N);
      
      for(int i=0;i<num_spacer_eqm;i++)
	{
	  R_uu=I-(array_r_up[i]*paras_t_up[cleavage]*array_l_up[i]*paras_t_up[cleavage].adjoint());
	  R_dd=I-(array_r_down[i]*paras_t_up[cleavage]*array_l_down[i]*paras_t_up[cleavage].adjoint());
	  R_ud=I-(array_r_up[i]*paras_t_up[cleavage]*array_l_down[i]*paras_t_up[cleavage].adjoint());
	  R_du=I-(array_r_down[i]*paras_t_up[cleavage]*array_l_up[i]*paras_t_up[cleavage].adjoint());
	  
	  R=R_uu*R_dd*R_ud.inverse()*R_du.inverse();
	  
	  arr[i]= 0.5*(real(log(R.determinant()))*pow(M_PI,-1)*(-KbT/(2*M_PI))); 
	}
      
      return 0;
    }
  else                  /*does ooe integrand*/
    {

      
      if(num_spacer==0)                           /*cant do all at the same time due to AI choosing unpredictable energies*/
	{
    }
  else
    {
      for(int i=1;i<=num_spacer;i++)
	{
	  g_l_up=t_spacer.inverse()*(((-1)*t_spacer.adjoint()*g_l_up+v_spacer*t_spacer.inverse()).inverse());
	  g_l_down=t_spacer.inverse()*(((-1)*t_spacer.adjoint()*g_l_down+v_spacer*t_spacer.inverse()).inverse());
	  g_r_up=(t_spacer.adjoint()).inverse()*(((-1)*t_spacer*g_r_up+v_spacer*(t_spacer.adjoint()).inverse()).inverse());
	  g_r_down=(t_spacer.adjoint()).inverse()*(((-1)*t_spacer*g_r_down+v_spacer*(t_spacer.adjoint()).inverse()).inverse());
	}
    }
  
  
  
      
	    

  complex<double> FL_L(right_fermi+voltage,0);                                                  /*Left fermi level- VARIABLE*/
  complex<double> FL_R(right_fermi,0);                                                    /*Right fermi level- VARIABLE*/
  

  MatrixXcd A(2*N,2*N),B(2*N,2*N),t_matrix(2*N,2*N),g_r_theta(2*N,2*N),g_l_matrix(2*N,2*N),g_r_matrix(2*N,2*N),T(2*N,2*N),Pauli(2*N,2*N);                                  /*Keldysh/Landauer stuff*/
  
  
  MatrixXcd IL(2*N,2*N);
  complex<double> IL1(1,0),IL2(0,0),IL3(0,0),IL4(0,0),IL5(0,0),IL6(1,0),IL7(0,0),IL8(0,0),IL9(0,0),IL10(0,0),IL11(1,0),IL12(0,0),IL13(0,0),IL14(0,0),IL15(0,0),IL16(1,0);                      /*Larger Identity- VARIABLE*/
  IL<<IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12,IL13,IL14,IL15,IL16;
  
  complex<double> z1(1,0);
  
  complex<double> f1,f2;                                                  /*fermi functions*/
  f1=z1/(z1+pow(M_E,(w2-FL_L)/KbT));
  f2=z1/(z1+pow(M_E,(w2-FL_R)/KbT));
  
  complex<double> Pauli_1(0,0),Pauli_2(0,-1),Pauli_3(0,1),Pauli_4(0,0);              /*Pauli y*/
  Pauli<<Pauli_1*I,Pauli_2*I,Pauli_3*I,Pauli_4*I;

  


  double s_theta=0.0;
  double n_tot_theta=0.0;

    for(theta=0;theta<=M_PI;theta=theta+M_PI/10)                      /*start of theta integration*/
	{
  
	  complex<double> T_1(cos(theta/2.0),0),T_2(sin(theta/2.0),0),T_3(-sin(theta/2.0),0),T_4(cos(theta/2.0),0);  
	  T<<T_1*I,T_2*I,T_3*I,T_4*I;
	  
	  
	  g_r_matrix<<g_r_up*I,Z,Z,g_r_down*I;
	  g_r_theta=T.inverse()*g_r_matrix*T;
	  
	  
	  t_matrix<<paras_t_up[cleavage]*I,Z,Z,paras_t_up[cleavage]*I;        
	  
	  
	  g_l_matrix<<g_l_up*I,Z,Z,g_l_down*I;
	  
	  
	  
	  A=(IL-(g_r_theta*t_matrix.adjoint()*g_l_matrix*t_matrix)).inverse();
	  B=(IL-(g_r_theta.adjoint()*t_matrix.adjoint()*g_l_matrix.adjoint()*t_matrix)).inverse();
	  
	  double integrand=real(f1-f2)*(1/(8*pow(M_PI,3)))*real(((g_l_matrix*t_matrix*A*B*g_r_theta.adjoint()*t_matrix.adjoint()-A*B+0.5*(A+B))*Pauli).trace());
	  
	  
	  n_tot_theta+=1;
	  
	  if(n_tot_theta<=1.1 && n_tot_theta>=0.9)
	    {
	      s_theta+=integrand;
	    }
	  else if(n_tot_theta<=11.1 && n_tot_theta>=10.9)
	    {
	      s_theta+=integrand;
	    }
	  else
	    {
	      s_theta+=2*integrand;
	    } 
	  
	  
	}                       /*end of theta loop*/
    
    s_theta=s_theta*0.5*(M_PI/(n_tot_theta-1));          /*trapezium rule for theta*/
    
    
    
  return s_theta;

      
    }
}






double  eqm_func()
{
  big_flag=0;                            /*flag to tell f to do eqm integrand*/
  
  for(int eqm_flag=0;eqm_flag<=1;eqm_flag+=1)     /*does left and right eqm parts seperately*/
    {
      
  double dq;     
  double conv[num_spacer_eqm];
  double s[num_spacer_eqm];
  double b[2][num_spacer_eqm];

  for(int i=0;i<num_spacer_eqm;i++)
    {
      s[i]=0.0;
    }
   


  int matsub=15;                                                     /*Number of terms in Matsubara sum- VARIABLE*/
  complex<double> E_s[matsub];

  
  if(eqm_flag==0)                                         /*left eqm matsub frequencies*/
    {
      for(int m=0;m<matsub;m++)
	{
	  complex<double> dummy(voltage,KbT*M_PI*((2.0*m)+1.0)+delta);              /*Fermi level on left*/
	  E_s[m]=dummy;                           
	}
    }
  else                                                 /*right eqm matsub frequencies*/
    {
      for(int m=0;m<matsub;m++)
	{
	  complex<double> dummy(0,KbT*M_PI*((2.0*m)+1.0)+delta);              /*Fermi level on left*/
	  E_s[m]=dummy;                                               
	}
    }

  int grid_num=20;                                      /*number of different grids- VARIABLE*/

  int c[grid_num];
  c[0]=2;
  for(int i=1;i<grid_num;i++)
    {
      c[i]=2*c[i-1];
    }

                                         

   for(int m=0;m<matsub;m++)                            /*mat loop start*/
	 {

	   for(int i=0;i<num_spacer_eqm;i++)
	     {
	       b[0][i]=10000;               	     
	       b[1][i]=0.0;
	     }
	   	   
	   for(int k=0;k<grid_num;k++)                                     /*start of convergence loop*/
	     {
	       
	       int n_tot=0;                                           /*total number of points in grid*/

	       for(int i=0;i<num_spacer_eqm;i++)
		 {
		   b[1][i]=0.0;
		 }
	       
	       dq=M_PI/(c[k]-1);                                /*grid point gap size*/
	       
	       double a[c[k]];                                 /*holds coordinate points*/
	       a[0]=0;
	       for(int i=1;i<c[k];i++)
		 {
		   a[i]=a[i-1]+dq;
		 }
	       
	       	       
	                                             /*lower left corner*/
		   int weight1=1;		     
		   n_tot=n_tot+weight1;
		   x=a[0];
		   y=a[0];
		   f(E_s[m]);
		   for(int i=0;i<num_spacer_eqm;i++)
		     {
		       b[1][i]=b[1][i]+(weight1* arr[i]);
		     }
	       
	     
	                                     /*upper right corner*/
		 
		   int weight2=1;
		   n_tot=n_tot+weight2;
		   x=a[c[k]-1];
		   y=a[c[k]-1];
		   f(E_s[m]);
		    for(int i=0;i<num_spacer_eqm;i++)
		      {
			b[1][i]=b[1][i]+(weight2* arr[i]);
		      }
	       
	     
	                                  /*lower right corner*/
		 
		   int weight3=2;
		   n_tot=n_tot+weight3;
		   x=a[c[k]-1];
		   y=a[0];
		   f(E_s[m]);
		    for(int i=0;i<num_spacer_eqm;i++)
		      {
			b[1][i]=b[1][i]+(weight3* arr[i]);
		      }
       	       
		 
		   for(int i=1;i<=c[k]-2;i++)                           /*lower edge*/
		     {
		       int weight4=4;
		       n_tot=n_tot+weight4;
		       x=a[i];
		       y=a[0];
		       f(E_s[m]);
		        for(int j=0;j<num_spacer_eqm;j++)
			  {
			    b[1][j]=b[1][j]+(weight4*arr[j]);
			  }
		     }
		 
	       


	       
	       
		 
		   for(int i=1;i<=c[k]-2;i++)                           /*right edge*/
		     {
		       int weight5=4;
		       n_tot=n_tot+weight5;
		       x=a[c[k]-1];
		       y=a[i];
		       f(E_s[m]);
		        for(int j=0;j<num_spacer_eqm;j++)
			  {
			    b[1][j]=b[1][j]+(weight5*arr[j]);
			  }
		       
		     }
		       
		 
		   for(int i=1;i<=c[k]-2;i++)                        /*diagonal edge*/
		     {
		       int weight6=4;
		       n_tot=n_tot+weight6;
		       x=a[i];
		       y=a[i];
		       f(E_s[m]);
		        for(int j=0;j<num_spacer_eqm;j++)
			  {
			    b[1][j]=b[1][j]+(weight6* arr[j]);
			  }
		       
		     }
	 	       
		 
		   for(int i=2;i<=c[k]-2;i++)                              /*inside*/
		     {
		       for(int j=1;j<=i-1;j++)
			 {
			   int weight7=8;
			   n_tot=n_tot+weight7;
			   x=a[i];
			   y=a[j];
			   f(E_s[m]);
			    for(int k=0;k<num_spacer_eqm;k++)
			      {
				b[1][k]=b[1][k]+(weight7* arr[k]);
			      }
			   
			 }
		     }
		 	       

	       for(int i=0;i<num_spacer_eqm;i++)
		 {
		   b[1][i]=((b[1][i])*pow(2*M_PI,2))/n_tot;
		   conv[i]=abs(b[1][i]-b[0][i]);	     	     	     	     
		   b[0][i]=b[1][i];
		 }

	       sort(conv,conv+(sizeof(conv)/sizeof(conv[0])));
	     	     
	     
             if(conv[num_spacer_eqm-1]<0.00000001)              /*precision of integral convergence- VARIABLE*/
	       {
		 for(int i=0;i<num_spacer_eqm;i++)
		   {
		     s[i]=s[i]+b[1][i];
		   }
	         
	         break;
               }
	     }                                                   /*end of convergence loop*/
	 }                                             /*end of matsub loop*/
   
   
      
   
   if(eqm_flag==0)              /*stores final eqmIEC values in global arrays*/
     {
       left_eqm[0]=s[0];
       left_eqm[1]=s[1];
       left_eqm[2]=s[2];
       left_eqm[3]=s[3];
       left_eqm[4]=s[4];
     }
   else
     {
      right_eqm[0]=s[0];
      right_eqm[1]=s[1];
      right_eqm[2]=s[2];
      right_eqm[3]=s[3];
      right_eqm[4]=s[4]; 
     }
   
    }                            /*end of eqm_flag loop*/

   return 0;
   
}





double sub_func(double variable)               /*adds delta to energy here so f can have complex number as input*/
{
  complex<double> complex_energy(variable,delta);

  return f(complex_energy);
}




static void NAG_CALL f(const double x[], Integer nx, double fv[],
                       Integer *iflag, Nag_Comm *comm) {
  Integer i;

  /* Set iflag negative to terminate execution for any reason. */
  *iflag = 0;

  for (i = 0; i < nx; i++) {
    if (x[i] == 1.0) {
      *iflag = -1;
      /* Store chosen value of iflag in iuser */
      comm->iuser[0] = *iflag;
    }
  }
  if (*iflag == 0) {
    for (i = 0; i < nx; i++)
      {
	fv[i] = sub_func(x[i]);
      }
  }
}





double ooe_func()
{
  big_flag=1;                          /*tells f to do ooe integrand*/

  
   double grid_sum=0.0;                          /*finds total number of points in k-space grid*/
  for(int i=1;i<=grid_size;i++)
    {
      grid_sum+=i;
    }


 
  
  Integer exit_status = 0;

   double total_excluded=0.0;                                 /*number of k points thrown out*/

  
   for(num_spacer=0;num_spacer<=4;num_spacer+=1)              /*runs through each spacer size seperately*/
    {
           

	    double n_tot=0.0;
	    double s=0.0;
	  
	  for(x=0;x<=M_PI+0.0001;x=x+M_PI/(grid_size-1))
	    {
	      for(y=0;y<=x+0.0001;y=y+M_PI/(grid_size-1))
		{

	        

		  int flag_nathan=0;                                     
		  
		  /*start of adaptive integration*/
		  double a, b, result, epsabs, epsrel, abserr;
		  Integer lrinfo, liinfo, maxsub;
		  Integer *iinfo = 0;
		  double *rinfo = 0;
		  
		  /* Nag Types */
		  NagError fail;
		  Nag_Comm comm;
		  
		  INIT_FAIL(fail);
		  
		  epsabs = 0.000001;                                                              /*absolute error- VARIABLE*/
		  epsrel = 0.000001;                                                              /*relative error- VARIABLE*/
		  
		  /*NOTE- Limits now fermi levels not direct energies*/
		  
		  double extra_bit=real(KbT)*log(1000000000);    /*Extra bit either side of fermi levels to be integrated over*/
		  
		  
		  if(voltage<0)                                      /*a lower limit, b upper limit*/
		    {
		      a=(right_fermi+voltage)-extra_bit;
		      b=right_fermi+extra_bit;
		    }
		  else
		    {
		      a=right_fermi-extra_bit;
		      b=(right_fermi+voltage)+extra_bit;
		    }
		  
		  
		  maxsub = 2000;                                         /*max number of intervals- VARIABLE*/
		  lrinfo = 4 * maxsub;
		  liinfo = MAX(maxsub, 4);
		  
		  /* Allocate memory */
		  if (!(rinfo = NAG_ALLOC(lrinfo, double)) ||
		      !(iinfo = NAG_ALLOC(liinfo, Integer))) {
		    printf("Allocation failure\n");
		    exit_status = -1;
		    goto END;
		  }
		  
		  
		  /* Evaluate the integral using the vectorized One-dimensional adaptive
		   * quadrature routine nag_quad_dim1_fin_general (d01rjc).
		   */
		  nag_quad_dim1_fin_general(f, a, b, epsabs, epsrel, maxsub, &result, &abserr,
					    rinfo, iinfo, &comm, &fail);
		  if (fail.code != NE_NOERROR)
		    {
		      flag_nathan=1;
		      /*cout<<"At x="<<x<<" ,y="<<y<<" :"<<endl;
			printf("Error or warning from nag_quad_dim1_fin_general (d01rjc) %s\n",fail.message);*/
		    }
		  
		  if (fail.code == NE_USER_STOP)
		    {
		      cout<<"At x="<<x<<" ,y="<<y<<" :"<<endl;
		      printf("Exit requested from f with iflag = %" NAG_IFMT "\n", comm.iuser[0]);
		    }
		  else
		    {
		      exit_status = 1;
		      goto END;
		    }
		  
		END:
		  
		  NAG_FREE(rinfo);
		  NAG_FREE(iinfo);                                                    /*end of adaptive integration*/

		  

		  
		  double incre=0.0000001;
		  
		  if(flag_nathan==1)                        
		    {
		      total_excluded+=1;
		      n_tot=n_tot;
		      s=s;
		    }
		  else if(x<=0+incre && x>=0-incre && y<=0+incre && y>=0-incre)                           /*k-space weights*/
		    { 
		      n_tot+=1;
		      s=s+result;
		    }
		  else if(x<=M_PI+incre && x>=M_PI-incre && y<=M_PI+incre && y>=M_PI-incre)
		    {
		      n_tot+=1;
		      s=s+result;
		    }
		  else if(x<=M_PI+incre && x>=M_PI-incre && y<=0+incre && y>=0-incre)
		    {
		      n_tot+=2;
		      s=s+(2*result);
		    }
		  else if(y<=0+incre && y>=0-incre && x>0 && x<M_PI)
		    {
		      n_tot+=4;
		      s=s+(4*result);
		    }
		  else if(x<=M_PI+incre && x>=M_PI-incre && y>0 && y<M_PI)
		    {
		      n_tot+=4;
		      s=s+(4*result);
		    }
		  else if(x<=y+incre && x>=y-incre && x>0 && x<M_PI && y>0 && y<M_PI)
		    {
		      n_tot+=4;
		      s=s+(4*result);
		    }
		  else
		    {
		      n_tot+=8;
		      s=s+(8*result);
		    }
		  
		  
		  
		}
	    }
	
      
	  s=s*((4*pow(M_PI,2))/n_tot);               /*final step of k|| integration*/

	  holder[num_spacer]=s;              /*stores in global array*/

	

    }                                                   /*end of num_spacer loop*/

   

      if(total_excluded>=grid_sum*11)                                     
	{
	  cout<<"ERROR!!!! Either major convergence or licence key issue"<<endl;
	  }
      
        
  
      return exit_status;
}




int main()
{

  auto start=high_resolution_clock::now();

  

  eqm_func();          /*calls to eqm and ooe funcs to find IEC components and store in global arrays*/
  ooe_func();

  double zero=left_eqm[0]+right_eqm[0]+holder[0];          /*adds all the components together*/
  double one=left_eqm[1]+right_eqm[1]+holder[1];
  double two=left_eqm[2]+right_eqm[2]+holder[2];
  double three=left_eqm[3]+right_eqm[3]+holder[3];
  double four=left_eqm[4]+right_eqm[4]+holder[4];
  

  
  cout<<"equilibrium terms"<<endl;
  cout<<"["<<voltage<<",["<<left_eqm[0]+right_eqm[0]<<","<<left_eqm[1]+right_eqm[1]<<","<<left_eqm[2]+right_eqm[2]<<","<<left_eqm[3]+right_eqm[3]<<","<<left_eqm[4]+right_eqm[4]<<"]],"<<endl<<endl;


   cout<<"out-of-equilibrium terms"<<endl;
  cout<<"["<<voltage<<",["<<holder[0]<<","<<holder[1]<<","<<holder[2]<<","<<holder[3]<<","<<holder[4]<<"]],"<<endl<<endl;
  

  cout<<"Total ooeIEC"<<endl;
  cout<<"["<<voltage<<",["<<zero<<","<<one<<","<<two<<","<<three<<","<<four<<"]],"<<endl<<endl;

 
  

    
   auto stop=high_resolution_clock::now();
   auto duration=duration_cast<seconds>(stop-start);
   cout<<"time= "<<duration.count()<<" seconds"<<endl<<endl;

   return 0;
}
    
