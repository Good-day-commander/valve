#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cmath>
	/*Care for streaming in tables and update x in the extrapolation */
#include <string.h>
#include <string>
#include <omp.h>
#include <cstdlib> 
#include <ctime>
#include <vector>
#include <time.h>
#include <complex>


// 0.4um당 1칸
// 8um짜리

//#define     M_PI        3.14159265358979


#define     s_stif      1.0//(0.01)
#define     b_stif     -0.2//(-0.2)

#define		r_stif		0.00007//0.00007

#define     v_stif     0//(0.00)
#define     al_stif     0.02//(0.02) 
#define     at_stif     0.001//(0.05)


//#define     v_stif     (0.5)
//#define     al_stif     (0.02) 
//#define     at_stif     (0.05)

#define     De     (-0.005)
#define     s_beta     (0.1) 
#define     bumwe     (1.5)
#define     safezone    (20.0)


#define		fakex     	    1

#define		Q2          27 //JungShin
#define         epsilon2     0.000001 // JungShin
#define         epsilon     0.000001 // JungShin
#define		slip_ratio    0.0

//#define		boxXLENGTH     	86//(20+1.5)*4 + 1.5*2
//#define		boxYLENGTH     	86 //22 //20 + (1.2)*2 -1
//#define		boxZLENGTH     	86// (6.6+1.f7)*10 + 2*2



#define		chLENGTH		1000//800 //채널의 길이
#define		chRADIUS		35//45//30	//채널의 직경


#define		XLENGTH			180//1000//boxXLENGTH //*2+chLENGTH
#define		YLENGTH			101//boxYLENGTH     
#define		ZLENGTH			101//boxZLENGTH


#define		ALENGTH			100
#define		CRADIUS			40


#define		peri_ref        (XLENGTH/2.0)
#define		tau           1.0

#define		pipe_r         chRADIUS //43

#define		uMax           (0.03)    //(0.00157) 줬을때 0.0006이 평균속도로 나옴
//#define     forcepowerAC    (uMax*(tau-0.5)/3.14/(pipe_r)/(pipe_r)*4.0*1.05)
//#define     forcepowerAC    (uMax*(tau-0.5)/3.14/(pipe_r)/(pipe_r)*4.0*1.00*2.0)
//#define     forcepowerAC  0.0

#define		rbcn   3//적혈구 갯수
//#define		p_num       500
//#define     f_main      996
//#define     l_main       1494

#define		p_num       3200//20000
#define     f_main      6396//39996
#define     l_main      9437//59597 //선의 수 = 6rc - 2r - 2c -3

#define		tp_num       (rbcn*p_num)
#define     tf_main       (rbcn*f_main)
#define     tl_main       (rbcn*l_main)



//#define     Radius  (8.5*2.0/7.7924025)
#define     Radius  1.0//(10.0)
#define     off_set_x  ((XLENGTH+1)/2.0)
#define     off_set_y  ((YLENGTH+1)/2.0)
#define     off_set_z  ((ZLENGTH+1)/2.0)
#define     offzz      7.5
#define     offxx      2.5
#define     offyy      10.0



#define     forcepower   0.0  //pN
#define     opforce2     (2.77/100000*forcepower/10)


#define     wallef       (0.1) 
#define     walldis       (1.1) 


#define		pulse_step	20
#define		pulse_time	1000//2500
#define		pulse_initial 20000

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





#define		Q           	19
#define		w0              1.0 / 3.0
#define		w1              1.0 / 18.0
#define		w2              1.0 / 36.0
#define		nu          XLENGTH*YLENGTH*ZLENGTH
#define		RHO_0	        1.0

#define		gamma           1.0


#define		vis_c			0.8 
#define		invis			2.0 
#define		outvis			1.0 

#define		T				100000//100000
#define		ptime           500//500//1000 // mesh
#define		pptime          500//500 //1000//  flow rate
#define		ppptime         500//500//1000//   data

#define		aftermoving        1//
#define		afterprint          1//


#define		reintime        0


#define		retime			1000000


#define		divstep			100
#define		pulsefrequency	20

/*
#define		seng	        0.3568
#define		mu0             0.056
#define		muinf           0.0035
#define		ramda           3.313
*/
#define	    i_main2         (XLENGTH*YLENGTH*ZLENGTH*Q+1)
#define     i_main3         (XLENGTH*YLENGTH*ZLENGTH+1)
#define     i_main4         (XLENGTH*YLENGTH*ZLENGTH*3*3+1)


#define		PI2				3.14159265

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void        tables                              ( void );
void        propagate                           ( void );
void        initialise                          ( void );
void        calc_obs                            ( void );
void        project                             ( void );
void        density_data                        ( void );
void		pressure_data						( void );
void		timestep_data						( void );
void        show_lattice                        ( void );
void		impedance							( void );
void		initial_density_data				( void );
void	    pulsatile                           ( void );
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void		setparticle							( void );
void		ParticleContour						( void );
void		TotalForce   						( void );
void		Rigidspring 						( void );
void		UpdatePosition 						( void );
void        initialise_2(void);
void	flowrate(void);
void    viscosity(void);
void    viscosityshear(void);
void   intial_viscosity(void);

void		normalvc(void);

void		oscillation							(void);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void		mesh           						( void );
void		mesh_nor     						( void );
void		mesh_com     						( void );
void		mesh_com2							(void);
void	    springcal                           ( void );
void		bending        						( void );
void		area_vol_con  						( void );
void		constructimage                      (void);
void		aggregation							(void);
void		recovery							(void);
void		wallefcal						    (void);
void		wallefcal_shear                     (void);

void		node_track                     (void);



void		chh(void);
void restartin(void);
void initial_vol_area(void);
void    restart_write1(void);
void    restart_write2(void);
void    restart_write3(void);
void    restart_write4(void);
void    restart_write5(void);




void	    pulsatile                           ( void );
void	    pulse_flow                          ( void );
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double		*b = new double[i_main2];
double		*back = new double[i_main2];
double		*beq = new double[i_main2];


int	       next_x         [XLENGTH+1][Q];
int	       next_y	      [YLENGTH+1][Q];
int	       next_z	      [ZLENGTH+1][Q];
int         cx                          [Q] = { 0,-1,-1,-1,-1,-1, 0, 0, 0, 0,1,1,1, 1, 1,0,0, 0,0};
int         cy                          [Q] = { 0,-1, 0, 0, 0, 1,-1,-1,-1, 0,1,0,0, 0,-1,1,1, 1,0};
int         cz                          [Q] = { 0, 0,-1, 0, 1, 0,-1, 0, 1,-1,0,1,0,-1, 0,1,0,-1,1};

int	         v	        [Q][3];
int	        *bnode = new int[i_main3];


double		t                           [Q] = {w0,w2,w2,w1,w2,w2,w2,w1,w2,w1,w2,w2,w1,w2,w2,w2,w1,w2,w1};

double      *rhoN = new double[i_main3];    
double	    *ux	= new double[i_main3];         
double	    *uy = new double[i_main3];          
double	    *uz = new double[i_main3];          
double	    *rho = new double[i_main3];        

double	    *inout_vis = new double[i_main3];
int	        *vis_check = new int[i_main3];



double	    *press = new double[i_main3];

double		stiff_s[rbcn+1];
double		stiff_b[rbcn+1];
double		stiff_at[rbcn+1];
double		stiff_al[rbcn+1];
double		stiff_v[rbcn+1];


double		cen_r[rbcn+1];
double		cen_rx[rbcn+1];
double		cen_ry[rbcn+1];
double		cen_rz[rbcn+1];
double		an_r[rbcn + 1];

//double      *tau = new double[i_main3];  
//double      *mu = new double[i_main3];  
//double      *summ = new double[i_main3];  

//double      *strate = new double[i_main4];

////////////////////////////////////////////////////////////////////////////////
//IBM
/*
int	        *image_particle = new int[i_main3];
double	    *FforceX = new double[i_main3];          
double	    *FforceY = new double[i_main3];          
double	    *FforceZ = new double[i_main3];        


double      um           [Xm_number+1][3];
double      S_Force      [Xm_number+1][3];
double      Xm_location  [Xm_number+1][3];
double      W_f          [Xm_number+1];
double      RXm_location [Xm_number+1][3];
double      pf           [Xm_number+1][3];
*/
//////////////////////////////////////////////////////////////////////
//Real IBM
int	        *image_particle = new int[i_main3];
int	        *rpoint = new int[i_main3];

double      um           [tp_num+1][3];
double      umD           [tp_num+1][3];

double	    *FforceX = new double[i_main3];          
double	    *FforceY = new double[i_main3];          
double	    *FforceZ = new double[i_main3];  

double		*point_x = new double[tp_num+1];
double		*point_y = new double[tp_num+1];
double		*point_z = new double[tp_num+1];

double		*point_xi = new double[tp_num+1];
double		*point_yi = new double[tp_num+1];
double		*point_zi = new double[tp_num+1];

double	    *tforceX = new double[tp_num+1];          
double	    *tforceY = new double[tp_num+1];          
double	    *tforceZ = new double[tp_num+1];  

double	    *sforceX = new double[tp_num+1];          
double	    *sforceY = new double[tp_num+1];          
double	    *sforceZ = new double[tp_num+1];  

double	    *rforceX = new double[tp_num+1];          
double	    *rforceY = new double[tp_num+1];          
double	    *rforceZ = new double[tp_num+1];  


double	    *bforceX = new double[tp_num+1];          
double	    *bforceY = new double[tp_num+1];          
double	    *bforceZ = new double[tp_num+1];  
	
double	    *vforceX = new double[tp_num+1];          
double	    *vforceY = new double[tp_num+1];          
double	    *vforceZ = new double[tp_num+1];  

double	    *aforceX = new double[tp_num+1];          
double	    *aforceY = new double[tp_num+1];          
double	    *aforceZ = new double[tp_num+1];  

double	    *agforceX = new double[tp_num + 1];
double	    *agforceY = new double[tp_num + 1];
double	    *agforceZ = new double[tp_num + 1];


double	    *aforceX2 = new double[tp_num + 1];
double	    *aforceY2 = new double[tp_num + 1];
double	    *aforceZ2 = new double[tp_num + 1];

double	    *aforceX3 = new double[tp_num + 1];
double	    *aforceY3 = new double[tp_num + 1];
double	    *aforceZ3 = new double[tp_num + 1];



double	    *Deff = new double[tp_num + 1];
double	    *defn = new double[tp_num + 1];
//double		*nor_vec = new double[fmain2+1]; ;



//////////////////////////////////////////////////////////
int    	    *vtx0 =			new int[tf_main+1];
int    	    *vtx1 =			new int[tf_main+1];
int    	    *vtx2 =			new int[tf_main+1];
int    	    *vtx3 = new int[tf_main + 1];


         
double		*n_vect_x =		new double[tf_main+1];
double		*n_vect_y =		new double[tf_main+1];
double		*n_vect_z =		new double[tf_main+1];

double		*area_c =		new double[tf_main+1];
double		*area_r =		new double[tf_main+1];


double		*vol_c =		new double[tf_main+1];
double		*vol_r =		new double[tf_main+1];



double		*p_angle =		new double[tl_main+1];

int		*line_twopoint1 = new int[tl_main+1];
int		*line_twopoint2 = new int[tl_main+1];

int		*line_face1 = new int[tl_main+1];
int		*line_face2 = new int[tl_main+1];

int		*line_verpoint1 = new int[tl_main+1];
int		*line_verpoint2 = new int[tl_main+1];

double		*r_angle =		new double[tl_main+1];
double		*line_ref = new double[tl_main+1];

double      CoMz[rbcn + 1];



///////////////////////////////////////////////////////////////




int	        *slip_node = new int[i_main3];
int	        *SN_x_2 = new int[i_main3];
int	        *SN_y_2 = new int[i_main3];
int	        *SN_z_2 = new int[i_main3];



int	       next_x2[XLENGTH + 1][Q2];
int	       next_y2[YLENGTH + 1][Q2];
int	       next_z2[ZLENGTH + 1][Q2];


double		*back_2 = new double[i_main2];



int         cx2[Q2] = { 0, 1, 0, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, -1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1 };
int         cy2[Q2] = { 0, 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1 };
int         cz2[Q2] = { 0, 0, 0, 1, 0, 0, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0 };


int opposite[Q] = { 0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9 };



double      pulse           [20];
double		u_pulse;

double forcepowerAC;

////////////////////////////////////////////////////////////////////////////////////
double	    OMEGA_B, opforce, vw_f, O_xl, O_yl, O_zl, flowori;
double	  r_tot_vol[rbcn + 1], r_tot_area[rbcn + 1];
double	  cen_x[rbcn + 1], cen_y[rbcn + 1], cen_z[rbcn + 1], time_cx, meanvel_t;
int         Timek, vvol;

double		uOsci;// oscillating wall speed

int p_phase;

using namespace std;




   
int	main	(void)
{
    //int count;
	int jj;
	p_phase = 0;
	vw_f = Radius * Radius * 4.0 * 3.14 / double(p_num);
		
	opforce = 0.0;
	Timek = 0;
	time_cx = 0.0;
	forcepowerAC = 0;
    tables		();



    OMEGA_B = 1.0;

//	pulse_step = 0;
	pulsatile();
	pulse_flow();
	oscillation();
    initialise	();
	constructimage();
	initialise_2();

	mesh ();
	//intial_viscosity ();
//	restartin();
	
	chh();
	normalvc();
	initial_vol_area();


	for (Timek = reintime; Timek < T; Timek++)
	{
		//if (time_cx > 0.99) break;
		pulse_flow();
		//oscillation();
		if (Timek > aftermoving)
		{
			normalvc();
			springcal();
			recovery();
			bending();
			aggregation();
			area_vol_con();


			timestep_data();
			wallefcal();
			//wallefcal_shear();
			TotalForce();
		}

		propagate();
		calc_obs();

		project();

		if (Timek == 1)
		{

			//ParticleContour();
			density_data();
			//			mesh_nor ();
			mesh_com();
			
			//chh();

		}

		if ( Timek > afterprint)
		{
			if (Timek >= 1) if (Timek % ppptime == 0)
		{

			//ParticleContour();
			density_data();
			//			mesh_nor ();

			
			//chh();

		}

		if (Timek >= 1) if (Timek % ptime == 0)
		{

			//ParticleContour();

			//			mesh_nor ();
			mesh_com();
			
			//chh();
		}
		}
		if (Timek >= 1) if (Timek % retime == 0)
		{


			restart_write1();
			restart_write2();
			restart_write3();
			restart_write4();
			restart_write5();

		}


		if (Timek > aftermoving)
		{
			UpdatePosition();
			node_track();
		}
		//viscosity();

		if (Timek >= 1) if (Timek % pptime == 0)
		{
			flowrate();
			viscosityshear();

		}
		
		if (Timek == afterprint)
		{
			flowrate();
		}
    }


    scanf("%d",&Timek);


	delete[]b;
	delete[]back;
	delete[]press;
	delete[]rhoN;

	delete[]rho;
	delete[]bnode;

	delete[]ux;
	delete[]uy;
	delete[]uz;
	//delete[]tau;
	//delete[]mu;
    delete[]beq;
//	delete[]summ;
	//delete[]strate;
    delete[]image_particle;
	
	//delete[]S_Force;
	delete[]FforceX;
	delete[]FforceY;
	delete[]FforceZ;
	delete[]um;
	//delete[]Xm_location;
	//delete[]W_f;
    //delete[]RXm_location;
	//delete[]pf;

	delete[]p_angle;
	
	delete[]point_x;
	delete[]point_y;
	delete[]point_z;

	delete[]Deff;
	delete[]defn;

	delete[]inout_vis;
	delete[]vis_check;
	

	delete[]vtx0;
	delete[]vtx1;
	delete[]vtx2;
	delete[]vtx3;
	
	delete[]n_vect_x;
	delete[]n_vect_y;
	delete[]n_vect_z;
	

	delete[]area_c;
	delete[]area_r;
	delete[]vol_c;
	delete[]vol_r;


	delete[]line_ref;

	delete[]r_angle;


	delete[]rpoint;

	delete[]CoMz;









	delete[]sforceX;
	delete[]sforceY;
	delete[]sforceZ;
	
	delete[]rforceX;
	delete[]rforceY;
	delete[]rforceZ;

	delete[]aforceX;
	delete[]aforceY;
	delete[]aforceZ;


	delete[]agforceX;
	delete[]agforceY;
	delete[]agforceZ;


	delete[]aforceX2;
	delete[]aforceY2;
	delete[]aforceZ2;

	delete[]aforceX3;
	delete[]aforceY3;
	delete[]aforceZ3;







	return  0;
}


void chh()
{
	int i_nd;

	char FileNameAngle[100];
	char string[100] = { '0' };
	sprintf(FileNameAngle, "angle_%d.txt", Timek);
	FILE *fP4;


	fP4 = fopen(FileNameAngle, "w");








	for (i_nd = 1; i_nd <= l_main; i_nd++)
	{


		fprintf(fP4, "%lf   %lf\n", *(r_angle + i_nd), *(p_angle + i_nd) );


	}
	fprintf(fP4, "\n");
	fprintf(fP4, "\n");





	fflush(fP4);
	fclose(fP4);

	return;

}

void initial_vol_area(void)
{

	int     i, i_nd, i_nd2, pn1, pn2, pn3;
	double  ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3;

	FILE *ff3;



	for (i = 1; i <= rbcn; i++)
	{
		cen_x[i] = cen_y[i] = cen_z[i] = 0.0;
	}


	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	{
		for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
		{
			cen_x[i_nd2] += *(point_x + i_nd);
			cen_y[i_nd2] += *(point_y + i_nd);
			cen_z[i_nd2] += *(point_z + i_nd);


		}
	}


	for (i = 1; i <= rbcn; i++)
	{
		cen_x[i] /= double(p_num);
		cen_y[i] /= double(p_num);
		cen_z[i] /= double(p_num);


	}


	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
	{
		pn1 = *(vtx0 + i_nd);
		pn2 = *(vtx1 + i_nd);
		pn3 = *(vtx2 + i_nd);

		ax1 = *(point_x + pn2) - *(point_x + pn1);
		ay1 = *(point_y + pn2) - *(point_y + pn1);
		az1 = *(point_z + pn2) - *(point_z + pn1);


		ax2 = *(point_x + pn3) - *(point_x + pn1);
		ay2 = *(point_y + pn3) - *(point_y + pn1);
		az2 = *(point_z + pn3) - *(point_z + pn1);

		ax3 = (*(point_x + pn3) + *(point_x + pn2) + *(point_x + pn1)) / 3.0;
		ay3 = (*(point_y + pn3) + *(point_y + pn2) + *(point_y + pn1)) / 3.0;
		az3 = (*(point_z + pn3) + *(point_z + pn2) + *(point_z + pn1)) / 3.0;




		*(area_r + i_nd) = (pow(((ay1*az2 - ay2*az1)*(ay1*az2 - ay2*az1) + (ax2*az1 - ax1*az2)*(ax2*az1 - ax1*az2) + (ax1*ay2 - ay1*ax2)*(ax1*ay2 - ay1*ax2)), 0.5)) / 2.0;

		//		ax3 = cen_x[i_nd2] - *(point_x + pn1);
		//		ay3 = cen_y[i_nd2] - *(point_y + pn1);
		//		az3 = cen_z[i_nd2] - *(point_z + pn1);
		//	*(vol_r + i_nd) = (pow(((ay1*az2 - ay2*az1) * ax3 + (ax2*az1 - ax1*az2) * ay3 + (ax1*ay2 - ay1*ax2) * az3)  *  ((ay1*az2 - ay2*az1) * ax3 + (ax2*az1 - ax1*az2) * ay3 + (ax1*ay2 - ay1*ax2) * az3), 0.5)) / 6.0;
		//Li
		*(vol_r + i_nd) = (ax3 * *(n_vect_x + i_nd) + ay3 * *(n_vect_y + i_nd) + az3 * *(n_vect_z + i_nd))* *(area_r + i_nd) / 3.0;


	}
	/////////////////////////////////                       배열로 만들어서 RBC 각각 선언해야함

	for (i = 1; i <= rbcn; i++)
	{
		r_tot_area[i] = 0.0;
		r_tot_vol[i] = 0.0;
	}
	
	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
	{

		r_tot_area[i_nd2] += *(area_r + i_nd);
		r_tot_vol[i_nd2] += *(vol_r + i_nd);

	}







	
	{


		ff3 = fopen("process.txt", "a");
		fprintf(ff3, "hello%lf\n%lf\n", r_tot_area, r_tot_vol);
		fflush(ff3);
		fclose(ff3);



	}







	return;
}

void area_vol_con(void)
{

	int     i_nd, i_nd2, pn1, pn2, pn3, i;
	double   ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, f_mag, aar;
	double   wx1, wx2, wx3, wy1, wy2, wy3, wz1, wz2, wz3, cen_fx, cen_fy, cen_fz, vec_mag;
	double	 c_tot_area[rbcn + 1], c_tot_vol[rbcn + 1], total_a_dif[rbcn + 1];


	FILE *ff3;

	aar = double(f_main);
	aar = sqrt(aar);

#pragma omp parallel
	{

#pragma omp	for
		for (i = 1; i <= rbcn; i++)
		{
			cen_x[i] = 0.0;
			cen_y[i] = 0.0;
			cen_z[i] = 0.0;
		}

#pragma omp for private(i_nd)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		{
			for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
			{
				cen_x[i_nd2] += *(point_x + i_nd);
				cen_y[i_nd2] += *(point_y + i_nd);
				cen_z[i_nd2] += *(point_z + i_nd);

				*(aforceX + i_nd) = 0.0;
				*(aforceY + i_nd) = 0.0;
				*(aforceZ + i_nd) = 0.0;

				*(aforceX2 + i_nd) = 0.0;
				*(aforceY2 + i_nd) = 0.0;
				*(aforceZ2 + i_nd) = 0.0;

				*(aforceX3 + i_nd) = 0.0;
				*(aforceY3 + i_nd) = 0.0;
				*(aforceZ3 + i_nd) = 0.0;

				*(vforceX + i_nd) = 0.0;
				*(vforceY + i_nd) = 0.0;
				*(vforceZ + i_nd) = 0.0;

			}

			cen_x[i_nd2] /= double(p_num);
			cen_y[i_nd2] /= double(p_num);
			cen_z[i_nd2] /= double(p_num);

		}




#pragma omp for private(i_nd, pn1, pn2, pn3, ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
		{


			pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);

			ax1 = *(point_x + pn2) - *(point_x + pn1);
			ay1 = *(point_y + pn2) - *(point_y + pn1);
			az1 = *(point_z + pn2) - *(point_z + pn1);


			ax2 = *(point_x + pn3) - *(point_x + pn1);
			ay2 = *(point_y + pn3) - *(point_y + pn1);
			az2 = *(point_z + pn3) - *(point_z + pn1);

			ax3 = (*(point_x + pn3) + *(point_x + pn2) + *(point_x + pn1)) / 3.0;
			ay3 = (*(point_y + pn3) + *(point_y + pn2) + *(point_y + pn1)) / 3.0;
			az3 = (*(point_z + pn3) + *(point_z + pn2) + *(point_z + pn1)) / 3.0;




			*(area_c + i_nd) = (pow(((ay1*az2 - ay2*az1)*(ay1*az2 - ay2*az1) + (ax2*az1 - ax1*az2)*(ax2*az1 - ax1*az2) + (ax1*ay2 - ay1*ax2)*(ax1*ay2 - ay1*ax2)), 0.5)) / 2.0;
			*(vol_c + i_nd) = (ax3 * *(n_vect_x + i_nd) + ay3 * *(n_vect_y + i_nd) + az3 * *(n_vect_z + i_nd))* *(area_c + i_nd) / 3.0;








			/*


			pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);

			ax1 = *(point_x + pn2) - *(point_x + pn1);
			ay1 = *(point_y + pn2) - *(point_y + pn1);
			az1 = *(point_z + pn2) - *(point_z + pn1);


			ax2 = *(point_x + pn3) - *(point_x + pn1);
			ay2 = *(point_y + pn3) - *(point_y + pn1);
			az2 = *(point_z + pn3) - *(point_z + pn1);



			ax3 = cen_x[i_nd2] - *(point_x + pn1);
			ay3 = cen_y[i_nd2] - *(point_y + pn1);
			az3 = cen_z[i_nd2] - *(point_z + pn1);



			*(area_c + i_nd) = (pow(((ay1*az2 - ay2*az1)*(ay1*az2 - ay2*az1) + (ax2*az1 - ax1*az2)*(ax2*az1 - ax1*az2) + (ax1*ay2 - ay1*ax2)*(ax1*ay2 - ay1*ax2)), 0.5)) / 2.0;
			*(vol_c + i_nd) = (pow(((ay1*az2 - ay2*az1) * ax3 + (ax2*az1 - ax1*az2) * ay3 + (ax1*ay2 - ay1*ax2) * az3)  *  ((ay1*az2 - ay2*az1) * ax3 + (ax2*az1 - ax1*az2) * ay3 + (ax1*ay2 - ay1*ax2) * az3), 0.5)) / 6.0;
			*/
		}

#pragma omp	for
		for (i = 1; i <= rbcn; i++)
		{
			c_tot_area[i] = 0.0;
			c_tot_vol[i] = 0.0;
			total_a_dif[i] = 0.0;

		}


#pragma omp for private(i_nd)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
		{

			c_tot_area[i_nd2] += *(area_c + i_nd);
			c_tot_vol[i_nd2] += *(vol_c + i_nd);

		}

#pragma omp	for
		for (i = 1; i <= rbcn; i++)
		{

			//total_a_dif[i] = -at_stif * (c_tot_area[i] - r_tot_area[i]) * fabs(c_tot_area[i] - r_tot_area[i]) / r_tot_area[i] / 2.0;
			total_a_dif[i] = -stiff_at[i] * (c_tot_area[i] - r_tot_area[i]) / r_tot_area[i];

		}



#pragma omp for private(i_nd, pn1, pn2, pn3, f_mag)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
		{
			pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);

			//f_mag = - v_stif * ( *(vol_c + i_nd) -  *(vol_r + i_nd) ) / *(vol_r + i_nd) * *(area_c + i_nd)/3.0;
			//f_mag = -v_stif * vw_f * (c_tot_vol[i_nd2] - r_tot_vol[i_nd2]) * fabs(c_tot_vol[i_nd2] - r_tot_vol[i_nd2]) / r_tot_vol[i_nd2] / 2.0 / 3.0;

			f_mag = -1.0*(stiff_v[i_nd2] * (c_tot_vol[i_nd2] - r_tot_vol[i_nd2]) * *(area_c + i_nd) / r_tot_vol[i_nd2] / 3.0);
			
			*(vforceX + pn1) += f_mag * *(n_vect_x + i_nd);
			*(vforceY + pn1) += f_mag * *(n_vect_y + i_nd);
			*(vforceZ + pn1) += f_mag * *(n_vect_z + i_nd);

			*(vforceX + pn2) += f_mag * *(n_vect_x + i_nd);
			*(vforceY + pn2) += f_mag * *(n_vect_y + i_nd);
			*(vforceZ + pn2) += f_mag * *(n_vect_z + i_nd);

			*(vforceX + pn3) += f_mag * *(n_vect_x + i_nd);
			*(vforceY + pn3) += f_mag * *(n_vect_y + i_nd);
			*(vforceZ + pn3) += f_mag * *(n_vect_z + i_nd);

		}


#pragma omp for private(i_nd, pn1, pn2, pn3, cen_fx, cen_fy, cen_fz, wx1, wy1, wz1, wx2, wy2, wz2, wx3, wy3, wz3, vec_mag, f_mag)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
		{
			pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);

			cen_fx = (*(point_x + pn1) + *(point_x + pn2) + *(point_x + pn3)) / 3.0;
			cen_fy = (*(point_y + pn1) + *(point_y + pn2) + *(point_y + pn3)) / 3.0;
			cen_fz = (*(point_z + pn1) + *(point_z + pn2) + *(point_z + pn3)) / 3.0;


			wx1 = *(point_x + pn1) - cen_fx;
			wy1 = *(point_y + pn1) - cen_fy;
			wz1 = *(point_z + pn1) - cen_fz;

			wx2 = *(point_x + pn2) - cen_fx;
			wy2 = *(point_y + pn2) - cen_fy;
			wz2 = *(point_z + pn2) - cen_fz;

			wx3 = *(point_x + pn3) - cen_fx;
			wy3 = *(point_y + pn3) - cen_fy;
			wz3 = *(point_z + pn3) - cen_fz;

			vec_mag = pow((wx1 * wx1 + wy1* wy1 + wz1 * wz1), 0.5);
			wx1 /= vec_mag;
			wy1 /= vec_mag;
			wz1 /= vec_mag;
			vec_mag = pow((wx2 * wx2 + wy2* wy2 + wz2 * wz2), 0.5);
			wx2 /= vec_mag;
			wy2 /= vec_mag;
			wz2 /= vec_mag;
			vec_mag = pow((wx3 * wx3 + wy3* wy3 + wz3 * wz3), 0.5);
			wx3 /= vec_mag;
			wy3 /= vec_mag;
			wz3 /= vec_mag;

			f_mag = (-stiff_al[i_nd2] * (*(area_c + i_nd) - *(area_r + i_nd)) / pow(*(area_r + i_nd), 0.5));



			*(aforceX2 + pn1) += f_mag * wx1;
			*(aforceY2 + pn1) += f_mag * wy1;
			*(aforceZ2 + pn1) += f_mag * wz1;

			*(aforceX2 + pn2) += f_mag * wx2;
			*(aforceY2 + pn2) += f_mag * wy2;
			*(aforceZ2 + pn2) += f_mag * wz2;

			*(aforceX2 + pn3) += f_mag * wx3;
			*(aforceY2 + pn3) += f_mag * wy3;
			*(aforceZ2 + pn3) += f_mag * wz3;



			f_mag = total_a_dif[i_nd2];



			*(aforceX3 + pn1) += f_mag * wx1;
			*(aforceY3 + pn1) += f_mag * wy1;
			*(aforceZ3 + pn1) += f_mag * wz1;

			*(aforceX3 + pn2) += f_mag * wx2;
			*(aforceY3 + pn2) += f_mag * wy2;
			*(aforceZ3 + pn2) += f_mag * wz2;

			*(aforceX3 + pn3) += f_mag * wx3;
			*(aforceY3 + pn3) += f_mag * wy3;
			*(aforceZ3 + pn3) += f_mag * wz3;





		}

#pragma omp for private(i_nd, pn1, pn2, pn3, cen_fx, cen_fy, cen_fz, wx1, wy1, wz1, wx2, wy2, wz2, wx3, wy3, wz3, vec_mag, f_mag)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
		{
			pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);



			*(aforceX + pn1) = *(aforceX2 + pn1) + *(aforceX3 + pn1);
			*(aforceY + pn1) = *(aforceY2 + pn1) + *(aforceY3 + pn1);
			*(aforceZ + pn1) = *(aforceZ2 + pn1) + *(aforceZ3 + pn1);

			*(aforceX + pn2) = *(aforceX2 + pn2) + *(aforceX3 + pn2);
			*(aforceY + pn2) = *(aforceY2 + pn2) + *(aforceY3 + pn2);
			*(aforceZ + pn2) = *(aforceZ2 + pn2) + *(aforceZ3 + pn2);


			*(aforceX + pn3) = *(aforceX2 + pn3) + *(aforceX3 + pn3);
			*(aforceY + pn3) = *(aforceY2 + pn3) + *(aforceY3 + pn3);
			*(aforceZ + pn3) = *(aforceZ2 + pn3) + *(aforceZ3 + pn3);


		}



	}

	if (Timek >= 1) if (Timek % pptime == 0)
	{

		ff3 = fopen("process.txt", "a");
		fprintf(ff3, "time %d    %lf   %lf\n", Timek, c_tot_vol[1], c_tot_area[1]);
		fflush(ff3);
		fclose(ff3);
	}








	return;
}




void normalvc(void)
{

	int     i_nd, i_nd2, pn1, pn2, pn3, pn4, pn5, pn6;
	double   nor_ax, nor_ay, nor_az, nor_bx, nor_by, nor_bz, nor_mag, PI;
	double   ax1, ay1, az1, ax2, ay2, az2, abc1, abc2;
	PI = 4.0 * atan(1.0);

#pragma omp parallel
	{

#pragma omp	for
		for (i_nd = 1; i_nd <= tf_main; i_nd++)
		{
			*(n_vect_x + i_nd) = 0.0;
			*(n_vect_y + i_nd) = 0.0;
			*(n_vect_z + i_nd) = 0.0;
		}




#pragma omp for private(pn1, pn2, pn3, nor_ax, nor_ay, nor_az, nor_bx, nor_by, nor_bz, nor_mag)	
		for (i_nd = 1; i_nd <= tf_main; i_nd++)
		{

			pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);


			nor_ax = *(point_x + pn2) - *(point_x + pn1);
			nor_ay = *(point_y + pn2) - *(point_y + pn1);
			nor_az = *(point_z + pn2) - *(point_z + pn1);

			nor_bx = *(point_x + pn3) - *(point_x + pn1);
			nor_by = *(point_y + pn3) - *(point_y + pn1);
			nor_bz = *(point_z + pn3) - *(point_z + pn1);

			*(n_vect_x + i_nd) = nor_ay*nor_bz - nor_az*nor_by;
			*(n_vect_y + i_nd) = nor_az*nor_bx - nor_ax*nor_bz;
			*(n_vect_z + i_nd) = nor_ax*nor_by - nor_ay*nor_bx;

			nor_mag = pow((*(n_vect_x + i_nd)**(n_vect_x + i_nd) + *(n_vect_y + i_nd)**(n_vect_y + i_nd) + *(n_vect_z + i_nd)**(n_vect_z + i_nd)), 0.5);

			*(n_vect_x + i_nd) /= nor_mag;
			*(n_vect_y + i_nd) /= nor_mag;
			*(n_vect_z + i_nd) /= nor_mag;


		}
	}

	return;
}



void bending(void)
{

	int     i_nd, i_nd2, pn1, pn2, pn3, pn4, pn5, pn6, ap1, ap2, ap3, ap4, ap5, ap6;
	double   nor_ax, nor_ay, nor_az, nor_bx, nor_by, nor_bz, nor_mag, PI, mx, my, mz;
	double   ax1, ay1, az1, ax2, ay2, az2, abc1, abc2, abc3, abc4, abc5;

	PI = 4.0 * atan(1.0);



#pragma omp parallel
	{




#pragma omp	for
		for (i_nd = 1; i_nd <= tl_main; i_nd++)
		{
			*(p_angle + i_nd) = 0.0;
		}




#pragma omp for private(i_nd,pn1, pn2, pn3, pn4, ap1, ap2, ap3, ap4, ap5, ap6, mx, my, mz, abc1, ax1, ay1, az1, abc3, abc4)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
			for (i_nd = ((i_nd2 - 1)*l_main + 1); i_nd <= l_main*i_nd2; i_nd++)
			{
			pn1 = *(line_face1 + i_nd);
			pn2 = *(line_face2 + i_nd);

			pn3 = *(line_verpoint1 + i_nd);
			pn4 = *(line_verpoint2 + i_nd);

			ap1 = *(vtx0 + pn1);
			ap2 = *(vtx1 + pn1);
			ap3 = *(vtx2 + pn1);

			ap4 = *(vtx0 + pn2);
			ap5 = *(vtx1 + pn2);
			ap6 = *(vtx2 + pn2);



			mx = (*(point_x + ap4) + *(point_x + ap5) + *(point_x + ap6)) / 3.0 - (*(point_x + ap1) + *(point_x + ap2) + *(point_x + ap3)) / 3.0;
			my = (*(point_y + ap4) + *(point_y + ap5) + *(point_y + ap6)) / 3.0 - (*(point_y + ap1) + *(point_y + ap2) + *(point_y + ap3)) / 3.0;
			mz = (*(point_z + ap4) + *(point_z + ap5) + *(point_z + ap6)) / 3.0 - (*(point_z + ap1) + *(point_z + ap2) + *(point_z + ap3)) / 3.0;

			abc1 = (mx * *(n_vect_x + pn1) + my * *(n_vect_y + pn1) + mz * *(n_vect_z + pn1));





			if (abc1 <= 0.0)
			{
				ax1 = *(n_vect_x + pn1) * *(n_vect_x + pn2);
				ay1 = *(n_vect_y + pn1) * *(n_vect_y + pn2);
				az1 = *(n_vect_z + pn1) * *(n_vect_z + pn2);
				abc3 = ax1 + ay1 + az1;
				abc4 = PI - acos(abc3);

			}


			if (abc1 > 0.0)
			{
				ax1 = *(n_vect_x + pn1) * *(n_vect_x + pn2);
				ay1 = *(n_vect_y + pn1) * *(n_vect_y + pn2);
				az1 = *(n_vect_z + pn1) * *(n_vect_z + pn2);
				abc3 = ax1 + ay1 + az1;
				abc4 = PI + acos(abc3);

			}


			abc3 = abc4 - *(r_angle + i_nd);


			
			if (abc3 < 0)
			{
				*(p_angle + i_nd) = -stiff_b[i_nd2] * (1.0 - cos(abc3));
			}

			if (abc3 >= 0)
			{
				*(p_angle + i_nd) = stiff_b[i_nd2] * (1.0 - cos(abc3));
			}
			

			/*
			if (abc3 < 0)
			{
				*(p_angle + i_nd) = -b_stif * abc3*abc3/2.0;
					
			}

			if (abc3 >= 0)
			{
				*(p_angle + i_nd) = b_stif * abc3* abc3/2.0;
			}
			*/
			//yamaguchi

			//*(p_angle + i_nd) = b_stif *tan((abc4 - *(r_angle + i_nd))/2.0);



			}






#pragma omp	for
		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{
			*(bforceX + i_nd) = 0.0;
			*(bforceY + i_nd) = 0.0;
			*(bforceZ + i_nd) = 0.0;
		}
	}


		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
			for (i_nd = ((i_nd2 - 1)*l_main + 1); i_nd <= l_main*i_nd2; i_nd++)
		{


			pn1 = *(line_verpoint1 + i_nd);
			pn2 = *(line_face1 + i_nd);
			pn3 = *(line_verpoint2 + i_nd);
			pn4 = *(line_face2 + i_nd);
			pn5 = *(line_twopoint1 + i_nd);
			pn6 = *(line_twopoint2 + i_nd);


			ax1 = *(p_angle + i_nd) *  *(n_vect_x + pn2);
			ay1 = *(p_angle + i_nd) *  *(n_vect_y + pn2);
			az1 = *(p_angle + i_nd) *  *(n_vect_z + pn2);
			*(bforceX + pn1) += ax1;
			*(bforceY + pn1) += ay1;
			*(bforceZ + pn1) += az1;

			*(bforceX + pn5) -= (ax1 / 2.0);
			*(bforceY + pn5) -= (ay1 / 2.0);
			*(bforceZ + pn5) -= (az1 / 2.0);

			*(bforceX + pn6) -= (ax1 / 2.0);
			*(bforceY + pn6) -= (ay1 / 2.0);
			*(bforceZ + pn6) -= (az1 / 2.0);


			ax2 = *(p_angle + i_nd) *  *(n_vect_x + pn4);
			ay2 = *(p_angle + i_nd) *  *(n_vect_y + pn4);
			az2 = *(p_angle + i_nd) *  *(n_vect_z + pn4);
			*(bforceX + pn3) += ax2;
			*(bforceY + pn3) += ay2;
			*(bforceZ + pn3) += az2;


			*(bforceX + pn5) -= (ax2 / 2.0);
			*(bforceY + pn5) -= (ay2 / 2.0);
			*(bforceZ + pn5) -= (az2 / 2.0);

			*(bforceX + pn6) -= (ax2 / 2.0);
			*(bforceY + pn6) -= (ay2 / 2.0);
			*(bforceZ + pn6) -= (az2 / 2.0);


		}


	




	return;
}



void restartin(void)
{


	int     i_in, i_nd;

	FILE *res1;
	FILE *res2;
	FILE *res3;
	FILE *res4;
	FILE *res5;



	char res_data1[100];
	char res_data2[100];
	char res_data3[100];
	char res_data4[100];
	char res_data5[100];

	sprintf(res_data1, "restart_1_2.txt");
	if ((res1 = fopen(res_data1, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_nd = 1; i_nd <= tf_main; i_nd++)
		{
			fscanf(res1, "%d %d %d", &*(vtx0 + i_nd), &*(vtx1 + i_nd), &*(vtx2 + i_nd));
		}

		
		fflush(res1);
		fclose(res1);
	}


	sprintf(res_data2, "restart_2_2.txt");
	if ((res2 = fopen(res_data2, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{
			fscanf(res2, "%lf %lf %lf", &*(point_x + i_nd), &*(point_y + i_nd), &*(point_z + i_nd));
		}


		fflush(res2);
		fclose(res2);
	}



	sprintf(res_data3, "restart_3_2.txt");
	if ((res3 = fopen(res_data3, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_nd = 1; i_nd <= tl_main; i_nd++)
		{
			fscanf(res3, "%lf %lf", &*(line_ref + i_nd), &*(r_angle + i_nd));
		}


		fflush(res3);
		fclose(res3);
	}


	sprintf(res_data4, "restart_4_2.txt");
	if ((res4 = fopen(res_data4, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_nd = 1; i_nd <= tl_main; i_nd++)
		{
			fscanf(res4, "%d %d %d %d %d %d", &*(line_face1 + i_nd), &*(line_face2 + i_nd), &*(line_verpoint1 + i_nd), &*(line_verpoint2 + i_nd), &*(line_twopoint1 + i_nd), &*(line_twopoint2 + i_nd));
		}


		fflush(res4);
		fclose(res4);
	}




	sprintf(res_data5, "restart_5_2.txt");
	if ((res5 = fopen(res_data5, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_nd = 1; i_nd <= i_main2; i_nd++)
		{
			fscanf(res5, "%lf", &*(b + i_nd));
		}


		fflush(res5);
		fclose(res5);
	}

}




void mesh(void)
{

	int     i_in, i_nd, i_nd2, i_nd3, i, j, k, k1, ii, jj, kk, kk1, iii, jjj, kkk, check_l1, check_l2, line_c;
	int		dummy, pn1, pn2, pn3, pn4, ap1, ap2, ap3, ap4, ap5, ap6;
	int		tvtx0, tvtx1, tvtx2, tvtx3;

	double  r_length, nor_ax, nor_ay, nor_az, nor_bx, nor_by, nor_bz, nor_mag, ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, abc1, abc2, abc3, abc4, PI, cen_x, cen_y, cen_z;
	double mx, my, mz;
	
	FILE *ff1;
	FILE *ff2;
	FILE *ff3;


	char mesh_data[100];
	char mesh_data2[100];
	char mesh_data3[100];


	PI = 4.0 * atan(1.0);

	sprintf(mesh_data, "point.txt");
	if ((ff2 = fopen(mesh_data, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_in = 1; i_in <= p_num; i_in++)
		{
			fscanf(ff2, "%lf %lf %lf", &*(point_z + i_in), &*(point_y + i_in), &*(point_x + i_in)); //적혈구 돌리기
		}
		/*
		for(i_in = 1; i_in <= p_num; i_in++)
		{

		*(point_x + i_in) /= 30.0;
		*(point_y + i_in) /= 30.0;
		*(point_z + i_in) /= 30.0;
		}
		*/
		for (i_in = 1; i_in <= p_num; i_in++)
		{

			*(point_x + i_in) *= Radius;
			*(point_y + i_in) *= Radius;
			*(point_z + i_in) *= Radius;
		}



		for (i_in = 1; i_in <= p_num; i_in++)
		{

			//*(point_x + i_in) += off_set_x;
			//*(point_y + i_in) += off_set_y;
			//*(point_z + i_in) += off_set_z;
		}
		fflush(ff2);
		fclose(ff2);
	}






	for (i_in = 1; i_in <= f_main; i_in++)
	{
		*(vtx0 + i_in) = 0;
		*(vtx1 + i_in) = 0;
		*(vtx2 + i_in) = 0;

	}

	sprintf(mesh_data2, "face.txt");
	if ((ff1 = fopen(mesh_data2, "r")) == NULL)
	{
		return;
	}
	else
	{
		for (i_in = 1; i_in <= f_main; i_in++)
		{
			fscanf(ff1, "%d %d %d \n", &tvtx0, &tvtx1, &tvtx2);
			*(vtx0 + i_in) = tvtx0;
			*(vtx1 + i_in) = tvtx1;
			*(vtx2 + i_in) = tvtx2;
		}
		fflush(ff1);
		fclose(ff1);
	}




	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////         평면    노멀구하기 평면 점순서정하기 ///////////////////////////////////////////////////////////////////////////////
	for (i_nd = 1; i_nd <= f_main; i_nd++)
	{
		*(n_vect_x + i_nd) = 0.0;
		*(n_vect_y + i_nd) = 0.0;
		*(n_vect_z + i_nd) = 0.0;
	}



	for (i_nd = 1; i_nd <= f_main; i_nd++)
	{
		pn1 = *(vtx0 + i_nd);
		pn2 = *(vtx1 + i_nd);
		pn3 = *(vtx2 + i_nd);


		nor_ax = *(point_x + pn2) - *(point_x + pn1);
		nor_ay = *(point_y + pn2) - *(point_y + pn1);
		nor_az = *(point_z + pn2) - *(point_z + pn1);

		nor_bx = *(point_x + pn3) - *(point_x + pn1);
		nor_by = *(point_y + pn3) - *(point_y + pn1);
		nor_bz = *(point_z + pn3) - *(point_z + pn1);

		*(n_vect_x + i_nd) = nor_ay*nor_bz - nor_az*nor_by;
		*(n_vect_y + i_nd) = nor_az*nor_bx - nor_ax*nor_bz;
		*(n_vect_z + i_nd) = nor_ax*nor_by - nor_ay*nor_bx;

		nor_mag = pow((*(n_vect_x + i_nd)**(n_vect_x + i_nd) + *(n_vect_y + i_nd)**(n_vect_y + i_nd) + *(n_vect_z + i_nd)**(n_vect_z + i_nd)), 0.5);

		*(n_vect_x + i_nd) /= nor_mag;
		*(n_vect_y + i_nd) /= nor_mag;
		*(n_vect_z + i_nd) /= nor_mag;
	}



	for (i_nd = 1; i_nd <= f_main; i_nd++)
	{
		pn1 = *(vtx0 + i_nd);
		pn2 = *(vtx1 + i_nd);
		pn3 = *(vtx2 + i_nd);

		nor_ax = (*(point_x + pn1) + *(point_x + pn2) + *(point_x + pn3)) / 3.0;
		nor_ay = (*(point_y + pn1) + *(point_y + pn2) + *(point_y + pn3)) / 3.0;
		nor_az = (*(point_z + pn1) + *(point_z + pn2) + *(point_z + pn3)) / 3.0;

		nor_bx = ((nor_ax - off_set_x)*(nor_ax - off_set_x) + (nor_ay - off_set_y)*(nor_ay - off_set_y) + (nor_az - off_set_z)*(nor_az - off_set_z));
		nor_by = ((nor_ax + *(n_vect_x + i_nd) - off_set_x)*(nor_ax + *(n_vect_x + i_nd) - off_set_x) + (nor_ay + *(n_vect_y + i_nd) - off_set_y)*(nor_ay + *(n_vect_y + i_nd) - off_set_y) + (nor_az + *(n_vect_z + i_nd) - off_set_z)*(nor_az + *(n_vect_z + i_nd) - off_set_z));

		// 노말이 바깥을 향할때
		if (nor_bx > nor_by)
		{
			*(n_vect_x + i_nd) *= -1.0;
			*(n_vect_y + i_nd) *= -1.0;
			*(n_vect_z + i_nd) *= -1.0;

			*(vtx1 + i_nd) = pn3;
			*(vtx2 + i_nd) = pn2;
		}
	}







	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////     라인을 이루는 점 2개 번호       /////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////     라인을 에 붙어있는 2개 face      /////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////     라인 제외한 face에 의 점 2개       /////////////////////////////////////////////////////////////////////////////////////////////


	for (i_nd = 1; i_nd <= l_main; i_nd++)
	{
		*(line_twopoint1 + i_nd) = 0;
		*(line_twopoint2 + i_nd) = 0;
		*(line_face1 + i_nd) = 0;
		*(line_face2 + i_nd) = 0;
		*(line_verpoint1 + i_nd) = 0;
		*(line_verpoint2 + i_nd) = 0;
		*(line_ref + i_nd) = 0.0;
	}





	i = 0;

	for (i_nd = 1; i_nd <= f_main - 1; i_nd++)
	{
		for (i_nd2 = i_nd + 1; i_nd2 <= f_main; i_nd2++)
		{
			{ k = 0;
			if (*(vtx0 + i_nd) == *(vtx0 + i_nd2)) k += 1;
			if (*(vtx0 + i_nd) == *(vtx1 + i_nd2)) k += 1;
			if (*(vtx0 + i_nd) == *(vtx2 + i_nd2)) k += 1;
			if (*(vtx1 + i_nd) == *(vtx0 + i_nd2)) k += 1;
			if (*(vtx1 + i_nd) == *(vtx1 + i_nd2)) k += 1;
			if (*(vtx1 + i_nd) == *(vtx2 + i_nd2)) k += 1;
			if (*(vtx2 + i_nd) == *(vtx0 + i_nd2)) k += 1;
			if (*(vtx2 + i_nd) == *(vtx1 + i_nd2)) k += 1;
			if (*(vtx2 + i_nd) == *(vtx2 + i_nd2)) k += 1;

			if (k == 2)
			{
				i += 1;

				*(line_face1 + i) = i_nd;
				*(line_face2 + i) = i_nd2;


				if (*(vtx0 + i_nd) == *(vtx0 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx0 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx0 + i_nd); }
				}

				if (*(vtx0 + i_nd) == *(vtx1 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx0 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx0 + i_nd); }
				}
				if (*(vtx0 + i_nd) == *(vtx2 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx0 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx0 + i_nd); }
				}
				if (*(vtx1 + i_nd) == *(vtx0 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx1 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx1 + i_nd); }
				}

				if (*(vtx1 + i_nd) == *(vtx1 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx1 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx1 + i_nd); }
				}
				if (*(vtx1 + i_nd) == *(vtx2 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx1 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx1 + i_nd); }
				}
				if (*(vtx2 + i_nd) == *(vtx0 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx2 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx2 + i_nd); }
				}
				if (*(vtx2 + i_nd) == *(vtx1 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx2 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx2 + i_nd); }

				}
				if (*(vtx2 + i_nd) == *(vtx2 + i_nd2))
				{
					if (*(line_twopoint1 + i) != 0) { *(line_twopoint2 + i) = *(vtx2 + i_nd); }
					if (*(line_twopoint1 + i) == 0) { *(line_twopoint1 + i) = *(vtx2 + i_nd); }
				}

				if ((*(vtx0 + i_nd) != *(line_twopoint1 + i)) && (*(vtx0 + i_nd) != *(line_twopoint2 + i)))
				{
					*(line_verpoint1 + i) = *(vtx0 + i_nd);
				}
				if ((*(vtx1 + i_nd) != *(line_twopoint1 + i)) && (*(vtx1 + i_nd) != *(line_twopoint2 + i)))
				{
					*(line_verpoint1 + i) = *(vtx1 + i_nd);
				}
				if ((*(vtx2 + i_nd) != *(line_twopoint1 + i)) && (*(vtx2 + i_nd) != *(line_twopoint2 + i)))
				{
					*(line_verpoint1 + i) = *(vtx2 + i_nd);
				}

				if ((*(vtx0 + i_nd2) != *(line_twopoint1 + i)) && (*(vtx0 + i_nd2) != *(line_twopoint2 + i)))
				{
					*(line_verpoint2 + i) = *(vtx0 + i_nd2);
				}
				if ((*(vtx1 + i_nd2) != *(line_twopoint1 + i)) && (*(vtx1 + i_nd2) != *(line_twopoint2 + i)))
				{
					*(line_verpoint2 + i) = *(vtx1 + i_nd2);
				}
				if ((*(vtx2 + i_nd2) != *(line_twopoint1 + i)) && (*(vtx2 + i_nd2) != *(line_twopoint2 + i)))
				{
					*(line_verpoint2 + i) = *(vtx2 + i_nd2);
				}

			}

			}

		}
	}


	///    i 프린트 했을때 main_l이 나와야함.



	//				ff3 = fopen("process.txt", "a");
	//	fprintf(ff3,"hello%d\n", i);
	//	fflush(ff3);
	//	fclose(ff3);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////    라인 ref정하기            /////////////////////////////////////////////////////////////////////////////////////////////
	for (i_nd = 1; i_nd <= l_main; i_nd++)
	{
		*(line_ref + i_nd) = 0.0;
	}


	for (i_nd = 1; i_nd <= l_main; i_nd++)
	{
		kk = *(line_twopoint1 + i_nd);
		kk1 = *(line_twopoint2 + i_nd);
		r_length = pow(*(point_x + kk1) - *(point_x + kk), 2) + pow(*(point_y + kk1) - *(point_y + kk), 2) + pow(*(point_z + kk1) - *(point_z + kk), 2);
		*(line_ref + i_nd) = pow(r_length, 0.5);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////ref angle////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	printf("test2\n");



	printf("test3\n");







	for (i_nd = 1; i_nd <= l_main; i_nd++)
	{


		pn1 = *(line_face1 + i_nd);
		pn2 = *(line_face2 + i_nd);

		pn3 = *(line_verpoint1 + i_nd);
		pn4 = *(line_verpoint2 + i_nd);

		ap1 = *(vtx0 + pn1);
		ap2 = *(vtx1 + pn1);
		ap3 = *(vtx2 + pn1);

		ap4 = *(vtx0 + pn2);
		ap5 = *(vtx1 + pn2);
		ap6 = *(vtx2 + pn2);



		mx = (*(point_x + ap4) + *(point_x + ap5) + *(point_x + ap6)) / 3.0 - (*(point_x + ap1) + *(point_x + ap2) + *(point_x + ap3)) / 3.0;
		my = (*(point_y + ap4) + *(point_y + ap5) + *(point_y + ap6)) / 3.0 - (*(point_y + ap1) + *(point_y + ap2) + *(point_y + ap3)) / 3.0;
		mz = (*(point_z + ap4) + *(point_z + ap5) + *(point_z + ap6)) / 3.0 - (*(point_z + ap1) + *(point_z + ap2) + *(point_z + ap3)) / 3.0;

		abc3 = (mx * *(n_vect_x + pn1) + my * *(n_vect_y + pn1) + mz * *(n_vect_z + pn1));




		if (abc3 <= 0.0)
		{
			ax1 = *(n_vect_x + pn1) * *(n_vect_x + pn2);
			ay1 = *(n_vect_y + pn1) * *(n_vect_y + pn2);
			az1 = *(n_vect_z + pn1) * *(n_vect_z + pn2);
			abc1 = ax1 + ay1 + az1;
			abc2 = PI - acos(abc1);

		}


		if (abc3 > 0.0)
		{
			ax1 = *(n_vect_x + pn1) * *(n_vect_x + pn2);
			ay1 = *(n_vect_y + pn1) * *(n_vect_y + pn2);
			az1 = *(n_vect_z + pn1) * *(n_vect_z + pn2);
			abc1 = ax1 + ay1 + az1;
			abc2 = PI + acos(abc1);

		}




		*(r_angle + i_nd) = abc2;
	}




	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////ref area////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////ref vol////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cen_x = cen_y = cen_z = 0.0;

	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{
		cen_x += *(point_x + i_nd);
		cen_y += *(point_y + i_nd);
		cen_z += *(point_z + i_nd);
	}

	cen_x /= double(p_num);
	cen_y /= double(p_num);
	cen_z /= double(p_num);



	/*
	for(i_nd = 1; i_nd <= f_main; i_nd++)
	{
	pn1 = *(vtx0 + i_nd);
	pn2 = *(vtx1 + i_nd);
	pn3 = *(vtx2 + i_nd);

	ax1 = *(point_x + pn2)  - *(point_x + pn1) ;
	ay1 = *(point_y + pn2)  - *(point_y + pn1) ;
	az1 = *(point_z + pn2)  - *(point_z + pn1) ;


	ax2 = *(point_x + pn3)  - *(point_x + pn1) ;
	ay2 = *(point_y + pn3)  - *(point_y + pn1) ;
	az2 = *(point_z + pn3)  - *(point_z + pn1) ;



	ax3 = cen_x  - *(point_x + pn1) ;
	ay3 = cen_y  - *(point_y + pn1) ;
	az3 = cen_z  - *(point_z + pn1) ;



	*(area_r + i_nd) = (pow(((ay1*az2-ay2*az1)*(ay1*az2-ay2*az1) + (ax2*az1 - ax1*az2)*(ax2*az1 - ax1*az2) + (ax1*ay2-ay1*ax2)*(ax1*ay2-ay1*ax2)),0.5))/2.0;
	*(vol_r + i_nd) = (pow(((ay1*az2-ay2*az1) * ax3 + (ax2*az1 - ax1*az2) * ay3 + (ax1*ay2-ay1*ax2) * az3)  *  ((ay1*az2-ay2*az1) * ax3 + (ax2*az1 - ax1*az2) * ay3 + (ax1*ay2-ay1*ax2) * az3 ),0.5))/6.0;



	/////////////////////////////////                       배열로 만들어서 RBC 각각 선언해야함
	r_tot_area = 0.0;
	r_tot_vol  = 0.0;

	for(i_nd = 1; i_nd <= f_main; i_nd++)
	{

	r_tot_area += *(area_r + i_nd) ;
	r_tot_vol += *(vol_r + i_nd) ;

	}








	ff3 = fopen("process.txt", "a");
	fprintf(ff3,"hello%lf\n%lf\n", r_tot_area, r_tot_vol);
	fflush(ff3);
	fclose(ff3);

	//!!병렬화할때 이부분 꼭 선언
	*/
	for (i_nd = 0; i_nd <1; i_nd++)
	{
		*(vtx0 + i_nd) = *(vtx0 + f_main);
		*(vtx1 + i_nd) = *(vtx1 + f_main);
		*(vtx2 + i_nd) = *(vtx2 + f_main);

		*(n_vect_x + i_nd) = *(n_vect_x + (f_main));
		*(n_vect_y + i_nd) = *(n_vect_y + (f_main));
		*(n_vect_z + i_nd) = *(n_vect_z + (f_main));
		//					*(area_r + i_nd) = *(area_r + (f_main));
		//*(vol_r + i_nd) = *(vol_r + (f_main)) ;
	}


	for (i_nd = 0; i_nd <1; i_nd++)
	{
		*(line_twopoint1 + i_nd) = *(line_twopoint1 + l_main);
		*(line_twopoint2 + i_nd) = *(line_twopoint2 + l_main);
		*(line_face1 + i_nd) = *(line_face1 + l_main);
		*(line_face2 + i_nd) = *(line_face2 + l_main);
		*(line_verpoint1 + i_nd) = *(line_verpoint1 + l_main);
		*(line_verpoint2 + i_nd) = *(line_verpoint2 + l_main);
		*(line_ref + i_nd) = *(line_ref + l_main);
		*(r_angle + i_nd) = *(r_angle + l_main);



	}


	for (i_nd = 0; i_nd <1; i_nd++)
	{

		*(point_x + i_nd) = *(point_x + p_num);
		*(point_y + i_nd) = *(point_y + p_num);
		*(point_z + i_nd) = *(point_z + p_num);


	}


	//for(i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	//for(i_nd = (i_nd2-1)*f_main+1; i_nd <= f_main*i_nd2; i_nd++)

	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
	{


		if (i_nd2 != 1)
		{
			*(vtx0 + i_nd) = (i_nd2 - 1) * p_num + *(vtx0 + (i_nd % f_main));
			*(vtx1 + i_nd) = (i_nd2 - 1) * p_num + *(vtx1 + (i_nd % f_main));
			*(vtx2 + i_nd) = (i_nd2 - 1) * p_num + *(vtx2 + (i_nd % f_main));

			*(n_vect_x + i_nd) = *(n_vect_x + (i_nd % f_main));
			*(n_vect_y + i_nd) = *(n_vect_y + (i_nd % f_main));
			*(n_vect_z + i_nd) = *(n_vect_z + (i_nd % f_main));
			//*(area_r + i_nd) = *(area_r + (i_nd % f_main));
			//   *(vol_r + i_nd) = *(vol_r + (i_nd % f_main)) ;


		}
	}


	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	for (i_nd = ((i_nd2 - 1)*l_main + 1); i_nd <= l_main*i_nd2; i_nd++)
	{

		if (i_nd2 != 1)
		{
			*(line_twopoint1 + i_nd) = (i_nd2 - 1) * p_num + *(line_twopoint1 + (i_nd % l_main));
			*(line_twopoint2 + i_nd) = (i_nd2 - 1) * p_num + *(line_twopoint2 + (i_nd % l_main));
			*(line_face1 + i_nd) = (i_nd2 - 1) * f_main + *(line_face1 + (i_nd % l_main));
			*(line_face2 + i_nd) = (i_nd2 - 1) * f_main + *(line_face2 + (i_nd % l_main));
			*(line_verpoint1 + i_nd) = (i_nd2 - 1) * p_num + *(line_verpoint1 + (i_nd % l_main));
			*(line_verpoint2 + i_nd) = (i_nd2 - 1) * p_num + *(line_verpoint2 + (i_nd % l_main));
			*(line_ref + i_nd) = *(line_ref + (i_nd % l_main));
			*(r_angle + i_nd) = *(r_angle + (i_nd % l_main));
		}
	}




	//////////////////////////////////////////////////////// 
	// setparticle 배열


	int rbcnx, rbcnz, rbcny, rbcindex;
	double a_dis, a_angle, a_angle2;
	double unit_y, unit_z;
	//int iii;

	rbcnz = 1;//10; //쌓는방향
	rbcnx = 1;//rbcn / rbcnz;
	rbcny = 1;
	rbcindex=0;
	//setparticle2
	//for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	for (i_nd = 1; i_nd <= rbcn; i_nd++) 
	{
		rbcindex ++;
		cen_rx[rbcindex] = XLENGTH/2+1;
		cen_ry[rbcindex] = YLENGTH/2+1;
		cen_rz[rbcindex] = ZLENGTH/2+1;
	}


	for (i_nd2 = 1; i_nd2 <= 1; i_nd2++)
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		*(point_x + i_nd) += cen_rx[i_nd2];
		*(point_y + i_nd) += cen_ry[i_nd2];
		*(point_z + i_nd) += cen_rz[i_nd2];
	}




	
	for(iii=1; iii<=p_num; iii++)
	{
		i_nd2=2;
		i_nd = p_num*(i_nd2-1)+iii;

		a_dis = 0;
		a_dis += (*(point_y + iii)-cen_ry[i_nd2])*(*(point_y + iii)-cen_ry[i_nd2]);
		a_dis += (*(point_z + iii)-cen_rz[i_nd2])*(*(point_z + iii)-cen_rz[i_nd2]);
		a_dis *=10000;
		a_dis = pow(a_dis, 0.5);
		a_dis /=100;

		unit_y = (*(point_y + iii)-cen_ry[i_nd2])/a_dis;
		unit_z = (*(point_z + iii)-cen_rz[i_nd2])/a_dis;
		//if (unit_z==0) a_angle = 0;
		a_angle = atan(unit_y/unit_z);

		a_angle2 = a_angle + (145*PI2/180);


		*(point_x + i_nd) = *(point_x + iii);
		*(point_y + i_nd) += a_dis * cos(a_angle2)+cen_ry[i_nd2];
		*(point_z + i_nd) += a_dis * sin(a_angle2)+cen_rz[i_nd2];

	}

	for(iii=1; iii<=p_num; iii++)
	{
		i_nd2=3;
		i_nd = p_num*(i_nd2-1)+iii;

		a_dis = 0;
		a_dis += (*(point_y + iii)-cen_ry[i_nd2])*(*(point_y + iii)-cen_ry[i_nd2]);
		a_dis += (*(point_z + iii)-cen_rz[i_nd2])*(*(point_z + iii)-cen_rz[i_nd2]);
		a_dis *=10000;
		a_dis = pow(a_dis, 0.5);
		a_dis /=100;

		unit_y = (*(point_y + iii)-cen_ry[i_nd2])/a_dis;
		unit_z = (*(point_z + iii)-cen_rz[i_nd2])/a_dis;
		//if (unit_z==0) a_angle = 0;
		a_angle = atan(unit_y/unit_z);

		a_angle2 = a_angle + (265*PI2/180);


		*(point_x + i_nd) = *(point_x + iii);
		*(point_y + i_nd) += a_dis * cos(a_angle2)+cen_ry[i_nd2];
		*(point_z + i_nd) += a_dis * sin(a_angle2)+cen_rz[i_nd2];

	}
	

	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	for(iii=1; iii<=p_num; iii++)
	{
		i_nd = p_num*(i_nd2-1)+iii;
			*(point_xi + i_nd) = *(point_x + i_nd);
			*(point_yi + i_nd) = *(point_y + i_nd);
			*(point_zi + i_nd) = *(point_z + i_nd);
	}




	return;
}




void   intial_viscosity(void)
{
	int    pn1, pn2, pn3, i_nd, i_nd2, i_nd3, x, y, z, x1, y1, z1, x2, y2, z2;
	double  nor_ax, nor_ay, nor_az, nor_bx, nor_by, xx, yy, zz, xxx, yyy, zzz;
	double  cen_x, cen_y, cen_z, cen_x2, cen_y2, cen_z2;
	
	cen_x = cen_y = cen_z = 0.0;

	

	for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
			for (z = 1; z <= ZLENGTH; z++)
			{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				*(inout_vis + i_nd) = outvis;
				*(vis_check+ i_nd) = 1;

			}



	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	{
		for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
		{
			cen_x += *(point_x + i_nd);
			cen_y += *(point_y + i_nd);
			cen_z += *(point_z + i_nd);
		}
		cen_x /= p_num;
		cen_y /= p_num;
		cen_z /= p_num;


		xx = cen_x + Radius+1.0;
		yy = cen_y + Radius+1.0;
		zz = cen_z + Radius+1.0;

		xxx = cen_x - Radius-1.0;
		yyy = cen_y - Radius-1.0;
		zzz = cen_z - Radius-1.0;

		x1 = abs(xxx);
		x2 = abs(xx);

		y1 = abs(yyy);
		y2 = abs(yy);

		z1 = abs(zzz);
		z2 = abs(zz);

		if (x1 < 0) { x1 = 0; }
		if (y1 < 0) { y1 = 0; }
		if (z1 < 0) { z1 = 0; }
		if (x2 > XLENGTH) { x2 = XLENGTH; }
		if (y2 > YLENGTH) { y2 = YLENGTH; }
		if (z2 > ZLENGTH) { z2 = ZLENGTH; }

		

		for (z = z1; z <= z2; z++)
		{
			for (y = y1; y <= y2; y++)
			{
				for (x = x1; x <= x2; x++)
				{
					i_nd3 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


					if (((x - cen_x) * (x - cen_x) + (y - cen_y) * (y - cen_y) + (z - cen_z) * (z - cen_z)) < Radius*Radius)
					{
						*(inout_vis + i_nd3) = invis;
						*(vis_check + i_nd3) = 2;
					}


				}
			}
		}


		
	}



	return;
}



void    viscosity(void)
{
	int    pn1, pn2, pn3, i_nd, i_nd2, i_nd3, x, y, z, x1, y1, z1, x2, y2, z2;
	double  nor_ax, nor_ay, nor_az, nor_bx, nor_by, xx, yy, zz, xxx, yyy, zzz;
	double  cen_x, cen_y, cen_z, cen_x2, cen_y2, cen_z2;

	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	{

		for (i_nd = ((i_nd2 - 1)*f_main + 1); i_nd <= f_main*i_nd2; i_nd++)
		{

			cen_x = (*(point_x + pn1) + *(point_x + pn2) + *(point_x + pn3)) / 3.0;
			cen_y = (*(point_y + pn1) + *(point_y + pn2) + *(point_y + pn3)) / 3.0;
			cen_z = (*(point_z + pn1) + *(point_z + pn2) + *(point_z + pn3)) / 3.0;

			xx = cen_x + vis_c;
			yy = cen_y + vis_c;
			zz = cen_z + vis_c;

			xxx = cen_x - vis_c;
			yyy = cen_y - vis_c;
			zzz = cen_z - vis_c;

			x1 = abs(xxx);
			x2 = abs(xx);

			y1 = abs(yyy);
			y2 = abs(yy);

			z1 = abs(zzz);
			z2 = abs(zz);

			if (x1 < 0) { x1 = 0; }
			if (y1 < 0) { y1 = 0; }
			if (z1 < 0) { z1 = 0; }
			if (x2 > XLENGTH) { x2 = XLENGTH; }
			if (y2 > YLENGTH) { y2 = YLENGTH; }
			if (z2 > ZLENGTH) { z2 = ZLENGTH; }


			for (z = z1; z <= z2; z++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (x = x1; x <= x2; x++)
					{
						i_nd3 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						if (((x - (cen_x + 0.1**(n_vect_x + i_nd))) * (x - (cen_x + 0.1**(n_vect_x + i_nd))) + (y - (cen_y + 0.1**(n_vect_y + i_nd))) * (y - (cen_y + 0.1**(n_vect_y + i_nd))) + (z - (cen_z + 0.1**(n_vect_z + i_nd))) * (z - (cen_z + 0.1**(n_vect_z + i_nd)))) > ((x - cen_x) * (x - cen_x) + (y - cen_y) * (y - cen_y) + (z - cen_z) * (z - cen_z)))
						{
							*(inout_vis + i_nd3) = invis;
							*(vis_check + i_nd3) = 2;

						}
						else
						{ 
							*(inout_vis + i_nd3) = outvis;
							*(vis_check + i_nd3) = 1;
						
						}

					}
				}
			}


		}
	}




	return;
}



// slide
/*
void    density_data (void)
{   int     x, y, z, i_nd;
	double  ux1, uy1, uz1, presw;	
    FILE    *f1;
    FILE    *f2;
    char FileNameVelocity[100];
    char string[100] = {'0'};

    sprintf(FileNameVelocity, "data2_%d.plt",Timek);




 //  Data for Tecplot

   f1 = fopen("movie.dat", "a");
   fprintf(f1,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH,YLENGTH, ZLENGTH);
      
       f2 = fopen(FileNameVelocity, "w");
   fprintf(f2,"Zone I=   %d,J=   %d,F=POINT\n", XLENGTH, YLENGTH);


   for (z = (1 + ZLENGTH) / 2; z <= (1 + ZLENGTH) / 2; z++)
  {
   for ( y = 1; y <= YLENGTH; y++ )
      {
        for ( x = 1; x <= XLENGTH; x++ )
         {
			                  
                  {
         i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;      
		   


		 ux1 = *(ux + i_nd);
		 uy1 = *(uy + i_nd);
		 uz1 = *(uz + i_nd);
		 presw = *(press + i_nd); 
		
		 
		 if(*(bnode+i_nd) != 0)
			{
		     ux1 = 0.0;
			 uy1 = 0.0;
			 uz1 = 0.0;
			 presw = 1.0/3.0;
			 }

//		 fprintf(f1, "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, z, vis_check[i_nd], ux1, uy1, uz1, *(FforceX + i_nd), *(FforceY + i_nd), *(FforceZ + i_nd));
		 //fprintf(f1, "\n");


		 fprintf(f2, "%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, *(bnode + i_nd), ux1, uy1, uz1, *(FforceX + i_nd), *(FforceY + i_nd), *(FforceZ + i_nd));
		 fprintf(f2, "\n");


//         fprintf(f1,"%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,z,vis_check[i_nd],ux1,uy1,uz1,*(FforceX + i_nd),*(FforceY + i_nd),*(FforceZ + i_nd));
//         fprintf(f1, "\n");

		 
		 //fprintf(f2, "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, z, vis_check[i_nd], ux1, uy1, uz1, *(FforceX + i_nd), *(FforceY + i_nd), *(FforceZ + i_nd));
         //fprintf(f2, "\n");

			 }
       }
   }
   }


   

   
    fflush(f1);
    fclose(f1);
        fflush(f2);
    fclose(f2);


    return;
}
*/


void    density_data(void)
{
	int     x, y, z, i_nd;
	double  ux1, uy1, uz1, presw;
	FILE    *f1;
	FILE    *f2;
	char FileNameVelocity[100];
	char string[100] = { '0' };

	sprintf(FileNameVelocity, "data2_%d.plt", Timek);




	//  Data for Tecplot

	f1 = fopen("movie.dat", "a");
	fprintf(f1, "Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH, YLENGTH, ZLENGTH);

	f2 = fopen(FileNameVelocity, "w");
	fprintf(f2, "Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH, YLENGTH, ZLENGTH);


	for (z = 1; z <= ZLENGTH; z++)
	{
		for (y = 1; y <= YLENGTH; y++)
		{
			for (x = 1; x <= XLENGTH; x++)
			{

				{
					i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;



					ux1 = *(ux + i_nd);
					uy1 = *(uy + i_nd);
					uz1 = *(uz + i_nd);
					presw = *(press + i_nd);


					if (*(bnode + i_nd) != 0)
					{
						ux1 = 0.0;
						uy1 = 0.0;
						uz1 = 0.0;
						presw = 1.0 / 3.0;
					}

					//		 fprintf(f1, "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, z, vis_check[i_nd], ux1, uy1, uz1, *(FforceX + i_nd), *(FforceY + i_nd), *(FforceZ + i_nd));
					//fprintf(f1, "\n");


					fprintf(f2, "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, z, *(bnode + i_nd), ux1, uy1, uz1, presw, *(FforceX + i_nd), *(FforceY + i_nd), *(FforceZ + i_nd));
					fprintf(f2, "\n");


					//         fprintf(f1,"%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,z,vis_check[i_nd],ux1,uy1,uz1,*(FforceX + i_nd),*(FforceY + i_nd),*(FforceZ + i_nd));
					//         fprintf(f1, "\n");


					//fprintf(f2, "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, z, vis_check[i_nd], ux1, uy1, uz1, *(FforceX + i_nd), *(FforceY + i_nd), *(FforceZ + i_nd));
					//fprintf(f2, "\n");

				}
			}
		}
	}





	fflush(f1);
	fclose(f1);
	fflush(f2);
	fclose(f2);


	return;
}



void    restart_write1(void)
{
	int     i_nd;


	FILE    	*f7;
	char 		FileNamerestart[10] = "restart";
	char 		string[10] = { '0' };
	sprintf(FileNamerestart, "restart_1_%d.txt", Timek);



	f7 = fopen(FileNamerestart, "w");


	for (i_nd = 1; i_nd <= tf_main; i_nd++)
	{
		fprintf(f7, "%d\t%d\t%d\t", *(vtx0 + i_nd), *(vtx1 + i_nd), *(vtx2 + i_nd));
	}

	fflush(f7);
	fclose(f7);
}

void    restart_write2(void)
{
	int     i_nd;


	FILE    	*f8;
	char 		FileNamerestart[10] = "restart";
	char 		string[10] = { '0' };
	sprintf(FileNamerestart, "restart_2_%d.txt", Timek);



	f8 = fopen(FileNamerestart, "w");


	for (i_nd = 1; i_nd <= tp_num; i_nd++)
	{
		fprintf(f8, "%lf\t%lf\t%lf\t", *(point_x + i_nd), *(point_y + i_nd), *(point_z + i_nd));
	}

	fflush(f8);
	fclose(f8);
}

void    restart_write3(void)
{
	int     i_nd;


	FILE    	*f9;
	char 		FileNamerestart[10] = "restart";
	char 		string[10] = { '0' };
	sprintf(FileNamerestart, "restart_3_%d.txt", Timek);



	f9 = fopen(FileNamerestart, "w");


	for (i_nd = 1; i_nd <= tl_main; i_nd++)
	{
		fprintf(f9, "%lf\t%lf\t", *(line_ref + i_nd), *(r_angle + i_nd));
	}

	fflush(f9);
	fclose(f9);
}


void    restart_write4(void)
{
	int     i_nd;


	FILE    	*f10;
	char 		FileNamerestart[10] = "restart";
	char 		string[10] = { '0' };
	sprintf(FileNamerestart, "restart_4_%d.txt", Timek);



	f10 = fopen(FileNamerestart, "w");


	for (i_nd = 1; i_nd <= tl_main; i_nd++)
	{
		fprintf(f10, "%d\t%d\t%d\t%d\t%d\t%d\t", *(line_face1 + i_nd), *(line_face2 + i_nd), *(line_verpoint1 + i_nd), *(line_verpoint2 + i_nd), *(line_twopoint1 + i_nd), *(line_twopoint2 + i_nd));
	}

	fflush(f10);
	fclose(f10);
}

void    restart_write5(void)
{
	int     i_nd, x, y, z, i;


	FILE    	*f11;
	char 		FileNamerestart[10] = "restart";
	char 		string[10] = { '0' };
	sprintf(FileNamerestart, "restart_5_%d.txt", Timek);



	f11 = fopen(FileNamerestart, "w");

	for (i_nd = 1; i_nd <= i_main2; i_nd++)
	{
		fprintf(f11, "%lf\t", *(b + i_nd));
	}

	fflush(f11);
	fclose(f11);
}


/*
void    initial_density_data (void)
{   int     x, y, z, i_nd;
	double  ux1, uy1, uz1, presw;	
    FILE    *f4;
    
    char FileNamePP[100];
    char string[100] = {'0'};

    sprintf(FileNamePP, "initial_data2_%d.plt",k);


 //  Data for Tecplot

      
       f4 = fopen(FileNamePP, "w");
   fprintf(f4,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH,YLENGTH, ZLENGTH);

 for (z = 1; z <= ZLENGTH; z++)
  {
   for ( y = 1; y <= YLENGTH; y++ )
      {
        for ( x = 1; x <= XLENGTH; x++ )
         {
			                  
                  {
         
		      i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;    


		 ux1 = *(ux + i_nd);
		 uy1 = *(uy + i_nd);
		 uz1 = *(uz + i_nd);
		 presw = *(press + i_nd); 
		
		 
		 if(*(bnode+i_nd) != 0)
			{
		     ux1 = 0.0;
			 uy1 = 0.0;
			 uz1 = 0.0;
			 presw = 1.0/3.0;
			 }

         		 
         fprintf(f4,"%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n",x,y,z,*(bnode + i_nd),ux1,uy1,uz1,presw);
         fprintf(f4, "\n");

			 }
       }
   }
   }



   
   
    fflush(f4);
    fclose(f4);
 

    return;
}
*/
/*
void pressure_data()
{
	FILE *fP0, *fP1, *fP2;
	int upoint0, upoint1, upoint2, U;

	//i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
	//upoint0 = XLENGTH * (YLENGTH/2) + 1; // inlet
	//upoint1 = XLENGTH * (YLENGTH/2) + (XLENGTH/2); //midpoint
	//upoint2 = XLENGTH * (YLENGTH/2) + (XLENGTH-1); //outlet

	upoint0 = (ZLENGTH/2) + ZLENGTH*(YLENGTH/2 + YLENGTH*(1)) +1; //inlet
	upoint1 = (ZLENGTH/2) + ZLENGTH*(YLENGTH/2 + YLENGTH*(XLENGTH/2)) +1;
	upoint2 = (ZLENGTH/2) + ZLENGTH*(YLENGTH/2 + YLENGTH*(XLENGTH-1)) +1;

	U= pow((pow(*(ux + upoint0),2.0)) + (pow(*(uy + upoint0),2.0)) + (pow(*(uz + upoint0),2.0)),0.5);
	fP0 = fopen("inlet_flowwrate_rho.dat", "a");
	fprintf(fP0,"%d   %f   %f\n", k, U, *(rho + upoint0));
	fflush(fP0);
	fclose(fP0);
	
	U= pow((pow(*(ux + upoint1),2.0)) + (pow(*(uy + upoint1),2.0)) + (pow(*(uz + upoint1),2.0)),0.5);
	fP1 = fopen("middle_flowwrate_rho.dat", "a");
	fprintf(fP1,"%d   %f    %f\n", k, U, *(rho + upoint1)); 
	fflush(fP1);
	fclose(fP1);

	U= pow((pow(*(ux + upoint2),2.0)) + (pow(*(uy + upoint2),2.0)) + (pow(*(uz + upoint2),2.0)),0.5);
	fP2 = fopen("outlet_Z_flowwrate_rho.dat", "a");
	fprintf(fP2, "%d %f %f %f\n", k, Z[k%N], U, *(rho + upoint2));
	fflush(fP2);
	fclose(fP2);
	
	return;

}
*/

void    mesh_nor (void)
{   int     x, y, z, i_nd,pn1,pn2,pn3;
	double  f_x1, f_y1, f_z1, px, py, pz;	

    FILE    *f4;
    char FileNamemesh[100];
    char string[100] = {'0'};

    sprintf(FileNamemesh, "adata_normal_%d.plt",Timek);


 //  Data for Tecplot

   f4 = fopen(FileNamemesh, "w");
    fprintf(f4,"VARIABLES = ""X"", ""Y"", ""Z"", ""u"", ""v"", ""w"" ");
	fprintf(f4, "\n");
	fprintf(f4,"ZONE  I = %d, DATAPACKING=POINT",tf_main);
	fprintf(f4, "\n");
   //fprintf(f4,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH,YLENGTH, ZLENGTH);



      	  for(i_nd = 1; i_nd <= tf_main; i_nd++)
		  {		  
		    pn1 = *(vtx0 + i_nd);
			pn2 = *(vtx1 + i_nd);
			pn3 = *(vtx2 + i_nd);

			  f_x1 = (*(point_x + pn1) + *(point_x + pn2) + *(point_x + pn3)) / 3.0;
			  f_y1 = (*(point_y + pn1) + *(point_y + pn2) + *(point_y + pn3)) / 3.0;
			  f_z1 = (*(point_z + pn1) + *(point_z + pn2) + *(point_z + pn3)) / 3.0;
		


			  px = fmod(f_x1, double(XLENGTH));
			  py = fmod(f_y1, double(YLENGTH));
			  pz = fmod(f_z1, double(ZLENGTH));

			  if (px < 0) { px = px + double(XLENGTH); }
			  if (py < 0) { py = py + double(YLENGTH); }
			  if (pz < 0) { pz = pz + double(ZLENGTH); }
			  fprintf(f4,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",px,py,pz,*(n_vect_x + i_nd), *(n_vect_y + i_nd), *(n_vect_z + i_nd));
              fprintf(f4, "\n");
				  




			

		  }




    fflush(f4);
    fclose(f4);


    return;
}

void    mesh_com(void)
{
	int     x, y, z, i_nd, i_nd2, pn1, pn2, pn3, elm;
	double  f_x1, f_y1, f_z1, px, py, pz, px1, py1, pz1, px2, py2, pz2, axi, transx, transy, axi2, transy2, transz, transz2;

	FILE    *f3;
	FILE    *ff7;
	
	char FileNamemesh[100];
	char string[100] = { '0' };

	sprintf(FileNamemesh, "data2_%d_1.plt", Timek);


	//  Data for Tecplot

	elm = tf_main;



#pragma omp parallel
	{


#pragma omp	for private(i_nd)
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		{


			//좌표 봐뀌는 지점
			if (*(point_x + ((i_nd2 - 1)*p_num + 1)) < (1.0 - Radius))
			{
				for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
				{

					*(point_x + i_nd) += XLENGTH;

				}


			}


			if (*(point_x + ((i_nd2 - 1)*p_num + 1)) > (XLENGTH + Radius))
			{
				for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
				{

					*(point_x + i_nd) -= XLENGTH;

				}


			}


			if (*(point_y + ((i_nd2 - 1)*p_num + 1)) < (1.0 - Radius))
			{
				for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
				{

					*(point_y + i_nd) += YLENGTH;

				}


			}


			if (*(point_y + ((i_nd2 - 1)*p_num + 1)) > (YLENGTH + Radius))
			{
				for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
				{

					*(point_y + i_nd) -= YLENGTH;

				}


			}



		}




	}






	f3 = fopen(FileNamemesh, "w");
	//fprintf(f3, "VARIABLES = ""X"", ""Y"", ""Z"", ""bFx"", ""bFy"", ""bFz"", ""sFx"", ""sFy"", ""sFz"" ");
	fprintf(f3, "VARIABLES = ""V1"", ""V2"", ""V3"" ""A1"" ""A2"" ""A3"" ");
	fprintf(f3, "\n");
	fprintf(f3, "ZONE N = %d, E = %d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE", tp_num*fakex, tf_main*fakex);
	fprintf(f3, "\n");
	//fprintf(f4,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH,YLENGTH, ZLENGTH);


	for (i_nd2 = 1; i_nd2 <= fakex; i_nd2++)
	for (i_nd = 1; i_nd <= tp_num; i_nd++)
	{
		//fprintf(f3, "%lf\t%lf\t%lf\t%lf\n", *(point_x + i_nd) + (i_nd2 - 1)*XLENGTH, *(point_y + i_nd), *(point_z + i_nd), *(Deff + i_nd));
		//fprintf(f3, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", *(point_x + i_nd), *(point_y + i_nd), *(point_z + i_nd), *(bforceX + i_nd), *(bforceY + i_nd), *(bforceZ + i_nd), *(sforceX + i_nd), *(sforceY + i_nd), *(sforceZ + i_nd));
		fprintf(f3, "%lf\t%lf\t%lf\t%.12lf\t%.12lf\t%.12lf\n", *(point_x + i_nd) + (i_nd2 - 1)*XLENGTH, *(point_y + i_nd), *(point_z + i_nd), um[i_nd][0],um[i_nd][1],um[i_nd][2] );
		fprintf(f3, "\n");

	}



	for (i_nd2 = 1; i_nd2 <= fakex; i_nd2++)
	for (i_nd = 1; i_nd <= tf_main; i_nd++)
	{

		fprintf(f3, "%d\t%d\t%d\n", *(vtx0 + i_nd) + (i_nd2 - 1) * tp_num, *(vtx1 + i_nd) + (i_nd2 - 1) * tp_num, *(vtx2 + i_nd) + (i_nd2 - 1) * tp_num);
			fprintf(f3, "\n");
	


	}





	/*
	for (i_nd = 1; i_nd <= tf_main; i_nd++)
	{

		*(vtx3 + i_nd) = 1;

		px = fmod(*(point_x + *(vtx0 + i_nd)), double(XLENGTH));
		py = fmod(*(point_y + *(vtx0 + i_nd)), double(YLENGTH));
		pz = fmod(*(point_z + *(vtx0 + i_nd)), double(ZLENGTH));
		if (px < 0) { px = px + double(XLENGTH); }
		if (py < 0) { py = py + double(YLENGTH); }
		if (pz < 0) { pz = pz + double(ZLENGTH); }

		px1 = fmod(*(point_x + *(vtx1 + i_nd)), double(XLENGTH));
		py1 = fmod(*(point_y + *(vtx1 + i_nd)), double(YLENGTH));
		pz1 = fmod(*(point_z + *(vtx1 + i_nd)), double(ZLENGTH));
		if (px1 < 0) { px1 = px1 + double(XLENGTH); }
		if (py1 < 0) { py1 = py1 + double(YLENGTH); }
		if (pz1 < 0) { pz1 = pz1 + double(ZLENGTH); }

		px2 = fmod(*(point_x + *(vtx2 + i_nd)), double(XLENGTH));
		py2 = fmod(*(point_y + *(vtx2 + i_nd)), double(YLENGTH));
		pz2 = fmod(*(point_z + *(vtx2 + i_nd)), double(ZLENGTH));
		if (px2 < 0) { px2 = px2 + double(XLENGTH); }
		if (py2 < 0) { py2 = py2 + double(YLENGTH); }
		if (pz2 < 0) { pz2 = pz2 + double(ZLENGTH); }



		if ((fabs(px - px1) > 10) || (fabs(px - px2) > 10) || (fabs(px2 - px1) > 10) || (fabs(py - py1) > 10) || (fabs(py - py2) > 10) || (fabs(py2 - py1) > 10) || (fabs(pz - pz1) > 10) || (fabs(pz - pz2) > 10) || (fabs(pz2 - pz1) > 10))
		{
			elm -= 1;
			*(vtx3 + i_nd) = 0;
		}







	}



	f3 = fopen(FileNamemesh, "w");
	fprintf(f3, "VARIABLES = ""X"", ""Y"", ""Z"", ""bFx"", ""bFy"", ""bFz"", ""sFx"", ""sFy"", ""sFz"" ");
	fprintf(f3, "\n");
	fprintf(f3, "ZONE N = %d, E = %d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE", tp_num, elm);
	fprintf(f3, "\n");
	//fprintf(f4,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH,YLENGTH, ZLENGTH);



	for (i_nd = 1; i_nd <= tp_num; i_nd++)
	{


		px = fmod(*(point_x + i_nd), double(XLENGTH));
		py = fmod(*(point_y + i_nd), double(YLENGTH));
		pz = fmod(*(point_z + i_nd), double(ZLENGTH));
		if (px < 0) { px = px + double(XLENGTH); }
		if (py < 0) { py = py + double(YLENGTH); }
		if (pz < 0) { pz = pz + double(ZLENGTH); }


		fprintf(f3, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", px, py, pz, *(bforceX + i_nd), *(bforceY + i_nd), *(bforceZ + i_nd), *(sforceX + i_nd), *(sforceY + i_nd), *(sforceZ + i_nd));
		fprintf(f3, "\n");

	}


	for (i_nd = 1; i_nd <= tf_main; i_nd++)
	{

		if (*(vtx3 + i_nd) == 1)
		{
			fprintf(f3, "%d\t%d\t%d\n", *(vtx0 + i_nd), *(vtx1 + i_nd), *(vtx2 + i_nd));
			fprintf(f3, "\n");
		}


	}
	*/
	axi = transy = transz = axi2 = transy2 = transz2 = 0.0;

	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{
		for (i_nd2 = i_nd + 1; i_nd2 <= p_num; i_nd2++)
		{
			axi = fabs(*(point_x + i_nd) - *(point_x + i_nd2));
			transy = fabs(*(point_y + i_nd) - *(point_y + i_nd2));
			transz = fabs(*(point_z + i_nd) - *(point_z + i_nd2));

			if (axi > axi2)
			{
				axi2 = axi;
			}
			if (transy >transy2)
			{
				transy2 = transy;
			}
			if (transz >transz2)
			{
				transz2 = transz;
			}


		}

	}
	
	//axi2 = axi2 / Radius/2.0;
	//transy2 = transy2 / Radius/2.0;
	//transz2 = transz2 / Radius/2.0;

	ff7 = fopen("length_x_y.txt", "a");
	fprintf(ff7, "time=%d X= %lf Y= %lf Z= %lf\n", Timek, axi2, transy2, transz2);


	fflush(ff7);
	fclose(ff7);



	fflush(f3);
	fclose(f3);




	return;
}

void    mesh_com2(void)
{
	int     x, y, z, i_nd, pn1, pn2, pn3;
	double  f_x1, f_y1, f_z1, px, py, pz;

	FILE    *f6;
	char FileNamemesh[100];
	char string[100] = { '0' };

	sprintf(FileNamemesh, "c_poinmesh_%d.plt", Timek);


	//  Data for Tecplot

	f6 = fopen(FileNamemesh, "w");
	fprintf(f6, "ZONE");
	fprintf(f6, "\n");



	for (i_nd = 1; i_nd <= tp_num; i_nd++)
	{


		px = fmod(*(point_x + i_nd), double(XLENGTH));
		py = fmod(*(point_y + i_nd), double(YLENGTH));
		pz = fmod(*(point_z + i_nd), double(ZLENGTH));
		if (px < 0) { px = px + double(XLENGTH); }
		if (py < 0) { py = py + double(YLENGTH); }
		if (pz < 0) { pz = pz + double(ZLENGTH); }


		fprintf(f6, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", px, py, pz, *(bforceX + i_nd), *(bforceY + i_nd), *(bforceZ + i_nd), *(sforceX + i_nd), *(sforceY + i_nd), *(sforceZ + i_nd));
		fprintf(f6, "\n");

	}





	fflush(f6);
	fclose(f6);


	return;
}






void timestep_data()
{
	FILE *fP3;

	fP3 = fopen("timestep.txt", "a");
	fprintf(fP3,"%d\n", Timek);
	fprintf(fP3,"\n");
	fflush(fP3);
	fclose(fP3);
	
	return;

}





void	tables	(void)
{
    int x;
    int y;
    int z;
    int i;

	for (i = 1; i <= tp_num; i++)
	{
		 um[i][0] = 0.0;
		um[i][1] = 0.0;
		um[i][2] = 0.0;
	}






    for (x = 1; x <= XLENGTH; x++)
    {
		 for (y = 1; y <= YLENGTH; y++)
      {  for (z = 1; z <= ZLENGTH; z++)
        for (i = 0; i < Q; i++)
	     {   next_x [x][i] = x + cx[i];
		 		if (next_x[x][i] > XLENGTH)
				next_x [x][i] = 1;
	 			if (next_x[x][i] < 1)
				next_x [x][i] = XLENGTH;

             next_y [y][i] = y + cy[i];
		  	  	if (next_y[y][i] > YLENGTH)
				next_y [y][i] = 1;
	 		  	if (next_y[y][i] < 1)
				next_y [y][i] = YLENGTH;
			
             next_z [z][i] = z + cz[i];
		 		if (next_z[z][i] > ZLENGTH)
				next_z [z][i] = 1;
	 			if (next_z[z][i] < 1)
				next_z [z][i] = ZLENGTH;	
 				}
     }
    }




	for (x = 1; x <= XLENGTH; x++)
	{
		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			for (i = 0; i < Q2; i++)
			{
				next_x2[x][i] = x + cx2[i];
				if (next_x2[x][i] > XLENGTH)
					next_x2[x][i] = 1;
				if (next_x2[x][i] < 1)
					next_x2[x][i] = XLENGTH;

				next_y2[y][i] = y + cy2[i];
				if (next_y2[y][i] > YLENGTH)
					next_y2[y][i] = 1;
				if (next_y2[y][i] < 1)
					next_y2[y][i] = YLENGTH;

				next_z2[z][i] = z + cz2[i];
				if (next_z2[z][i] > ZLENGTH)
					next_z2[z][i] = 1;
				if (next_z2[z][i] < 1)
					next_z2[z][i] = ZLENGTH;
			}
		}
	}












    return;
}

/*
void	setparticle (void)
{   int i,j,ii;
    double     PI,pp,radi;
	   FILE *ff2;
      PI = 4.0 * atan(1.0);
	  char mesh_data [100];
	  sprintf(mesh_data, "point.txt");
     
	  if((ff2 = fopen(mesh_data, "r")) == NULL)
		  {
			  return;
		  }		
	     else
		  {
			  for(i = 1; i <= Xm_number; i++)
			  {


					fscanf( ff2, "%lf %lf %lf", &Xm_location[i][0], &Xm_location[i][1], &Xm_location[i][2]);


			  }

			  for(i = 1; i <= Xm_number; i++)
			  {

				  Xm_location[i][0] /= 30.0;
				  Xm_location[i][1] /= 30.0;
				  Xm_location[i][2] /= 30.0;
			  }

			  for(i = 1; i <= Xm_number; i++)
			  {
		          Xm_location[i][0] *= 6.5;
				  Xm_location[i][1] *= 6.5;
				  Xm_location[i][2] *= 6.5;

 			  }



			  for(i = 1; i <= Xm_number; i++)
			  {

		          Xm_location[i][0] += obstX;
				  Xm_location[i][1] += obstY;
				  Xm_location[i][2] += obstZ;


			  }

      fflush( ff2 );
      fclose( ff2 );
		  }





	  for (i = 1; i <= Xm_number; i++)
    {
        RXm_location[i][0] =  Xm_location[i][0];
		RXm_location[i][1] =  Xm_location[i][1];
		RXm_location[i][2] =  Xm_location[i][2];
	}


   return;
}
*/

void	initialise (void)
{   int        x, y, z, i, x0, y0, z0, a, i_nd, i_in, a1, a2, a3, a4, a5;
   


                  for (i = 0; i < Q; i++)
               {
                   v[0][i] = cx[i];
                   v[1][i] = cy[i];                   
                   v[2][i] = cz[i];                   
               }
   

				  for (i = 1; i < rbcn + 1; i++)
				  {

				
						  stiff_b[i] = b_stif;
						  stiff_s[i] = s_stif;
						  stiff_al[i] = al_stif;
						  stiff_at[i] = at_stif;
						  stiff_v[i] = v_stif;
					




				  }
				  
#pragma omp parallel
		{

			/*
				for (x = 1; x <= XLENGTH; x++)
				for (y = 1; y <= YLENGTH; y++)
				for (z = 1; z <= ZLENGTH; z++)
				{
				i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
				*(tau + i_nd) = 1.0;
				}
				*/

#pragma omp	for private(i_nd, y, z)
			for (x = 1; x <= XLENGTH; x++)
				for (y = 1; y <= YLENGTH; y++)
					for (z = 1; z <= ZLENGTH; z++)
					{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				*(bnode + i_nd) = 0;
					}

#pragma omp	for private(i_nd, y, z,i,i_in)
			for (x = 1; x <= XLENGTH; x++)
				for (y = 1; y <= YLENGTH; y++)
					for (z = 1; z <= ZLENGTH; z++)
					{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				//*(bnode + i_nd) = 0.0;

				for (i = 0; i < Q; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					*(b + i_in) = RHO_0 * t[i];
				}

					}
		}

				  /*
				  for (x = 1; x <= XLENGTH; x++)
					  for (y = 1; y <= YLENGTH; y++)
						  for (z = 1; z <= 1; z++)
						  {
					  i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
					  *(bnode + i_nd) = 1;
						  }
				  for (x = 1; x <= XLENGTH; x++)
					  for (y = 1; y <= YLENGTH; y++)
						  for (z = ZLENGTH; z <= ZLENGTH; z++)
						  {
					  i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
					  *(bnode + i_nd) = 1;
						  }

						  */


    return;
}
void	calc_obs ( void )
{   int         x, y, z, i, i_nd, i_in;
    double      sum1, sum2, sum3, sum4;
	double vis_sum, vis_top1, vis_bot1, vis_top2, vis_bot2, vis_top3, vis_bot3;

	vis_top1 = vis_bot1 = vis_top2 = vis_bot2 = vis_top3 = vis_bot3 = 0.0;
#pragma omp parallel
		{
#pragma omp	for private(i_nd, y, z, i, i_in, sum1, sum2, sum3, sum4)
    for (x = 1; x <= XLENGTH; x++)
      for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= ZLENGTH; z++)
      {

		     i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
*(rho + i_nd) =  *(ux + i_nd) = *(uy + i_nd) = *(uz + i_nd) = 0.0;
    
if (*(bnode + i_nd) == 0)
         {
			     sum1 = 0.0;
				 sum2 = 0.0;
				 sum3 = 0.0;
				 sum4 = 0.0;
				 

 	  for (i = 0; i < Q; i++)
          {
			i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
	        sum1	 += *(b + i_in);
            sum2	 += *(b + i_in)*double(cx[i]);
            sum3	 += *(b + i_in)*double(cy[i]);
            sum4	 += *(b + i_in)*double(cz[i]);
	      }
	  *(rho + i_nd) = sum1;
      *(ux + i_nd) = sum2+0.5*(*(FforceX + i_nd) + forcepowerAC);
	  *(uy + i_nd) = sum3+0.5**(FforceY + i_nd);
	  *(uz + i_nd) = sum4+0.5**(FforceZ + i_nd);
         }

         if ( *(rho + i_nd) !=0 )
          {
            *(rhoN + i_nd) = - *(rho + i_nd) / *(rho + i_nd);
	        *(ux + i_nd) /= *(rho + i_nd);
 	        *(uy + i_nd) /= *(rho + i_nd);
 	        *(uz + i_nd) /= *(rho + i_nd);
          }
	}
	


		}
	/*
	if (Timek >= 1) if (Timek % pptime == 0)
	{
		FILE *ff4;
		FILE *ff5;
		FILE *ff6;


	
		for (x = 1; x <= XLENGTH; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{

				vis_top1 += fabs(*(ux + ZLENGTH - 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + ZLENGTH - 1 - 2 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
				vis_bot1 += fabs(*(ux + 3 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + 2 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));

				vis_top2 += fabs(*(ux + ZLENGTH - 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + ZLENGTH - 1 - 3 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
				vis_bot2 += fabs(*(ux + 4 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + 2 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));

				vis_top3 += fabs(*(ux + ZLENGTH - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + ZLENGTH - 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
				vis_bot3 += fabs(*(ux + 2 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));

				                 
			}
		}

		vis_top2 /= 2.0;
		vis_bot2 /= 2.0;

		//vis_top1 *= (ZLENGTH / 4.0 / uMax / XLENGTH / YLENGTH);
		//vis_top2 *= (ZLENGTH / 4.0 / uMax / XLENGTH / YLENGTH);
		//vis_bot1 *= (ZLENGTH / 4.0 / uMax / XLENGTH / YLENGTH);
		//vis_bot2 *= (ZLENGTH / 4.0 / uMax / XLENGTH / YLENGTH);
		
		vis_top3 *= (ZLENGTH / 4.0 / uMax / XLENGTH / YLENGTH);
		vis_bot3 *= (ZLENGTH / 4.0 / uMax / XLENGTH / YLENGTH);

		ff4 = fopen("viscosity1.txt", "a");
		fprintf(ff4, "time=%d top= %lf bot= %lf sum= %lf \n", Timek, vis_top1, vis_bot1, (vis_top1 + vis_bot1) / 2.0);



		ff5 = fopen("viscosity2.txt", "a");
		fprintf(ff5, "time=%d top= %lf bot= %lf sum= %lf \n", Timek, vis_top2, vis_bot2, (vis_top2 + vis_bot2) / 2.0);

		ff6 = fopen("viscosity3.txt", "a");
		fprintf(ff6, "time=%d top= %lf bot= %lf sum= %lf \n", Timek, vis_top3, vis_bot3, (vis_top3 + vis_bot3) / 2.0);






		fflush(ff4);
		fclose(ff4);

		fflush(ff5);
		fclose(ff5);

		fflush(ff6);
		fclose(ff6);


		
	}

	*/





    return;
}
void	project ( void )
{   int         x, y, z, i, i_nd, i_in2, i_in, as, bs;
    double	    udotc, f0, u2, u_mag1, u_mag2, du,LF;

#pragma omp parallel
  {
#pragma omp	for private(y, z, i_in, i_nd, u2, udotc, f0, u_mag1, u_mag2, du,i)
	for (x = 1; x <= XLENGTH; x++)
	  for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= ZLENGTH; z++)
           {   
             i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
     if(*(bnode + i_nd ) != 15)
     	 {
             *(press + i_nd) = *(rho + i_nd)/3.0;
             u2    = (*(ux + i_nd)) * (*(ux + i_nd)) + (*(uy + i_nd)) * (*(uy + i_nd)) + (*(uz + i_nd)) * (*(uz + i_nd));

                if (*(bnode + i_nd) == 0)
                   {	
                  	  for (i = 0; i < Q; i++)
                       {
                         udotc  = (*(ux + i_nd)) * double(cx[i]) + (*(uy + i_nd)) * double(cy[i]) + (*(uz + i_nd)) * double(cz[i]);
			       	     i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
                         *(beq + i_in) =  t[i] * *(rho + i_nd) * (1.0 + 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2);
                        }
                    }
            }        
		}
/*
	#pragma omp	for private(y, z, as, bs, i_nd, i_in)
	for (x = 0; x < XLENGTH; x++)
	  for (y = 0; y < YLENGTH; y++)
        for (z = 0; z < ZLENGTH; z++)
           {   
             i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(bnode + i_nd ) != 15)
		  {
               *(summ + i_nd) =0.0;
          for (as = 0; as < 3; as++)
            for (bs = 0; bs < 3; bs++)
              {  
     			i_in = 3*(as+3*(z+ZLENGTH*(y+YLENGTH*x)))+bs+1;          			
                *(strate + i_in) = 0.0;
              }
          }
		}

       
        
#pragma omp	for private(y, z, i, as, bs, i_in, i_in2, i_nd)
	for (x = 0; x < XLENGTH; x++)
	  for (y = 0; y < YLENGTH; y++)
        for (z = 0; z < ZLENGTH; z++)
           {   
             i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
			      if(*(bnode + i_nd ) != 15)
				 {
          for (as = 0; as < 3; as++)
            for (bs = 0; bs < 3; bs++)
              {  
     			i_in2 = 3*(as+3*(z+ZLENGTH*(y+YLENGTH*x)))+bs+1;          			
              for (i = 0; i < Q; i++)
                   { 
                   i_in = Q*(z+ZLENGTH*(y+YLENGTH*x))+i+1;
                       *(strate + i_in2) = *(strate + i_in2) - 3.0/2.0/ *(tau + i_nd)*(*(b + i_in) - *(beq + i_in))*v[i][as]*v[i][bs];   
                    }
                   *(summ + i_nd) = *(summ + i_nd) +(*(strate + i_in2) * *(strate + i_in2)); 
              }
              *(mu + i_nd) = muinf + (mu0 - muinf) * pow((1.0+  (ramda * ramda * 4.0 * *(summ + i_nd))),(seng-1.0)/2.0);
              *(tau + i_nd) = 3.0 * *(mu + i_nd) / *(rho + i_nd) + 0.5;

            }
		}


*/

				

    #pragma omp	for private(y, z, i_in, i_nd, i,LF)
	for (x = 1; x <= XLENGTH; x++)
	  for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= ZLENGTH; z++)
        {   i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;                         
          for (i = 0; i < Q; i++)
               {
                   if (*(bnode + i_nd) == 0)
                      {

                      		i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
							LF = (1.0 - 0.5*tau)*t[i] * ((3.0 * (double(cx[i]) - *(ux + i_nd)) + 9.0 * (double(cx[i]) * *(ux + i_nd) + double(cy[i]) * *(uy + i_nd) + double(cz[i]) * *(uz + i_nd)) * double(cx[i])) * ((*(FforceX + i_nd) + forcepowerAC)) + (3.0 * (double(cy[i]) - *(uy + i_nd)) + 9.0 * (double(cx[i]) * *(ux + i_nd) + double(cy[i]) * *(uy + i_nd) + double(cz[i]) * *(uz + i_nd)) * double(cy[i])) * *(FforceY + i_nd) + (3.0 * (double(cz[i]) - *(uz + i_nd)) + 9.0 * (double(cx[i]) * *(ux + i_nd) + double(cy[i]) * *(uy + i_nd) + double(cz[i]) * *(uz + i_nd)) * double(cz[i])) * *(FforceZ + i_nd));
							//LF = ((*(FforceX + i_nd) + forcepowerAC) * cx[i] + *(FforceY + i_nd) * cy[i] + *(FforceZ + i_nd) * cz[i]) * 1.5 * t[i];
							*(b + i_in) = (1.0-1.0 / tau) * *(b + i_in) + 1.0 / tau * *(beq + i_in) + LF;
							
                      }
                }
        } 

	/*
	여기하면 벽뚫는다


   #pragma omp	for private(y, z, i_in, i_nd, i)
	for (x = 1; x <= XLENGTH; x++)
	  for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= 1; z++)
        {   i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;                         
          for (i = 0; i < Q; i++)
               {
                   if (*(bnode + i_nd) == 0)
                      {
                  	  i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
                      *(b + i_in) = t[i]*(1.0+3.0*cx[i]*uMax);
							
                      }
                }
        } 
   #pragma omp	for private(y, z, i_in, i_nd, i)
	for (x = 1; x <= XLENGTH; x++)
	  for (y = 1; y <= YLENGTH; y++)
        for (z = ZLENGTH; z <= ZLENGTH; z++)
        {   i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;                         
          for (i = 0; i < Q; i++)
               {
                   if (*(bnode + i_nd) == 0)
                      {
                  	  i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
                      *(b + i_in) = t[i]*(1.0-3.0*cx[i]*uMax);
							
                      }
                }
        } 

  
  

    
    
    */
    	


   }

  
     return;

}

void	propagate(void)
{
	int	x, y, z, i, j, new_x, new_y, new_z, l, h, m, i_in, i_in2, i_in3, i_in4, i_nd;
	double aux, aux1, aux2, rho0, ru, u_Poi, out_pressure, yDiff, zDiff, rr, rrr, xxx, yyy, zzz;
	double out_pressure2;
	int newstep, newstep2, new_x_2, new_y_2, new_z_2, oi, i_2, i_in_2;
	double vec_l_1, vec_l_2, vec_l_0;




	out_pressure = 1.1;

#pragma omp parallel
	{
#pragma omp	for private(y, z, i, new_x, new_y, new_z,i_in,i_in2,i_in3)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		for (i = 0; i < Q; i++)
		{
			i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;

			i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			if (*(bnode + i_in2) != 15)
			{
				new_x = next_x[x][i];
				new_y = next_y[y][i];
				new_z = next_z[z][i];
				i_in3 = Q*(new_z - 1 + ZLENGTH*(new_y - 1 + YLENGTH*(new_x - 1))) + i + 1;
				*(back + i_in3) = *(b + i_in);

			}
		}

#pragma omp	for private(y, z, i_in,i_nd,i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			if (*(bnode + i_nd) != 15)
			{

				for (i = 0; i < Q; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					*(b + i_in) = *(back + i_in);
				}
			}
		}

#pragma omp	for private(y, z, i, i_nd, i_in, new_x_2, new_y_2, new_z_2, vec_l_1, vec_l_2, vec_l_0, oi, i_2, new_x, new_y, new_z, i_in_2)   
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

			if (slip_node[i_nd] == 1)
			{
				for (i = 0; i < Q; i++)
				{


					new_x_2 = cx[i];
					new_y_2 = cy[i];
					new_z_2 = cz[i];
					vec_l_1 = float((cx[i] + SN_x_2[i_nd]) * (cx[i] + SN_x_2[i_nd]) + (cy[i] + SN_y_2[i_nd]) * (cy[i] + SN_y_2[i_nd]) + (cz[i] + SN_z_2[i_nd]) * (cz[i] + SN_z_2[i_nd]));
					vec_l_2 = float((cx[i] - SN_x_2[i_nd]) * (cx[i] - SN_x_2[i_nd]) + (cy[i] - SN_y_2[i_nd]) * (cy[i] - SN_y_2[i_nd]) + (cz[i] - SN_z_2[i_nd]) * (cz[i] - SN_z_2[i_nd]));
					vec_l_0 = float((cx[i]) * (cx[i]) + (cy[i]) * (cy[i]) + (cz[i]) * (cz[i]));

					if ((vec_l_1 > (vec_l_0 + epsilon)) && (vec_l_2 < (vec_l_0 - epsilon)))
					{
						new_x_2 = cx[i] - SN_x_2[i_nd] * 2;
						new_y_2 = cy[i] - SN_y_2[i_nd] * 2;
						new_z_2 = cz[i] - SN_z_2[i_nd] * 2;
					}
					if ((vec_l_1 >(vec_l_0 + epsilon)) && ((vec_l_2 > (vec_l_0 - epsilon)) && (vec_l_2 < (vec_l_0 + epsilon))))
					{
						new_x_2 = cx[i] - SN_x_2[i_nd];
						new_y_2 = cy[i] - SN_y_2[i_nd];
						new_z_2 = cz[i] - SN_z_2[i_nd];

					}

					if ((vec_l_2 >(vec_l_0 + epsilon)) && (vec_l_1 < (vec_l_0 - epsilon)))
					{
						new_x_2 = cx[i] + SN_x_2[i_nd] * 2;
						new_y_2 = cy[i] + SN_y_2[i_nd] * 2;
						new_z_2 = cz[i] + SN_z_2[i_nd] * 2;

					}
					if ((vec_l_2 >(vec_l_0 + epsilon)) && ((vec_l_1 > (vec_l_0 - epsilon)) && (vec_l_1 < (vec_l_0 + epsilon))))
					{
						new_x_2 = cx[i] + SN_x_2[i_nd];
						new_y_2 = cy[i] + SN_y_2[i_nd];
						new_z_2 = cz[i] + SN_z_2[i_nd];

					}

					if ((vec_l_2 >(vec_l_0 + epsilon)) && (vec_l_1 > (vec_l_0 + epsilon)))
					{
						new_x_2 = cx[i];
						new_y_2 = cy[i];
						new_z_2 = cz[i];

					}
					if ((vec_l_2 < (vec_l_0 - epsilon)) && (vec_l_1 < (vec_l_0 - epsilon)))
					{
						new_x_2 = cx[i];
						new_y_2 = cy[i];
						new_z_2 = cz[i];

					}

					oi = -1;
					for (i_2 = 0; i_2 < Q; i_2++)
					{




						if ((new_x_2 == cx[i_2]) && (new_y_2 == cy[i_2]) && (new_z_2 == cz[i_2]))
						{
							oi = i_2;

						}

					}

					if (oi == -1)
					{
						oi = opposite[i];

					}


					new_x = next_x[x][i];
					new_y = next_y[y][i];
					new_z = next_z[z][i];

					//   i_3 = new_z+ZLENGTH[0]*(new_y+YLENGTH[0]*(new_x + XLENGTH[0] * i ));  //Q marker



					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + oi + 1;
					i_in_2 = i + Q * (z - 1 + ZLENGTH * (y - 1 + YLENGTH * (x - 1))) + 1; //Q marker
					back[i_in_2] = b[i_in];



				} // for i


				//sync




			}

		}


#pragma omp	for private(y, z, i, i_nd, i_in, oi, i_in_2)   
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

			if (slip_node[i_nd] == 1)
			{
				for (i = 0; i < Q; i++)
				{


					oi = opposite[i];

					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + oi + 1;
					i_in_2 = i + Q * (z - 1 + ZLENGTH * (y - 1 + YLENGTH * (x - 1))) + 1; //Q marker
					back_2[i_in_2] = b[i_in];



				} // for i


				//sync




			}

		}


#pragma omp	for private(y, z, i, i_nd, i_in)   
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

			if (slip_node[i_nd] == 1)
			{
				for (i = 0; i < Q; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					b[i_in] = back_2[i_in] * (1.0 - slip_ratio) + back[i_in] * slip_ratio;
				} // for i


				//sync
			}
		}

		//wall speed code
		/*
		#pragma omp	for private(y, i_in, i_nd, i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		{
		i_nd = 1 + ZLENGTH*(y + YLENGTH*x) + 1;
		for (i = 0; i < Q; i++)
		{
		if (*(bnode + i_nd) == 0)
		{

		i_in = Q*(1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
		*(b + i_in) = t[i] * (1.0 + 3.0*cx[i] * uOsci);
		}
		}
		}

		#pragma omp	for private(y, z, i_in, i_nd, i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		{
		z = (ZLENGTH - 1);
		i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
		for (i = 0; i < Q; i++)
		{
		if (*(bnode + i_nd) == 0)
		{
		i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
		*(b + i_in) = t[i] * (1.0 - 3.0*cx[i] * uOsci);
		}
		}
		}
		*/


		//pulse_flow inlet 조건
		/*
		#pragma omp	for private(z, x, i, i_in, i_in2,rho0,yDiff,zDiff)          
		for (y = 1; y <= YLENGTH; y++) for (z = 1; z <= ZLENGTH; z++) for (x = 1; x <= 1; x++)
		{
			i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

			if (*(bnode + i_in2) == 0) for (i = 0; i < 1; i++)
			{
				i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;


				//rrr = (double(z) - (double(ZLENGTH) + 1.0) / 2.0)*(double(z) - (double(ZLENGTH) + 1.0) / 2.0) + (double(y) - (double(YLENGTH) + 1.0) / 2.0)*(double(y) - (double(YLENGTH) + 1.0) / 2.0);
				//rr = pow(rrr, 0.5);
				//*(ux + i_in2) = uMax*pulse[pulse_step]*(1.0 + (rr / pipe_r)* (rr / pipe_r));
				*(ux + i_in2) = u_pulse;//*pulse[pulse_step];

				rho0 = (*(b + i_in) + *(b + i_in + 16) + *(b + i_in + 7) + *(b + i_in + 9) + *(b + i_in + 18) + *(b + i_in + 17) + *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 6) + 2.0 *(*(b + i_in + 3) + *(b + i_in + 1) + *(b + i_in + 5) + *(b + i_in + 4) + *(b + i_in + 2))) / (1.0 + *(ux + i_in2));
				yDiff = *(b + i_in + 15) + *(b + i_in + 17) + *(b + i_in + 16) - (*(b + i_in + 6) + *(b + i_in + 7) + *(b + i_in + 8));
				zDiff = *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 18) - (*(b + i_in + 6) + *(b + i_in + 9) + *(b + i_in + 17));

				*(b + i_in + 12) = *(b + i_in + 3) + (1.0 / 3.0)**(ux + i_in2)*rho0;
				*(b + i_in + 10) = *(b + i_in + 1) + (1.0 / 6.0)*rho0**(ux + i_in2) - (1.0 / 2.0) * (yDiff);
				*(b + i_in + 14) = *(b + i_in + 5) + (1.0 / 6.0)*rho0**(ux + i_in2) + (1.0 / 2.0) * (yDiff);
				*(b + i_in + 13) = *(b + i_in + 4) + (1.0 / 6.0)*rho0**(ux + i_in2) + (1.0 / 2.0) * (zDiff);
				*(b + i_in + 11) = *(b + i_in + 2) + (1.0 / 6.0)*rho0**(ux + i_in2) - (1.0 / 2.0) * (zDiff);
			}
					
		}
		*/
		/*
		//free outlet
		#pragma omp	for private(z, x, i, i_in, i_in2, i_in3)          
		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			{
				for (x = XLENGTH; x <= XLENGTH; x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
					
					if (*(bnode + i_in2) == 0)
					{

						//Normal Outlet
						
						{
							i=1;
							*(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+1+1)= 2* *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-2)))+1+1)- *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-3)))+1+1) ;
							*(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+2+1)= 2* *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-2)))+2+1)- *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-3)))+2+1) ;
							*(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+3+1)= 2* *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-2)))+3+1)- *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-3)))+3+1) ;
							*(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+4+1)= 2* *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-2)))+4+1)- *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-3)))+4+1) ;
							*(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+5+1)= 2* *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-2)))+5+1)- *(b + Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-3)))+5+1) ;
						}
						
					}
				}
			}
		}
		*/


	} //End of omp


	return;

}	/* End of procedure - propagate	*/



void	wallefcal(void)
{
	int        x, y, z, i, j, i_in2, x1, x2, y1, y2, z1, z2;
	double     a, b, c, a1, b1, c1, PI, px, py, pz, disg;
	PI = 4.0 * atan(1.0);

	for (i = 1; i <= tp_num; i++)
	{
		px = fmod(*(point_x + i), double(XLENGTH));
		py = fmod(*(point_y + i), double(YLENGTH));
		pz = fmod(*(point_z + i), double(ZLENGTH));
		if (px < 0) { px = px + double(XLENGTH); }
		if (py < 0) { py = py + double(YLENGTH); }
		if (pz < 0) { pz = pz + double(ZLENGTH); }

		//Z 로는 피리오딕 아니다.
		x1 = int(floor(abs(px))) - 1;
		x2 = int(floor(abs(px))) + 2;
		y1 = int(floor(abs(py))) - 1;
		y2 = int(floor(abs(py))) + 2;
		z1 = int(floor(abs(pz))) - 1;
		z2 = int(floor(abs(pz))) + 2;

		if (z1 < 1) { z1 = 1; }
		if (z2 > ZLENGTH) { z2 = ZLENGTH; }



		if ((x1 >= 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))  // 정 정
		{
			for (x = x1; x <= x2; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;




						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;



						if (*(bnode + i_in2) == 1)
						{

							a1 = (a*a + b*b + c*c);
							a1 = sqrt(a1);

							if (a1 <walldis)
							{
								b1 = wallef * (1.0 - a1*a1 / 4.0);


								if (fabs(a) < walldis)
								{
									*(sforceX + i) -= (a / a1 *b1);
								}


								if (fabs(b) < walldis)
								{
									*(sforceY + i) -= (b / a1 *b1);
								}


								if (fabs(c) < walldis)
								{
									*(sforceZ + i) -= (c / a1 *b1);
								}
							}


						}
					
					
					
					
					
					
					}
				}
			}
		}






		if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))     // 좌 정
		{
			for (x = 1; x <= x2; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = x - double(px);
						b = y - double(py);
						c = z - double(pz);


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }






						if (*(bnode + i_in2) == 1)
						{

							a1 = (a*a + b*b + c*c);
							a1 = sqrt(a1);

							if (a1 <walldis)
							{
								b1 = wallef * (1.0 - a1*a1 / 4.0);


								if (fabs(a) < walldis)
								{
									*(sforceX + i) -= (a / a1 *b1);
								}


								if (fabs(b) < walldis)
								{
									*(sforceY + i) -= (b / a1 *b1);
								}


								if (fabs(c) < walldis)
								{
									*(sforceZ + i) -= (c / a1 *b1);
								}
							}


						}









					}
				}
			}

			for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = x - double(px);
						b = y - double(py);
						c = z - double(pz);


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }




						if (*(bnode + i_in2) == 1)
						{

							a1 = (a*a + b*b + c*c);
							a1 = sqrt(a1);

							if (a1 <walldis)
							{
								b1 = wallef * (1.0 - a1*a1 / 4.0);


								if (fabs(a) < walldis)
								{
									*(sforceX + i) -= (a / a1 *b1);
								}


								if (fabs(b) < walldis)
								{
									*(sforceY + i) -= (b / a1 *b1);
								}


								if (fabs(c) < walldis)
								{
									*(sforceZ + i) -= (c / a1 *b1);
								}
							}


						}


					}
				}
			}

		}










		if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))   // 우 정
		{
			for (x = x1; x <= XLENGTH; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
						a = x - double(px);
						b = y - double(py);
						c = z - double(pz);


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }



						if (*(bnode + i_in2) == 1)
						{

							a1 = (a*a + b*b + c*c);
							a1 = sqrt(a1);

							if (a1 <walldis)
							{
								b1 = wallef * (1.0 - a1*a1 / 4.0);


								if (fabs(a) < walldis)
								{
									*(sforceX + i) -= (a / a1 *b1);
								}


								if (fabs(b) < walldis)
								{
									*(sforceY + i) -= (b / a1 *b1);
								}


								if (fabs(c) < walldis)
								{
									*(sforceZ + i) -= (c / a1 *b1);
								}
							}


						}




					}
				}
			}

			for (x = 1; x <= (x2 - XLENGTH); x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = x - double(px);
						b = y - double(py);
						c = z - double(pz);


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }



						if (*(bnode + i_in2) == 1)
						{

							a1 = (a*a + b*b + c*c);
							a1 = sqrt(a1);

							if (a1 <walldis)
							{
								b1 = wallef * (1.0 - a1*a1 / 4.0);


								if (fabs(a) < walldis)
								{
									*(sforceX + i) -= (a / a1 *b1);
								}


								if (fabs(b) < walldis)
								{
									*(sforceY + i) -= (b / a1 *b1);
								}


								if (fabs(c) < walldis)
								{
									*(sforceZ + i) -= (c / a1 *b1);
								}
							}


						}




					}
				}
			}

		}


















	}







	return;




}


void	initialise_2(void)
{
	int        x, y, z, i, i_nd, new_x, new_y, new_z, i_nd2, x2, y2, z2, x4, y4, z4, i_nd3;
	double x3, y3, z3;


#pragma omp parallel
	{

#pragma omp	for private(i_nd, y, z)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			*(slip_node + i_nd) = 0;
		}


#pragma omp	for private(i_nd, y, z, i, new_x, new_y, new_z, i_nd2)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

			if (*(bnode + i_nd) == 0)
			{
				for (i = 0; i<Q2; i++)
				{

					new_x = next_x2[x][i];
					new_y = next_y2[y][i];
					new_z = next_z2[z][i];


					i_nd2 = new_z - 1 + ZLENGTH*(new_y - 1 + YLENGTH*(new_x - 1)) + 1;



					if (*(bnode + i_nd2) == 1)
					{


						*(slip_node + i_nd2) = 1;

					}


				}



			}



		}




#pragma omp	for private(i_nd, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4, i, new_x, new_y, new_z, i_nd3)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			*(SN_x_2 + i_nd) = 0;
			*(SN_y_2 + i_nd) = 0;
			*(SN_z_2 + i_nd) = 0;


			if (slip_node[i_nd] == 1)//if 1
			{









				x2 = y2 = z2 = 0;

				x3 = y3 = z3 = 0.0;
				x4 = y4 = z4 = 0;
				for (i = 1; i<Q2; i++)
				{
					new_x = next_x2[x][i];
					new_y = next_y2[y][i];
					new_z = next_z2[z][i];






					//periodic?

					i_nd3 = new_z - 1 + ZLENGTH*(new_y - 1 + YLENGTH*(new_x - 1)) + 1;
					if (bnode[i_nd3] != 0)
					{

						if ((new_x == 1) && (x == XLENGTH))
						{

							new_x = XLENGTH + 1;


						}
						if ((new_x == XLENGTH) && (x == 1))
						{

							new_x = 0;


						}
						if ((new_y == 1) && (y == YLENGTH))
						{

							new_y = YLENGTH + 1;


						}
						if ((new_y == YLENGTH) && (y == 1))
						{

							new_y = 0;


						}
						if ((new_z == 1) && (z == ZLENGTH))
						{

							new_z = ZLENGTH + 1;


						}
						if ((new_z == ZLENGTH) && (z == 1))
						{

							new_z = 0;


						}

						x2 += x - new_x;
						y2 += y - new_y;
						z2 += z - new_z;

					}

				}
				x3 = float(x2) / pow((float(x2) * float(x2) + float(y2) * float(y2) + float(z2) * float(z2)), float(0.5));
				y3 = float(y2) / pow((float(x2) * float(x2) + float(y2) * float(y2) + float(z2) * float(z2)), float(0.5));
				z3 = float(z2) / pow((float(x2) * float(x2) + float(y2) * float(y2) + float(z2) * float(z2)), float(0.5));


				//*(kin_3[i_nd] = x3;
				//*(kin_4[i_nd] = y3;
				//*(kin_5[i_nd] = z3;



				if ((fabs(fabs(x3) - fabs(y3)) > 0.1) && (fabs(fabs(y3) - fabs(z3)) > 0.1) && (fabs(fabs(z3) - fabs(x3)) > 0.1))
				{

					x2 = y2 = z2 = 0;

					x3 = y3 = z3 = 0.0;

					for (i = 1; i<7; i++)
					{
						new_x = x + cx2[i];
						new_y = y + cy2[i];
						new_z = z + cz2[i];
						if (new_x > XLENGTH)
						{
							new_x = 1;     //sym 
						}
						if (new_x < 1)
						{
							new_x = XLENGTH;
						}
						if (new_y > YLENGTH)
						{
							new_y = 1;     //sym 
						}
						if (new_y < 1)
						{
							new_y = YLENGTH;
						}
						if (new_z > ZLENGTH)
						{
							new_z = 1;     //sym 
						}
						if (new_z < 1)
						{
							new_z = ZLENGTH;
						}






						i_nd3 = new_z - 1 + ZLENGTH*(new_y - 1 + YLENGTH*(new_x - 1)) + 1;
						if (bnode[i_nd3] != 0)
						{

							if ((new_x == 1) && (x == XLENGTH))
							{

								new_x = XLENGTH + 1;


							}
							if ((new_x == XLENGTH) && (x == 1))
							{

								new_x = 0;


							}
							if ((new_y == 1) && (y == YLENGTH))
							{

								new_y = YLENGTH + 1;


							}
							if ((new_y == YLENGTH) && (y == 1))
							{

								new_y = 0;


							}
							if ((new_z == 1) && (z == ZLENGTH))
							{

								new_z = ZLENGTH + 1;


							}
							if ((new_z == ZLENGTH) && (z == 1))
							{

								new_z = 0;


							}



							x2 += x - new_x;
							y2 += y - new_y;
							z2 += z - new_z;

						}

					}
					x3 = float(x2) / pow((float(x2) * float(x2) + float(y2) * float(y2) + float(z2) * float(z2)), float(0.5));
					y3 = float(y2) / pow((float(x2) * float(x2) + float(y2) * float(y2) + float(z2) * float(z2)), float(0.5));
					z3 = float(z2) / pow((float(x2) * float(x2) + float(y2) * float(y2) + float(z2) * float(z2)), float(0.5));



				}

				if (x3 > epsilon2)
				{
					x4 = 1;
					SN_x_2[i_nd] = x4;
					//	   x5 = 1.0;
				}
				else if (x3 < -epsilon2)
				{
					x4 = -1;
					SN_x_2[i_nd] = x4;


					//		   x5 = -1.0;
				}
				else if ((x3 <= epsilon2) && (x3 >= -epsilon2))
				{
					x4 = 0;
					SN_x_2[i_nd] = x4;


					//		   x5 = 0.0;
				}
				if (y3 > epsilon2)
				{
					y4 = 1;
					SN_y_2[i_nd] = y4;
					//		   y5 = 1.0;
				}
				else if (y3 < -epsilon2)
				{
					y4 = -1;
					SN_y_2[i_nd] = y4;
					//		   y5 = -1.0;
				}
				else if ((y3 <= epsilon2) && (y3 >= -epsilon2))
				{
					y4 = 0;
					SN_y_2[i_nd] = y4;
					//			   y5 = 0.0;
				}
				if (z3 > epsilon2)
				{
					z4 = 1;
					SN_z_2[i_nd] = z4;
					//				   z5 = 1.0;
				}
				else if (z3 < -epsilon2)
				{
					z4 = -1;
					SN_z_2[i_nd] = z4;
					//					   z5 = -1.0;
				}
				else if ((z3 <= epsilon2) && (z3 >= -epsilon2))
				{
					z4 = 0;
					SN_z_2[i_nd] = z4;
					//							   z5 = 0.0;
				}


			}//if kin_2






		}





	}

}




void	wallefcal_shear(void)
{


	int        x, y, z, i, j, i_in2, x1, x2, y1, y2, z1, z2;
	double     a, b, c, a1, b1, c1, PI, px, py, pz, disg;

	PI = 4.0 * atan(1.0);
#pragma omp parallel
	{


#pragma omp	for private(x1,x2,y1,y2,z1,z2,x,y,z,i_in2,px,py,pz,a,b,c,a1,b1,c1, disg)
		for (i = 1; i <= tp_num; i++)
		{
			px = fmod(*(point_x + i), XLENGTH);
			py = fmod(*(point_y + i), YLENGTH);
			pz = fmod(*(point_z + i), ZLENGTH);
			if (px < 0) { px = px + XLENGTH; }
			if (py < 0) { py = py + YLENGTH; }
			if (pz < 0) { pz = pz + ZLENGTH; }

			//Z 로는 피리오딕 아니다.
			x1 = int(floor(abs(px))) - 1;
			x2 = int(floor(abs(px))) + 2;
			y1 = int(floor(abs(py))) - 1;
			y2 = int(floor(abs(py))) + 2;
			z1 = int(floor(abs(pz))) - 1;
			z2 = int(floor(abs(pz))) + 2;
			if (z1 < 1) { z1 = 1; }
			if (z2 > ZLENGTH) { z2 = ZLENGTH; }



			if ((x1 >= 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))  // 정 정
			{
				for (x = x1; x <= x2; x++)
				{
					for (y = y1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{
						
							
							

							if ((z == 1) || (z == ZLENGTH))
							{
								
								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}




			}


			if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))     // 좌 정
			{
				for (x = 1; x <= x2; x++)
				{
					for (y = y1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{
						



							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

				for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
				{
					for (y = y1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}
						}
					}
				}

			}



			if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))   // 우 정
			{
				for (x = x1; x <= XLENGTH; x++)
				{
					for (y = y1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

				for (x = 1; x <= (x2 - XLENGTH); x++)
				{
					for (y = y1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

			}




			if ((y1 < 1) && (y2 <= YLENGTH) && (x1 >= 1) && (x2 <= XLENGTH))    // 정 아
			{
				for (y = 1; y <= y2; y++)
				{
					for (x = x1; x <= x2; x++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

				for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
				{
					for (x = x1; x <= x2; x++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

			}



			if ((y1 >= 1) && (y2 > YLENGTH) && (x1 >= 1) && (x2 <= XLENGTH))      // 정 위
			{
				for (y = y1; y <= YLENGTH; y++)
				{
					for (x = x1; x <= x2; x++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

				for (y = 1; y <= (y2 - YLENGTH); y++)
				{
					for (x = x1; x <= x2; x++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

			}



			if ((x1 < 1) && (x2 <= XLENGTH) && (y1 < 1) && (y2 <= YLENGTH))     // 좌 아
			{

				// 남서
				for (x = 1; x <= x2; x++)
				{
					for (y = 1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}


				//  북서
				for (x = 1; x <= x2; x++)
				{
					for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}



				// 북동
				for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
				{
					for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}




				// 남동
				for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
				{
					for (y = 1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

			}




			if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 > YLENGTH))     // 좌 위
			{

				// 남서
				for (x = 1; x <= x2; x++)
				{
					for (y = 1; y <= (y2 - YLENGTH); y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}


				//  북서
				for (x = 1; x <= x2; x++)
				{
					for (y = y1; y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}



				// 북동
				for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
				{
					for (y = y1; y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}




				// 남동
				for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
				{
					for (y = 1; y <= (y2 - YLENGTH); y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

			}



			if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 > YLENGTH))     // 우 위
			{

				// 남서
				for (x = 1; x <= (x2 - XLENGTH); x++)
				{
					for (y = 1; y <= (y2 - YLENGTH); y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}


				//  북서
				for (x = 1; x <= (x2 - XLENGTH); x++)
				{
					for (y = y1; y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}



				// 북동
				for (x = x1; x <= XLENGTH; x++)
				{
					for (y = y1; y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}




				// 남동
				for (x = x1; x <= XLENGTH; x++)
				{
					for (y = 1; y <= (y2 - YLENGTH); y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}

			}






			if ((x1 >= 1) && (x2 > XLENGTH) && (y1 < 1) && (y2 <= YLENGTH))     // 우 아
			{

				// 남서
				for (x = 1; x <= (x2 - XLENGTH); x++)
				{
					for (y = 1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}
						}
					}
				}


				//  북서
				for (x = 1; x <= (x2 - XLENGTH); x++)
				{
					for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}



				// 북동
				for (x = x1; x <= XLENGTH; x++)
				{
					for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}

						}
					}
				}




				// 남동
				for (x = x1; x <= XLENGTH; x++)
				{
					for (y = 1; y <= y2; y++)
					{
						for (z = z1; z <= z2; z++)
						{

							if ((z == 1) || (z == ZLENGTH))
							{

								c = double(z) - pz;
								if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								a1 = (c*c);
								a1 = sqrt(a1);

								if (a1 < 2.0)
								{
									b1 = wallef * (1.0 - a1*a1 / 4.0);


									if (fabs(c) < 2.0)
									{
										*(sforceZ + i) -= (c / a1 *b1);
									}
								}


							}
						}
					}
				}

			}




		}










	}


	return;
}




void	TotalForce(void)
{
	int        x, y, z, i, j, i_in2, x1, x2, y1, y2, z1, z2;
	double     a, b, c, a1, b1, c1, PI, px, py, pz, disg;
	PI = 4.0 * atan(1.0);

#pragma omp parallel
	{

#pragma omp	for private(y, z, i_in2)
		for (x = 1; x <= XLENGTH; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{
				for (z = 1; z <= ZLENGTH; z++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
					*(FforceX + i_in2) = 0.0;
					*(FforceY + i_in2) = 0.0;
					*(FforceZ + i_in2) = 0.0;

				}
			}
		}

	}




	for (i = 1; i <= tp_num; i++)
	{
		px = fmod(*(point_x + i), double(XLENGTH));
		py = fmod(*(point_y + i), double(YLENGTH));
		pz = fmod(*(point_z + i), double(ZLENGTH));
		if (px < 0) { px = px + double(XLENGTH); }
		if (py < 0) { py = py + double(YLENGTH); }
		if (pz < 0) { pz = pz + double(ZLENGTH); }

		//Z 로는 피리오딕 아니다.
		x1 = int(floor(abs(px))) - 1;
		x2 = int(floor(abs(px))) + 2;
		y1 = int(floor(abs(py))) - 1;
		y2 = int(floor(abs(py))) + 2;
		z1 = int(floor(abs(pz))) - 1;
		z2 = int(floor(abs(pz))) + 2;

		if (z1 < 1) { z1 = 1; }
		if (z2 > ZLENGTH) { z2 = ZLENGTH; }



		if ((x1 >= 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))  // 정 정
		{
			for (x = x1; x <= x2; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;


		
						//walleffect 
						/*
						if (*(bnode + i_in2) == 1)
						{
						disg = *(point_x + i)- x;
						if (fabs(disg) < 2.0)
						{
						*(sforceX + i) -= wallef*(1 - disg*disg / 4.0) * disg / fabs(disg) ;
						}



						disg = *(point_y + i) - y;
						if (fabs(disg) < 2.0)
						{

						*(sforceY + i) -= wallef*(1 - disg*disg / 4.0) * disg / fabs(disg);

						}

						disg = *(point_z + i) - z;
						if (fabs(disg) < 2.0)
						{

						*(sforceZ + i) -= wallef*(1 - disg*disg / 4.0)* disg / fabs(disg);
						}



						}
						*/


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						//*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						//*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						//*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i)+ *(rforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i)+ *(rforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i)+ *(rforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}




		}


		if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))     // 좌 정
		{
			for (x = 1; x <= x2; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;



						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

			for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}



		if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))   // 우 정
		{
			for (x = x1; x <= XLENGTH; x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

			for (x = 1; x <= (x2 - XLENGTH); x++)
			{
				for (y = y1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;

						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}




		if ((y1 < 1) && (y2 <= YLENGTH) && (x1 >= 1) && (x2 <= XLENGTH))    // 정 아
		{
			for (y = 1; y <= y2; y++)
			{
				for (x = x1; x <= x2; x++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;

						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

			for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
			{
				for (x = x1; x <= x2; x++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}



		if ((y1 >= 1) && (y2 > YLENGTH) && (x1 >= 1) && (x2 <= XLENGTH))      // 정 위
		{
			for (y = y1; y <= YLENGTH; y++)
			{
				for (x = x1; x <= x2; x++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

			for (y = 1; y <= (y2 - YLENGTH); y++)
			{
				for (x = x1; x <= x2; x++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;

						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}



		if ((x1 < 1) && (x2 <= XLENGTH) && (y1 < 1) && (y2 <= YLENGTH))     // 좌 아
		{

			// 남서
			for (x = 1; x <= x2; x++)
			{
				for (y = 1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}


			//  북서
			for (x = 1; x <= x2; x++)
			{
				for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}



			// 북동
			for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
			{
				for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}




			// 남동
			for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
			{
				for (y = 1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;


						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}




		if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 > YLENGTH))     // 좌 위
		{

			// 남서
			for (x = 1; x <= x2; x++)
			{
				for (y = 1; y <= (y2 - YLENGTH); y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}


			//  북서
			for (x = 1; x <= x2; x++)
			{
				for (y = y1; y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}



			// 북동
			for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
			{
				for (y = y1; y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}




			// 남동
			for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
			{
				for (y = 1; y <= (y2 - YLENGTH); y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}



		if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 > YLENGTH))     // 우 위
		{

			// 남서
			for (x = 1; x <= (x2 - XLENGTH); x++)
			{
				for (y = 1; y <= (y2 - YLENGTH); y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}


			//  북서
			for (x = 1; x <= (x2 - XLENGTH); x++)
			{
				for (y = y1; y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}



			// 북동
			for (x = x1; x <= XLENGTH; x++)
			{
				for (y = y1; y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}




			// 남동
			for (x = x1; x <= XLENGTH; x++)
			{
				for (y = 1; y <= (y2 - YLENGTH); y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);

						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}






		if ((x1 >= 1) && (x2 > XLENGTH) && (y1 < 1) && (y2 <= YLENGTH))     // 우 아
		{

			// 남서
			for (x = 1; x <= (x2 - XLENGTH); x++)
			{
				for (y = 1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}


			//  북서
			for (x = 1; x <= (x2 - XLENGTH); x++)
			{
				for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);


						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}



			// 북동
			for (x = x1; x <= XLENGTH; x++)
			{
				for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}




			// 남동
			for (x = x1; x <= XLENGTH; x++)
			{
				for (y = 1; y <= y2; y++)
				{
					for (z = z1; z <= z2; z++)
					{
						i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

						a = double(x) - px;
						b = double(y) - py;
						c = double(z) - pz;
						if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
						if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
						if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

						if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
						if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
						if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

						a1 = fabs(a);
						b1 = fabs(b);
						c1 = fabs(c);
						*(FforceX + i_in2) += (*(sforceX + i) + *(bforceX + i) + *(vforceX + i) + *(agforceX + i)+ *(aforceX + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceY + i_in2) += (*(sforceY + i) + *(bforceY + i) + *(vforceY + i)+ *(agforceY + i) + *(aforceY + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
						*(FforceZ + i_in2) += (*(sforceZ + i) + *(bforceZ + i) + *(vforceZ + i)+ *(agforceZ + i) + *(aforceZ + i))*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

					}
				}
			}

		}




	}


	//totalforce 계산 여기에


	return;
}



		   


		 
	  









void	UpdatePosition(void)
	   {


		   int        x, y, z, i, j, i_in2, x1, x2, y1, y2, z1, z2;
		   double     a, b, c, a1, b1, c1, PI, px, py, pz;


		   PI = 4.0 * atan(1.0);
#pragma omp parallel
		   {

#pragma omp	for
			   for (i = 1; i <= tp_num; i++)
			   {
				   um[i][0] = 0.0;
				   um[i][1] = 0.0;
				   um[i][2] = 0.0;
			   }




#pragma omp	for private(x1,x2,y1,y2,z1,z2,x,y,z,i_in2,px,py,pz,a,b,c,a1,b1,c1)
			   for (i = 1; i <= tp_num; i++)
			   {
				   px = fmod(*(point_x + i), double(XLENGTH));
				   py = fmod(*(point_y + i), double(YLENGTH));
				   pz = fmod(*(point_z + i), double(ZLENGTH));
				   if (px < 0) { px = px + double(XLENGTH); }
				   if (py < 0) { py = py + double(YLENGTH); }
				   if (pz < 0) { pz = pz + double(ZLENGTH); }

				   //Z 로는 피리오딕 아니다.
				   x1 = int(floor(abs(px))) - 1;
				   x2 = int(floor(abs(px))) + 2;
				   y1 = int(floor(abs(py))) - 1;
				   y2 = int(floor(abs(py))) + 2;
				   z1 = int(floor(abs(pz))) - 1;
				   z2 = int(floor(abs(pz))) + 2;
				   if (z1 < 1) { z1 = 1; }
				   if (z2 > ZLENGTH) { z2 = ZLENGTH; }



				   if ((x1 >= 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))  // 정 정
				   {
					   for (x = x1; x <= x2; x++)
					   {
						   for (y = y1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }




				   }


				   if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))     // 좌 정
				   {
					   for (x = 1; x <= x2; x++)
					   {
						   for (y = y1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

					   for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
					   {
						   for (y = y1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
							   }
						   }
					   }

				   }



				   if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 <= YLENGTH))   // 우 정
				   {
					   for (x = x1; x <= XLENGTH; x++)
					   {
						   for (y = y1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

					   for (x = 1; x <= (x2 - XLENGTH); x++)
					   {
						   for (y = y1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

				   }




				   if ((y1 < 1) && (y2 <= YLENGTH) && (x1 >= 1) && (x2 <= XLENGTH))    // 정 아
				   {
					   for (y = 1; y <= y2; y++)
					   {
						   for (x = x1; x <= x2; x++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

					   for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
					   {
						   for (x = x1; x <= x2; x++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

				   }



				   if ((y1 >= 1) && (y2 > YLENGTH) && (x1 >= 1) && (x2 <= XLENGTH))      // 정 위
				   {
					   for (y = y1; y <= YLENGTH; y++)
					   {
						   for (x = x1; x <= x2; x++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

					   for (y = 1; y <= (y2 - YLENGTH); y++)
					   {
						   for (x = x1; x <= x2; x++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);
								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

				   }



				   if ((x1 < 1) && (x2 <= XLENGTH) && (y1 < 1) && (y2 <= YLENGTH))     // 좌 아
				   {

					   // 남서
					   for (x = 1; x <= x2; x++)
					   {
						   for (y = 1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }


					   //  북서
					   for (x = 1; x <= x2; x++)
					   {
						   for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }



					   // 북동
					   for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
					   {
						   for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }




					   // 남동
					   for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
					   {
						   for (y = 1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

				   }




				   if ((x1 < 1) && (x2 <= XLENGTH) && (y1 >= 1) && (y2 > YLENGTH))     // 좌 위
				   {

					   // 남서
					   for (x = 1; x <= x2; x++)
					   {
						   for (y = 1; y <= (y2 - YLENGTH); y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);
								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }


					   //  북서
					   for (x = 1; x <= x2; x++)
					   {
						   for (y = y1; y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }



					   // 북동
					   for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
					   {
						   for (y = y1; y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }




					   // 남동
					   for (x = (x1 + XLENGTH); x <= XLENGTH; x++)
					   {
						   for (y = 1; y <= (y2 - YLENGTH); y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

				   }



				   if ((x1 >= 1) && (x2 > XLENGTH) && (y1 >= 1) && (y2 > YLENGTH))     // 우 위
				   {

					   // 남서
					   for (x = 1; x <= (x2 - XLENGTH); x++)
					   {
						   for (y = 1; y <= (y2 - YLENGTH); y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }


								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }


					   //  북서
					   for (x = 1; x <= (x2 - XLENGTH); x++)
					   {
						   for (y = y1; y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }



					   // 북동
					   for (x = x1; x <= XLENGTH; x++)
					   {
						   for (y = y1; y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }




					   // 남동
					   for (x = x1; x <= XLENGTH; x++)
					   {
						   for (y = 1; y <= (y2 - YLENGTH); y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }

				   }






				   if ((x1 >= 1) && (x2 > XLENGTH) && (y1 < 1) && (y2 <= YLENGTH))     // 우 아
				   {

					   // 남서
					   for (x = 1; x <= (x2 - XLENGTH); x++)
					   {
						   for (y = 1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);
								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }


					   //  북서
					   for (x = 1; x <= (x2 - XLENGTH); x++)
					   {
						   for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);

								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }



					   // 북동
					   for (x = x1; x <= XLENGTH; x++)
					   {
						   for (y = (y1 + YLENGTH); y <= YLENGTH; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }

								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);


								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;

							   }
						   }
					   }




					   // 남동
					   for (x = x1; x <= XLENGTH; x++)
					   {
						   for (y = 1; y <= y2; y++)
						   {
							   for (z = z1; z <= z2; z++)
							   {
								   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


								   a = x - double(px);
								   b = y - double(py);
								   c = z - double(pz);


								   if (a > (double(XLENGTH) / 2.0)) { a = a - double(XLENGTH); }
								   if (b > (double(YLENGTH) / 2.0)) { b = b - double(YLENGTH); }
								   if (c > (double(ZLENGTH) / 2.0)) { c = c - double(ZLENGTH); }

								   if (a < (-double(XLENGTH) / 2.0)) { a = a + double(XLENGTH); }
								   if (b < (-double(YLENGTH) / 2.0)) { b = b + double(YLENGTH); }
								   if (c < (-double(ZLENGTH) / 2.0)) { c = c + double(ZLENGTH); }
								   a1 = fabs(a);
								   b1 = fabs(b);
								   c1 = fabs(c);



								   um[i][0] += *(ux + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][1] += *(uy + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
								   um[i][2] += *(uz + i_in2)*(1.0 + cos(a1*PI / 2.0))*(1.0 + cos(b1*PI / 2.0))*(1.0 + cos(c1*PI / 2.0)) / 64.0;
							   }
						   }
					   }

				   }




			   }


			   // particle moves by velocity
			   //좌표 바뀌는 지점
#pragma omp	for
			   for (i = 1; i <= tp_num; i++)
			   {
				   *(point_x + i) += um[i][0];
				   *(point_y + i) += um[i][1];
				   *(point_z + i) += um[i][2];
			   }





			   
			   





		   }


		   return;
	   }


void	ParticleContour (void)
{  
	int        x, y, z, i, j, i_in2;
    double     a, b, c, a1, b1, c1, cx, cy, px, py, pz;

    



#pragma omp parallel
	{






#pragma omp for private(y,z,i_in2)	
   for (x = 1; x <= XLENGTH; x++)
      {
         for (y = 1; y <= YLENGTH; y++)
		 {         
			 for (z = 1; z <= ZLENGTH; z++)
          {
             	i_in2 = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1; 
		     	image_particle[i_in2] = 0 ;




		  }
          }
          }






#pragma omp for private(y,z,i_in2,i,px,py,pz,a,b,c,a1,b1,c1)	
   for (x = 1; x <= XLENGTH; x++)
      {
	   for (y = 1; y <= YLENGTH; y++)
	   {
		   for (z = 1; z <= ZLENGTH; z++)
		   {
			   i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			   for (i = 1; i <= tp_num; i++)
			   {

				   px = fmod(*(point_x + i), XLENGTH);
				   py = fmod(*(point_y + i), YLENGTH);
				   pz = fmod(*(point_z + i), ZLENGTH);

				   if (px < 0) { px = px + XLENGTH; }
				   if (py < 0) { py = py + YLENGTH; }
				   if (pz < 0) { pz = pz + ZLENGTH; }

				   a = x - px;
				   b = y - py;
				   c = z - pz;

				   if (a > (XLENGTH / 2.0)) { a = a - XLENGTH; }
				   if (b > (XLENGTH / 2.0)) { b = b - YLENGTH; }
				   if (c > (XLENGTH / 2.0)) { c = c - ZLENGTH; }

				   if (a < (-XLENGTH / 2.0)) { a = a + XLENGTH; }
				   if (b < (-YLENGTH / 2.0)) { b = b + YLENGTH; }
				   if (c < (-ZLENGTH / 2.0)) { c = c + ZLENGTH; }

				   a1 = fabs(a);
				   b1 = fabs(b);
				   c1 = fabs(c);





				   if ((a1 <= 0.6) && (b1 <= 0.6) && (c1 <= 0.6))
				   {

					   image_particle[i_in2] = 3;

				   }




			   }
		   }
	   }
   }

   }


   /*
   
   for (x = 1; x <= XLENGTH; x++)
      {
         for (y = 1; y <= YLENGTH; y++)
          {
            for (z = 1; z <= ZLENGTH; z++)
			{    	i_in2 = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1; 
   			for (i = p_num; i <= p_num; i++)       
             {
				px = fmod(*(point_x + i),XLENGTH);
				py = fmod(*(point_y + i),YLENGTH);
				pz = fmod(*(point_z + i),ZLENGTH);

				if(px < 0) {px = px + XLENGTH ;}
				if(py < 0) {py = py + YLENGTH ;}
				if(pz < 0) {pz = pz + ZLENGTH ;}

	            a = x -  px;
                b = y -  py;
				c = z -  pz;
				
				if(a > (XLENGTH/2.0)) {a = a - XLENGTH; }
				if(b > (XLENGTH/2.0)) {b = b - YLENGTH; }
				if(c > (XLENGTH/2.0)) {c = c - ZLENGTH; }

				if(a < (-XLENGTH/2.0)) {a = a + XLENGTH; }
				if(b < (-YLENGTH/2.0)) {b = b + YLENGTH; }
				if(c < (-ZLENGTH/2.0)) {c = c + ZLENGTH; }

				a1 =  fabs(a);
                b1 =  fabs(b);
                c1 =  fabs(c);     				 
				 
                 
             if (( a1 <= 1.5) && ( b1 <= 1.5) && ( c1 <= 1.5))
              {

                   image_particle[i_in2] = 3 ;

              }
                   

           
           
             }
          }
        }
   }


   */



   



    return;
}





void	springcal(void)
{
	int     i_nd, i_nd2, ii, kk, kk1, i;
	double  c_length, aa, Fmag, vecx, vecy, vecz, bb, cc, disp, xp, yp, zp, shearmo, shearde;


#pragma omp parallel
	{
#pragma omp	for	
		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{
			*(sforceX + i_nd) = 0.0;
			*(sforceY + i_nd) = 0.0;
			*(sforceZ + i_nd) = 0.0;

		}



#pragma omp for private(i_nd,kk,kk1,aa,c_length,Fmag,vecx,vecy,vecz,bb,cc,shearde, shearmo)	
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		{

			for (i_nd = ((i_nd2 - 1)*l_main + 1); i_nd <= l_main*i_nd2; i_nd++)
			{


				kk = *(line_twopoint1 + i_nd);
				kk1 = *(line_twopoint2 + i_nd);
				//aa는 두 점 사이의 거리
				aa = pow(*(point_x + kk1) - *(point_x + kk), 2) + pow(*(point_y + kk1) - *(point_y + kk), 2) + pow(*(point_z + kk1) - *(point_z + kk), 2);
				c_length = pow(aa, 0.5);
				//line_ref는 기준 거리, shearde는 변한 거리의 비율
				shearde = c_length / *(line_ref + i_nd);
				//shearmo는 ??
				shearmo = (pow(shearde, 0.5) + pow(shearde, -2.5)) / (shearde + pow(shearde, -3.0));
				//Fmag = shearmo * (c_length - *(line_ref + i_nd));


				Fmag = shearmo * stiff_s[i_nd2] * (c_length - *(line_ref + i_nd)) ; //실제 force
				vecx = *(point_x + kk) - *(point_x + kk1); //벡터 방향
				vecy = *(point_y + kk) - *(point_y + kk1);
				vecz = *(point_z + kk) - *(point_z + kk1);

				bb = vecx * vecx + vecy * vecy + vecz * vecz; 
				cc = pow(bb, 0.5);//유닛벡터

				vecx /= cc;//유닛방향벡터
				vecy /= cc;
				vecz /= cc;

				*(sforceX + kk) -= Fmag * vecx;
				*(sforceY + kk) -= Fmag * vecy;
				*(sforceZ + kk) -= Fmag * vecz;

				*(sforceX + kk1) += Fmag * vecx;
				*(sforceY + kk1) += Fmag * vecy;
				*(sforceZ + kk1) += Fmag * vecz;







			}
		}
	}
	/*
	for (i = 1; i <= tp_num; i++)
	{

	if (*(rpoint + i) == 1)
	{
	*(sforceX + i) -= opforce;

	}


	if (*(rpoint + i) == 2)
	{
	*(sforceX + i) += opforce;

	}


	}
	
	*/

	/*
	for(i_nd = 1; i_nd <= l_main; i_nd++)
	{
	kk    = *(line_twopoint1 + i_nd);
	kk1   = *(line_twopoint2 + i_nd);
	aa  =   pow(*(point_x + kk1) - *(point_x + kk),2) + pow(*(point_y + kk1) - *(point_y + kk),2) + pow(*(point_z + kk1) - *(point_z + kk),2);
	c_length = pow(aa,0.5);

	SFmag = s_stif * (c_length- *(line_ref + i_nd));

	vecx = *(point_x + kk) - *(point_x + kk1);
	vecy = *(point_y + kk) - *(point_y + kk1);
	vecz = *(point_z + kk) - *(point_z + kk1);

	bb = vecx * vecx + vecy * vecy + vecz * vecz;
	cc = pow(bb, 0.5);

	vecx /= cc;
	vecy /= cc;
	vecz /= cc;

	*(sforceX + kk) -= SFmag * vecx;
	*(sforceY + kk) -= SFmag * vecy;
	*(sforceZ + kk) -= SFmag * vecz;

	*(sforceX + kk1) += SFmag * vecx;
	*(sforceY + kk1) += SFmag * vecy;
	*(sforceZ + kk1) += SFmag * vecz;


	}








	*/

	if (Timek >= 1) if (Timek % ptime == 0)
	{
		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{
			*(Deff + i_nd) = 0.0;
			*(defn + i_nd) = 0.0;
			
		}

	
		for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
		{

			for (i_nd = ((i_nd2 - 1)*l_main + 1); i_nd <= l_main*i_nd2; i_nd++)
			{


				kk = *(line_twopoint1 + i_nd);
				kk1 = *(line_twopoint2 + i_nd);
				aa = pow(*(point_x + kk1) - *(point_x + kk), 2) + pow(*(point_y + kk1) - *(point_y + kk), 2) + pow(*(point_z + kk1) - *(point_z + kk), 2);
				c_length = pow(aa, 0.5);
				
				*(defn + kk) += 1.0;
				*(defn + kk1) += 1.0;

				Fmag = (c_length / *(line_ref + i_nd)-1.0) * (c_length / *(line_ref + i_nd)-1.0);


				*(Deff + kk) += Fmag;
				*(Deff + kk1) += Fmag;
				
				
				







			}
		}


		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{
			Fmag = *(Deff + i_nd) / *(defn + i_nd);
			*(Deff + i_nd) = pow(Fmag, 0.5);

		}



	
	}






	return;
}



void	constructimage(void)
{
	int i, ii, i_nd, i_nd2, new_x, new_y, new_z;
	int x, y, z, x1;
	double c_x1, c_y1, c_z1, c_x2, c_y2, c_z2, c_y3, ycenn, xx;
	double p_x, p_y, p_z, pipe_distance, chradius2, chradius3;

	p_x = XLENGTH/2+1;
	p_y = YLENGTH/2+1;
	p_z = ZLENGTH/2+1;
	chradius2 = chRADIUS * chRADIUS;

	FILE *ff12;


	for (z = 1; z <= ZLENGTH; z++)
	{
		for (x = 1; x <= XLENGTH; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{

				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

				*(bnode + i_nd) = 1;
			}
		}
	}


	for (z = 1; z <= ZLENGTH; z++)
	{
		for (x = 1; x <= XLENGTH; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{

				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				if(z!=1 && z!=ZLENGTH && y!=1 && y!=YLENGTH) *(bnode + i_nd) = 0;

			}
		}
	}

	for (z = 1; z <= ZLENGTH; z++)
	{
		for (x = 1; x <= XLENGTH; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				pipe_distance = (y-p_y)*(y-p_y) + (z-p_z)*(z-p_z);
				if (pipe_distance >= chradius2) *(bnode + i_nd) = 1;

			}
		}
	}


	for (z = 1; z <= ZLENGTH; z++)
	{
		for (x = (XLENGTH-ALENGTH)/2; x <= (XLENGTH+ALENGTH)/2; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				i = x - (XLENGTH-ALENGTH)/2;
				pipe_distance = (y-p_y)*(y-p_y) + (z-p_z)*(z-p_z);
				chradius3 = sin(i*PI2/ALENGTH)*(YLENGTH-chRADIUS-12)/4 + chRADIUS;
				chradius3 = chradius3*chradius3;

				if (pipe_distance < chradius3) *(bnode + i_nd) = 0;

			}
		}
	}




	vvol = 0;

	for (z = 1; z <= ZLENGTH; z++)
	{
		for (x = 1; x <= 1; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{

				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;



				if (*(bnode + i_nd) == 0)
				{
					vvol += 1;
				}

			}
		}
	}


	ff12 = fopen("areavol.txt", "a");
	fprintf(ff12, "%d  \n", vvol);




	fflush(ff12);
	fclose(ff12);



}


void	flowrate(void)
{
	int i, ii, i_nd, i_nd2, x1, i_in2;
	int x, y, z;
	double c_x1, c_y1, c_z1, c_x2, c_y2, c_z2, c_y3, ycenn, xx, yy, zz, uup, udown, sumb, sumb2, sums, sums2, sumv, sumv2, sumdef, suma, suma2, sumga, sumga2, meanvel, relviscosity;


	FILE *ff6;
	FILE *ff8;
	FILE *ff9;

	FILE *ff10;
	FILE *ff11;
	FILE *ff12;



	if (Timek < (aftermoving - 1))
	{
		meanvel_t = 0.0;

		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			{
				for (x = 1; x <= XLENGTH; x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

					if (*(bnode + i_in2) == 0)
					{

						meanvel_t += *(ux + i_in2);

					}
				}
			}
		}


		meanvel_t /= vvol;

	}


	meanvel = 0.0;
	for (z = 1; z <= ZLENGTH; z++)
	{
		for (x = 1; x <= XLENGTH; x++)
		{
			for (y = 1; y <= YLENGTH; y++)
			{

				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;



				if (*(bnode + i_nd) == 0)
				{
					meanvel += *(ux + i_nd);
				}

			}
		}
	}

	meanvel /= double(vvol);
	relviscosity = meanvel_t / meanvel;


	ff12 = fopen("areavol.txt", "a");
	fprintf(ff12, " time = %d meanvel = %lf  \n", Timek, meanvel);




	fflush(ff12);
	fclose(ff12);


	sumb = sumb2 = sums = sums2 = sumv = sumv2 = sumdef = 0.0;
	suma = suma2 = sumga = sumga2 = 0.0;
	if (Timek < (aftermoving - 1))
	{
		flowori = 0.0;

		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			{
				for (x = (XLENGTH - 10); x <= (XLENGTH - 10); x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

					if (*(bnode + i_in2) == 0)
					{

						flowori += *(ux + i_in2);

					}
				}
			}
		}

	}

	uup = 0.0;
	udown = 0.0;
	for (y = 1; y <= YLENGTH; y++)
	{
		for (z = 1; z <= ZLENGTH; z++)
		{
			for (x = (XLENGTH - 10); x <= (XLENGTH - 10); x++)
			{
				i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

				if (*(bnode + i_in2) == 0)
				{

					udown += *(ux + i_in2);

				}
			}
		}
	}


	for (y = (YLENGTH + 1) / 2; y <= YLENGTH; y++)
	{
		for (z = 1; z <= ZLENGTH; z++)
		{
			for (x = 1; x <= 1; x++)
			{
				i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

				if (*(bnode + i_in2) == 0)
				{


					for (i = 0; i < 1; i++)
					{

						uup += *(ux + i_in2);
					}
				}
			}
		}
	}

	xx = yy = zz = 0.0;


	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{
		xx += *(point_x + i_nd);
		yy += *(point_y + i_nd);
		zz += *(point_z + i_nd);
	}

	xx /= double(p_num);
	yy /= double(p_num);
	zz /= double(p_num);
	xx /= double(XLENGTH);
	time_cx = xx;

	ff6 = fopen("flowrate.txt", "a");
	//fprintf(ff6, "time = %d up = %.9lf down = %.9lf dif = %.9lf \n",Timek, udown, udown, udown/(uup+udown));

	fprintf(ff6, "time = %d rel_vis = %lf flowrate = %.9lf viscosity = %.9lf \n", Timek, relviscosity, udown, udown / flowori);

	ff8 = fopen("pointXYZ.txt", "a");


	fprintf(ff8, "time = %lf X = %lf Y = %lf Z = %lf \n", (Timek - aftermoving), time_cx, yy, zz);


	ff9 = fopen("pointX1.txt", "a");



	fprintf(ff9, "%lf \n", xx);




	fflush(ff6);
	fclose(ff6);


	fflush(ff8);
	fclose(ff8);


	fflush(ff9);
	fclose(ff9);







	sumb = sumb2 = sums = sums2 = sumv = sumv2 = sumdef = 0.0;
	suma = suma2 = sumga = sumga2 = 0.0;

	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		sumdef += *(Deff + i_nd); //길이변화율-1의 누적값

	}
	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(sforceX + i_nd);
		yy += *(sforceY + i_nd);
		zz += *(sforceZ + i_nd);

	}
	sums = xx * xx + yy * yy + zz * zz;

	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(sforceX + i_nd) * *(sforceX + i_nd);
		yy += *(sforceY + i_nd) * *(sforceY + i_nd);
		zz += *(sforceZ + i_nd) * *(sforceZ + i_nd);

	}
	sums2 = xx + yy + zz;

	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(aforceX2 + i_nd);
		yy += *(aforceY2 + i_nd);
		zz += *(aforceZ2 + i_nd);

	}
	suma = xx * xx + yy * yy + zz * zz;



	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(aforceX2 + i_nd) * *(aforceX2 + i_nd);
		yy += *(aforceY2 + i_nd) * *(aforceY2 + i_nd);
		zz += *(aforceZ2 + i_nd) * *(aforceZ2 + i_nd);

	}
	suma2 = xx + yy + zz;


	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(aforceX3 + i_nd);
		yy += *(aforceY3 + i_nd);
		zz += *(aforceZ3 + i_nd);

	}
	sumga = xx * xx + yy * yy + zz * zz;



	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(aforceX3 + i_nd) * *(aforceX3 + i_nd);
		yy += *(aforceY3 + i_nd) * *(aforceY3 + i_nd);
		zz += *(aforceZ3 + i_nd) * *(aforceZ3 + i_nd);

	}
	sumga2 = xx + yy + zz;


	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(bforceX + i_nd);
		yy += *(bforceY + i_nd);
		zz += *(bforceZ + i_nd);

	}
	sumb = xx * xx + yy * yy + zz * zz;



	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(bforceX + i_nd) * *(bforceX + i_nd);
		yy += *(bforceY + i_nd) * *(bforceY + i_nd);
		zz += *(bforceZ + i_nd) * *(bforceZ + i_nd);

	}
	sumb2 = xx + yy + zz;




	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(vforceX + i_nd);
		yy += *(vforceY + i_nd);
		zz += *(vforceZ + i_nd);

	}
	sumv = xx * xx + yy * yy + zz * zz;



	xx = yy = zz = 0.0;
	for (i_nd = 1; i_nd <= p_num; i_nd++)
	{

		xx += *(vforceX + i_nd) * *(vforceX + i_nd);
		yy += *(vforceY + i_nd) * *(vforceY + i_nd);
		zz += *(vforceZ + i_nd) * *(vforceZ + i_nd);

	}
	sumv2 = xx + yy + zz;

	sumv = pow(sumv, 0.5)*1000.0;;
	suma = pow(suma, 0.5)*1000.0;;
	sumga = pow(sumga, 0.5)*1000.0;;
	sumb = pow(sumb, 0.5)*1000.0;;
	sums = pow(sums, 0.5)*1000.0;;
	sumdef = pow(sumdef, 0.5)*1000.0;;

	sumv2 = pow(sumv2, 0.5)*1000.0;;
	suma2 = pow(suma2, 0.5)*1000.0;;
	sumga2 = pow(sumga2, 0.5)*1000.0;;
	sumb2 = pow(sumb2, 0.5)*1000.0;;
	sums2 = pow(sums2, 0.5)*1000.0;;



	if (Timek < (aftermoving - 1))
	{

		ff10 = fopen("slip_vel.txt", "a");

		for (i_nd = 1; i_nd <= YLENGTH; i_nd++)
		{


			i_in2 = (ZLENGTH + 1) / 2 - 1 + ZLENGTH*(i_nd - 1 + YLENGTH*((XLENGTH + 1) / 2 - 1)) + 1;
			fprintf(ff10, "%d   %d  %lf \n", Timek, i_nd, *(ux + i_in2));




		}
		fprintf(ff10, "\n");

		fflush(ff10);
		fclose(ff10);
	}


	ff11 = fopen("sum2.txt", "a");
	fprintf(ff11, "%lf  %lf  %lf  %lf  %lf \n", sumv2, suma2, sumga2, sumb2, sums2);




	fflush(ff11);
	fclose(ff11);

}


void	viscosityshear(void)
{
	int i, ii, i_nd, i_nd2, x1, i_in2;
	int x, y, z;
	double c_x1, c_y1, c_z1, c_x2, c_y2, c_z2, c_y3, ycenn, xx, yy, zz, uup, udown;
	double vis_sum, vis_top1, vis_bot1, vis_top2, vis_bot2, vis_top3, vis_bot3;
	FILE *ff4;
	FILE *ff5;
	FILE *ff6;

	vis_top1 = vis_bot1 = vis_top2 = vis_bot2 = vis_top3 = vis_bot3 = 0.0;

	for (x = 1; x <= XLENGTH; x++)
	{
		for (y = 1; y <= YLENGTH; y++)
		{

			vis_top1 += fabs(*(ux + ZLENGTH - 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + ZLENGTH - 1 - 2 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
			vis_bot1 += fabs(*(ux + 3 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + 2 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));

			vis_top2 += fabs(*(ux + ZLENGTH - 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + ZLENGTH - 1 - 3 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
			vis_bot2 += fabs(*(ux + 4 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + 2 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));

			vis_top3 += fabs(*(ux + ZLENGTH - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + ZLENGTH - 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
			vis_bot3 += fabs(*(ux + 2 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1) - *(ux + 1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1));
		}
	}

	vis_top2 /= 2.0;
	vis_bot2 /= 2.0;

	vis_top1 *= (1.0 / XLENGTH / YLENGTH);
	vis_top2 *= (1.0 / XLENGTH / YLENGTH);
	vis_bot1 *= (1.0 / XLENGTH / YLENGTH);
	vis_bot2 *= (1.0 / XLENGTH / YLENGTH);

	vis_top3 *= (1.0 / XLENGTH / YLENGTH);
	vis_bot3 *= (1.0 / XLENGTH / YLENGTH);

	ff4 = fopen("viscosity1.txt", "a");
	fprintf(ff4, "%d\t%lf\t%lf\t%lf\t%lf\n", Timek, vis_top1, vis_bot1, (vis_top1 + vis_bot1) / 2.0, uOsci);

	ff5 = fopen("viscosity2.txt", "a");
	fprintf(ff5, "%d\t%lf\t%lf\t%lf\t%lf\n", Timek, vis_top2, vis_bot2, (vis_top2 + vis_bot2) / 2.0, uOsci);

	ff6 = fopen("viscosity3.txt", "a");
	fprintf(ff6, "%d\t%lf\t%lf\t%lf\t%lf\n", Timek, vis_top3, vis_bot3, (vis_top3 + vis_bot3) / 2.0, uOsci);






	fflush(ff4);
	fclose(ff4);

	fflush(ff5);
	fclose(ff5);

	fflush(ff6);
	fclose(ff6);







}

void	flowprofile(void)
{

		int i, ii, i_nd, i_nd2, x1, i_in2;
		int x, y, z;


	FILE *ff15;


	ff15 = fopen("viscosity3.txt", "a");

    fprintf(ff15, "time=%d \n", Timek);

	for (x = 1; x <=1; x++)
	{
		for (y = ((YLENGTH + 1) / 2 - int(pipe_r)); y <= ((YLENGTH + 1) / 2 + int(pipe_r)); y++)
		{

			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;


			fprintf(ff15, "D = %d U = %lf V = %lf W = %lf \n", y, *(ux + i_nd), *(uy + i_nd), *(uz + i_nd));




		}
	}


	fflush(ff15);
	fclose(ff15);


}


void aggregation(void)
{

	int     i_nd, i_nd2, i_nd3, i_nd4, pn1, pn2, pn3, pn4, pn5, pn6, ap1, ap2, ap3, ap4, ap5, ap6, i,j;
	double   nor_ax, nor_ay, nor_az, nor_bx, nor_by, nor_bz, nor_mag, PI, mx, my, mz;
	double   ax1, ay1, az1, ax2, ay2, az2, abc1, abc2, abc3, abc4, abc5;
	double 	aa, bb, cc, c_length, Fmag, vecx, vecy, vecz, Morse;


	
	double xdis, ydis, zdis, px1, px2, py1, py2;

	PI = 4.0 * atan(1.0);


	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	{



		if (*(point_x + ((i_nd2 - 1)*p_num + 1)) < (1.0 - Radius))
		{
			for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
			{

				*(point_x + i_nd) += XLENGTH;

			}


		}


		if (*(point_x + ((i_nd2 - 1)*p_num + 1)) > (XLENGTH + Radius))
		{
			for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
			{

				*(point_x + i_nd) -= XLENGTH;

			}


		}


		if (*(point_y + ((i_nd2 - 1)*p_num + 1)) < (1.0 - Radius))
		{
			for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
			{

				*(point_y + i_nd) += YLENGTH;

			}


		}


		if (*(point_y + ((i_nd2 - 1)*p_num + 1)) > (YLENGTH + Radius))
		{
			for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
			{

				*(point_y + i_nd) -= YLENGTH;

			}


		}



	}






	for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)
	{
		
		    cen_rx[i_nd2] = 0.0;
			cen_ry[i_nd2] = 0.0;
			cen_rz[i_nd2] = 0.0;

			for (i_nd = ((i_nd2 - 1)*p_num + 1); i_nd <= p_num*i_nd2; i_nd++)
			{

			cen_rx[i_nd2] += *(point_x + i_nd);
			cen_ry[i_nd2] += *(point_y + i_nd);
			cen_rz[i_nd2] += *(point_z + i_nd);
			}
		
		cen_rx[i_nd2] /= p_num;
		cen_ry[i_nd2] /= p_num;
		cen_rz[i_nd2] /= p_num;

	}




	for (i_nd = 1; i_nd <= tp_num; i_nd++)
	{


		*(agforceX + i_nd) = 0.0;
		*(agforceY + i_nd) = 0.0;
		*(agforceZ + i_nd) = 0.0;
	}


	for (i_nd = 1; i_nd <= rbcn - 1; i_nd++)
	{
		for (i_nd2 = i_nd + 1; i_nd2 <= rbcn; i_nd2++)
		{



			xdis = cen_rx[i_nd] - cen_rx[i_nd2];
			ydis = cen_ry[i_nd] - cen_ry[i_nd2];
			zdis = cen_rz[i_nd] - cen_rz[i_nd2];


			if (fabs(xdis) > (XLENGTH - safezone )) //safezone 20으로 설정되어 있음
			{

				xdis = fabs(xdis) - XLENGTH; //safezone보다 거리가 클 경우 dis는 반대쪽으로
			}


			if (fabs(ydis) > (YLENGTH - safezone))
			{

				ydis = fabs(ydis) - YLENGTH;
			}

			


			if (fabs(xdis) < safezone)
			{
				if (fabs(ydis) < safezone)
				{
					if (fabs(zdis) < safezone)
					{

						for (i_nd3 = ((i_nd - 1)*p_num + 1); i_nd3 <= p_num*i_nd; i_nd3++)
						{
							for (i_nd4 = ((i_nd2 - 1)*p_num + 1); i_nd4 <= p_num*i_nd2; i_nd4++)
							{


								vecx = *(point_x + i_nd3) - *(point_x + i_nd4); //두 적혈구 포인트끼리의 거리
								if (vecx > safezone)
								{
									vecx -= XLENGTH;
								}
								if (vecx < (-safezone))
								{
									vecx += XLENGTH;
								}
								vecy = *(point_y + i_nd3) - *(point_y + i_nd4);
								
								if (vecy > safezone)
								{
									vecy -= YLENGTH;
								}
								if (vecy < (-safezone))
								{
									vecy += YLENGTH;
								}
								vecz = *(point_z + i_nd3) - *(point_z + i_nd4);





								aa = (vecx*vecx + vecy*vecy + vecz*vecz);
								c_length = pow(aa, 0.5);
								if (c_length < 3.5)
								{

									Fmag = 2.0 * De * s_beta * ((pow(2.748, 2.0*s_beta*(bumwe - c_length)) - pow(2.748, s_beta*(bumwe - c_length))));



									vecx /= c_length;
									vecy /= c_length;
									vecz /= c_length;

									//*(agforceX + i_nd3) -= Fmag * vecx;
									//*(agforceY + i_nd3) -= Fmag * vecy;
									//*(agforceZ + i_nd3) -= Fmag * vecz;

									//*(agforceX + i_nd4) += Fmag * vecx;
									//*(agforceY + i_nd4) += Fmag * vecy;
									//*(agforceZ + i_nd4) += Fmag * vecz;


									*(agforceX + i_nd3) =0;
									*(agforceY + i_nd3) =0;
									*(agforceZ + i_nd3) =0;
									*(agforceX + i_nd4) =0;
									*(agforceY + i_nd4) =0;
									*(agforceZ + i_nd4) =0;

								}

							}

						}
					}
				}
			}
		}
	}



	return;
}

void recovery(void)
{
	
	int     i_nd, i_nd2, ii, kk, kk1, i;
	double  c_length, aa, Fmag, vecx, vecy, vecz, bb, cc, disp, xp, yp, zp, shearmo, shearde;


#pragma omp parallel
{
#pragma omp	for	
		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{
			*(rforceX + i_nd) = 0.0;
			*(rforceY + i_nd) = 0.0;
			*(rforceZ + i_nd) = 0.0;

		}



#pragma omp for private(i_nd,kk,kk1,aa,c_length,Fmag,vecx,vecy,vecz,bb,cc,shearde, shearmo)	
		//for (i_nd2 = 1; i_nd2 <= rbcn; i_nd2++)			//for (i_nd = ((i_nd2 - 1)*l_main + 1); i_nd <= l_main*i_nd2; i_nd++)
		for (i_nd = 1; i_nd <= tp_num; i_nd++)
		{

				*(rforceX + i_nd) -= (*(point_x + i_nd) - *(point_xi + i_nd))*r_stif;
				*(rforceY + i_nd) -= (*(point_y + i_nd) - *(point_yi + i_nd))*r_stif;
				*(rforceZ + i_nd) -= (*(point_z + i_nd) - *(point_zi + i_nd))*r_stif;
			
		}
}//병렬화 종료

	return;
}


void	node_track (void)
{
	FILE *f_track;
	int i_nd, i_nd2, iii;
	double x, y,z, x0, y0, z0, x1, y1, z1, x2, y2, z2;
	double dis, vel, acc;

	f_track = fopen("f_track.txt", "a");

	if(Timek>=5) {x2=x1; y2=y1; z2=z1;}
	else {x2=0; y2=0; z2=0;}
	
	if(Timek>=3) {x1=x; y1=y; z1=z;}
	else {x1=0; y1=0; z1=0;}

	iii = 420;
	i_nd = iii;
	x = *(point_x + i_nd);
	y = *(point_y + i_nd);
	z = *(point_z + i_nd);

	x0 = *(point_xi + i_nd);
	y0 = *(point_yi + i_nd);
	z0 = *(point_zi + i_nd);

	
	//dis = ((x-x0)*100) * ((x-x0)*100);
	//dis += ((y-y0)*100) * ((y-y0)*100);
	//dis += ((z-z0)*100) * ((z-z0)*100);

	//dis = pow(dis, 0.5);
	//dis/=100;

	//dis2 = ((x-x1)*100) * ((x-x1)*100);
	//dis2 += ((y-y1)*100) * ((y-y1)*100);
	//dis2 += ((z-z1)*100) * ((z-z1)*100);
	
	//dis2 = pow(dis2, 0.5);
	//dis2/=100;

	fprintf(f_track, "%d\t%lf\t%lf\t%lf\n", Timek, x, y, z);

	fflush(f_track);
	fclose(f_track);


}




void oscillation(void)
{

	double height;
	double uFreq;
	int uT; //timestep for 1 cycle
	double pulseT;
	double uOffset; // 변경해야 하는 임시용
	height = ZLENGTH;
	uT = 10000;
	pulseT = Timek % uT;
	uFreq = 6.283 / (double) uT;
	uOffset = 15.38;

	uOsci = height * uMax * uFreq * cos(uFreq*pulseT);
	uOsci = uOsci * uOffset;
}

void	pulsatile	(void)
{
	/*
     pulse[0] =0.100904;
     pulse[1] =0.100904;
     pulse[2] =0.096386;
     pulse[3]=0.0512050;
     pulse[4] =0.0165660;
     pulse[5] =0.0662650;
     pulse[6] =0.632530;
     pulse[7] =1.0;
     pulse[8] =0.573795;
     pulse[9] =0.394578;
     pulse[10] =0.144578;
     pulse[11] =0;//=-0.03012;
     pulse[12] =0.085843;
     pulse[13] =0.114458;
     pulse[14] =0.138554;
     pulse[15] =0.11747;
     pulse[16] =0.054217;
     pulse[17] =0.069277;
     pulse[18] =0.103916;
     pulse[19] =0.100904;

	 */

	 pulse[0] =1;
     pulse[1] =0.975;
     pulse[2] =0.904;
     pulse[3] =0.793;
     pulse[4] =0.654;
     pulse[5] =0.5;
     pulse[6] =0.345;
     pulse[7] =0.206;
     pulse[8] =0.095;
     pulse[9] =0.024;
     pulse[10] =0;
     pulse[11] =-0.03012;
     pulse[12] =0.095;
     pulse[13] =0.206;
     pulse[14] =0.345;
     pulse[15] =0.5;
     pulse[16] =0.654;
     pulse[17] =0.793;
     pulse[18] =0.904;
     pulse[19] =0.975;
    return;
}





void	pulse_flow	(void)
{
	FILE *f_pulse;

	f_pulse = fopen("f_pulse.txt", "a");
	//int p_phase = 0;
	int pulse_s;
	double pulse_d;
	if (Timek%pulse_time==0) p_phase++; //pulse_time = 5000
	if (p_phase==pulse_step) p_phase = 0;
	


	//if(Timek<pulse_initial) u_pulse = uMax; 
	//else u_pulse = pulse[p_phase] * uMax;


	pulse_s = Timek%pulse_time;
	
	
	if (p_phase!=19) pulse_d = pulse[p_phase] + (pulse[p_phase+1]-pulse[p_phase])*(double(pulse_s)/double(pulse_time));
	else pulse_d = pulse[p_phase] + (pulse[0]-pulse[p_phase])*(double(pulse_s)/double(pulse_time));
	u_pulse = pulse_d * uMax;
	if(pulse_s%100 == 0) fprintf(f_pulse, "gg\t%d\t%d\t%lf\n", Timek, p_phase, u_pulse);
	
	/*
	if(Timek%20000<5000) u_pulse = uMax;
	else if(Timek%20000>10000 && Timek%20000<12000) u_pulse = uMax*0.2;
	else u_pulse =0;
	*/
	
	forcepowerAC = (u_pulse*(tau-0.5)/3.14/(pipe_r)/(pipe_r)*4.0*1.00*2.0);

	
	fflush(f_pulse);
	fclose(f_pulse);


    return;
}