/* Identify the largest cluster in each frame of a movie */
/* Revisions log:
 * 6.22.18 DM repurposing into a velocity field calculating.  removing clusters.
 * 6.02.16 include snapshot output of particles in the system
 * 5.13.16 default input name "smtest"
 * 7.17.15 Replacing ancient Stuart/Jared lookup table with my newly
           written version.
 * 7.15.15 Find overall largest cluster also.
 * 6.19.15 Writing a read_ascii subroutine to test float->double issues
 * 9.6.13 Well, that didn't work.  Next, I try using Hermann's method.
 * 8.29.13 Written */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define PI 3.14159265359
#define DEBUG 0
#define FILE_TYPE 0
#define ADRIAN 0

FILE *coarse_grain;

struct vortex{
  int id;
  int color;
  double x;
  double y;
  double fx; //values not in smtest.  measure needs to be done on ascii files
  double fy; //or during simulation 
  double radius;
  int clusterid;
};

struct parameters{

  int nV;
  int nP;
  int nV1, nV2;
  int maxnum;
  
  int maxnum_small;
  int maxnum_large;
  
  double radius_small;
  double radius_large;
  double runforce;
  int runtime;
  //
  double density;
  double phi2_small;
  double phi1_big;
  double pdensity;
  //
  double dt;
  int maxtime;
  int starttime;           //useful for restarting simulation from config file
  int writemovietime;
  double potential_radius; //change to pin_radius
  double potential_mag;    //change to pin_max_force;
  double kspring;
  double drive_mag;
  double drive_frq;
  int decifactor;

  //change the driving force at a given time by a given amount.
  int drive_step_time;
  double drive_step_force;

  int restart;
};

struct syssize{
  double SX;
  double SY;
  double SX2;
  double SY2;
};


void coarse_grain_field(struct vortex *vortex, int num_vor,
			struct syssize syssize, int time,
			struct parameters parameters,
			FILE *fpt_enstrophy);
  
//DM read in Pa0 
void get_parameters_file(struct parameters *parameters, struct syssize *syssize);
 
int read_frame(FILE *in,int time, int *nV, struct vortex *vortex);


int main(int argc,char *argv[])
{
  FILE *in,*out;
  FILE *enstrophy_file;
 
  struct syssize syssize;

  int completeframe;

  struct vortex *vortex;
  struct parameters parameters;    //nearly everything from Pa0

  int maxnum;

  //number of vortices
  int nV,time,nVold;
  int i,j,k,id,color;
  double x,y,rad;

  double SX,SY,SX2,SY2;
  double dx,dy,dist;

  int frame = 0;
  int avgmax;
  int maxcount;

  get_parameters_file(&parameters, &syssize);
  
  vortex=malloc((parameters.maxnum+1)*sizeof(struct vortex));

  //open and close file to clear it:
  if ((enstrophy_file = fopen("enstrophy_file", "w")) == NULL){
    printf("Can't open %s.\n","enstrophy_file");
    exit(1);
  }

  
  //loop through all of the printed times from writemovietime
  for(time=parameters.writemovietime;
      time < parameters.maxtime;
      time+= parameters.writemovietime){
	
    frame++;
    nVold=nV;

    //read in the ascii file
    completeframe=read_frame(in,time,&nV,vortex);

    //accumulate the velocities into a field
    if(DEBUG == 2) printf("time %d\n",time);
    coarse_grain_field(vortex,nV,syssize,time,parameters,enstrophy_file);

  }//end while loop
  
  printf("Reached frame %d\n",frame);

  fclose(enstrophy_file);

  return 0;
  
}//end main function

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
int read_frame(FILE *in,int time, int *nV, struct vortex *vortex)
{
  int i;
  float xin,yin,zin;
  float fxin, fyin, speedin;
  
  char ascii_file[120]; //="velocity_data/XV_data_t=";
  char str_time[10];

  int n_scan;
  char str1[120], str2[120], str3[120];
  char str4[120], str5[120], str6[120];
  char str7[120], str8[120];
  int int_junk;
  
  //set a new filename
  strcpy(ascii_file,"velocity_data/XV_data_t=");
  sprintf(str_time,"%08d",time); //convert current to a string   
  strcat(ascii_file,str_time);
  if(DEBUG) printf("%s\n",ascii_file);

  if((in=fopen(ascii_file,"r"))==NULL){
    printf("Error opening file %s\n",ascii_file);
    exit(-1);
  }
  else if(DEBUG == 2){
    printf("Opened %s\n",ascii_file);
  }

  //read header
  /*
    #Number of Vortices: 1075
    #N_{small}: 531
    #N_{big}: 545
    #time frame: 50
    #id type x y fx fy radius speed
   */

  n_scan = fscanf(in,"%s %s %s %d\n",str1,str2,str3,nV);
  if(FILE_TYPE == 1){
    n_scan = fscanf(in,"%s %d\n",str1,&int_junk);
    n_scan = fscanf(in,"%s %d\n",str1,&int_junk);
  }
  n_scan = fscanf(in,"%s %s %d\n",str1,str2,&time);
  n_scan = fscanf(in,"%s %s %s %s %s %s %s %s\n",str1,str2,str3,str4,str5,str6,str7,str8);

    if(DEBUG == 2){
      printf("%d %d \n",*nV,time);
      fflush(stdout);
  }
    
  if(feof(in)) return 0;
  for(i=0;i<*nV;i++){

    if(DEBUG == 2) printf("%d %d\n",i,*nV);
    
    n_scan = fscanf(in,"%d %d %f %f %f %f %f %f\n",&(vortex[i].id),
		    &(vortex[i].color), &xin, &yin, &fxin, &fyin, &zin, &speedin);
    
    //error checking - spot check these three values, should be enough
    //more extensive checking would ensure the forces aren't unusally large
    //and the radius is about the same size
    if( vortex[i].id < 0 || vortex[i].id > *nV){
      printf("Wrong number of vortices.  Vortex id: %d.",vortex[i].id);
      printf("Check your file types (FILE_TYPE = %d  \n",FILE_TYPE);
      printf("There are different headers in velocity_data/XV... files\n");
      exit(-1);
    }
    if( xin < 0.0 || xin > 100000.0){
      printf("x value too large.  x: %f.",xin);
      printf("Check your file types.\n");
      exit(-1);
    }
    if( yin < 0.0 || yin > 100000.0){
      printf("y value too large.  y: %f.",yin);
      printf("Check your file types.\n");
      exit(-1);
    }
    
    //assign values to vortex array
    vortex[i].x=(double)xin;
    vortex[i].y=(double)yin;
    vortex[i].fx=(double)fxin;
    vortex[i].fy=(double)fyin;
    vortex[i].radius=(double)zin;
    vortex[i].clusterid=-1;

    if(DEBUG == 2){
      printf("%f %f %f %f %f \n",vortex[i].x, vortex[i].y, vortex[i].fx, vortex[i].fy, vortex[i].radius);
      fflush(stdout);
    }
      
  }

  fclose(in);
  return 1;
  
}
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

//-----------------------------------------------------------------
void get_parameters_file(struct parameters *parameters,
			 struct syssize *syssize)
{
  FILE *in;
  char trash[120];
  double cellsize,length_scale;
  double resolution;

  int num_scan;
  
  resolution=1e-6;

  if((in=fopen("Pa0","r"))==NULL){
    printf("Input file Pa0 not found\n");
    exit(-1);
  }
  else if (DEBUG) printf("reading Pa0.\n");
  fflush(stdout);

  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).density));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).phi2_small));
  if(ADRIAN == 0){
    num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).phi1_big));
  }
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).pdensity));

  //---------------------------------------------------------------------
  
  if (DEBUG){
    printf("reading parameters \n density %f, \n pin density %f.\n", (*parameters).density, (*parameters).pdensity);
    fflush(stdout);
  }

  //-
  num_scan = fscanf(in,"%s %lf\n",trash,&((*syssize).SX));
  (*syssize).SX2=(*syssize).SX*0.5;
  //--
  num_scan = fscanf(in,"%s %lf\n",trash,&((*syssize).SY));
  (*syssize).SY2=(*syssize).SY*0.5;
  
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).radius_small));
  
  if(ADRIAN == 0){    
    num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).radius_large));
  }
  else{
    (*parameters).radius_large = 1.0;
  }

  if (DEBUG){
    printf("SX %f SY %f, \nr_s %f r_b %f\n", (*syssize).SX,
	   (*syssize).SY,
	   (*parameters).radius_small,(*parameters).radius_large);
    fflush(stdout);
  }
    
  //--where it is appropriate to accurately calculate maxnum
  //--based on particle sizes

  double prefix = (*syssize).SX*(*syssize).SY/(PI*(*parameters).radius_small*(*parameters).radius_small);
   
  (*parameters).maxnum= (int) ceil(prefix * (*parameters).density);
  if(DEBUG) printf("\n maxnum is: %d\n", (*parameters).maxnum);
  //--
  
  num_scan = fscanf(in,"%s %d\n",trash,&((*parameters).runtime));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).runforce));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).dt));
  num_scan = fscanf(in,"%s %d\n",trash,&((*parameters).maxtime));
  num_scan = fscanf(in,"%s %d\n",trash,&((*parameters).writemovietime));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).kspring));
  
  num_scan = fscanf(in,"%s %lf\n",trash,&cellsize); // size of lookup cell
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).potential_radius));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).potential_mag));
  
  num_scan = fscanf(in,"%s %lf\n",trash,&length_scale);
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).drive_mag));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).drive_frq));
  num_scan = fscanf(in,"%s %d\n",trash,&((*parameters).decifactor));
  
  //starting from a file?
  num_scan = fscanf(in,"%s %d\n",trash,&((*parameters).restart));

  //ramping the drive rate?
  num_scan = fscanf(in,"%s %d\n",trash,&((*parameters).drive_step_time));
  num_scan = fscanf(in,"%s %lf\n",trash,&((*parameters).drive_step_force));

  //don't bother setting lookupdata - no particle in cell method here.
  //the velocity grid is hardcoded
  
  fclose(in);
}
 
//not used
void distance(double *dr,double *dx,double *dy,double x1,double y1,
	      double x2,double y2,struct syssize syssize)
{
  double locdx,locdy;

  locdx=x1-x2;
  if(locdx>syssize.SX2) locdx-=syssize.SX;
  if(locdx<=-syssize.SX2) locdx+=syssize.SX;
  locdy=y1-y2;
  if(locdy>syssize.SY2) locdy-=syssize.SY;
  if(locdy<=-syssize.SY2) locdy+=syssize.SY;
  *dr=sqrt(locdx*locdx+locdy*locdy);
  *dx=locdx;
  *dy=locdy;
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/



/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

void coarse_grain_field(struct vortex *vortex,
			int num_vor,
			struct syssize syssize,
			int time,
			struct parameters parameters,
			FILE *enstrophy_file){

  
  //Make a grid that simply divides up the space:
  float length_scale = 1.2;  //large bins for testing
  int i, j, k;               //for loop control variable

  int m, n, p, q;            //for grid counting
  
  //file to print to, hardwired name
  char str_time[10];
  char coarse_grain_file[120] = "field_data_";
  
  if(time == (parameters.maxtime - parameters.writemovietime)){
    sprintf(str_time,"%08d",time); //convert current to a string  
    strcat(coarse_grain_file,str_time);
  }

  //for readability
  double xmax = syssize.SX;
  double ymax = syssize.SY;

  //calculate dvx/dy and dvy/dx for the curl
  double diff_vx_dy;
  double diff_vy_dx;

  //calculate how many bins in x and y based on the length_scale (hardwired)
  int x_bins = (int)(xmax / length_scale);
  int y_bins = (int)(ymax / length_scale);

  //if (DEBUG){
  //printf("%d \t %d \n",x_bins,y_bins);
  //exit(0);
  //}

  //define two arrays to hold the net vx and vy in each bin
  double vx_field[x_bins][y_bins];
  double vy_field[x_bins][y_bins];
  double curl_field[x_bins][y_bins];
  double enstrophy_field[x_bins][y_bins];
  double system_enstrophy = 0.0;

  //iterate over total number of bins, zeroing the values
   for(j=0;j<x_bins;j++){
     for(k=0;k<y_bins;k++){
       
       vx_field[j][k]=0.0;
       vy_field[j][k]=0.0;
       
       curl_field[j][k]=0.0;
       
       enstrophy_field[j][k]=0.0;

     }
   }

   //now iterate over the individual particle
   //i.e. values of x,vx and y,vy, 
   //writing them to the fields value by value

   for(i=0;i<num_vor;i++){

     //identify x-bin location
     j = (vortex[i].x/xmax)*x_bins;
     
     //identify y-bin location
     k = (vortex[i].y/ymax)*y_bins;     

     //make sure the bin is within the array vx_field
     if (j >= 0 && j <= x_bins && k >= 0 && k <= y_bins) {
       vx_field[j][k] += vortex[i].fx;
       vy_field[j][k] += vortex[i].fy;
     }
     else{
       printf("j=%d,xbins=%d,k=%d,ybins=%d,x=%f,y=%f",j, 
	      x_bins,k,y_bins,j*length_scale,k*length_scale);
       
       printf("\n Something bug in your data or coarse graining method\n");
       exit(-1); 
     }
   }

   //fields written, calculate the curl over the coarse grained grid, 
   //taking periodic boundary conditions into account

     for(j=0;j<x_bins;j++){
       for(k=0;k<y_bins;k++){

	 m = j + 1;  //high value for x, if not edge
	 n = j - 1;  //low value for x, if not edge

	 p = k + 1;  //high value for y, if not edge
	 q = k - 1;  //low value for y, if not edge

	 if (j == 0)
	   n = x_bins-1; //low value is a periodic wrap
	 else if ( j == x_bins-1)
	   m = 0;       //high value is periodic wrap

	 if (k == 0)
	   p = y_bins-1; //low value is a periodic wrap
	 else if ( k == y_bins-1)
	   q = 0;       //high value is periodic wrap

	 //Calculate the spatial derivatives via the two grid points 
	 //nearest that grid point in question
	 diff_vx_dy = (vx_field[j][p] - vx_field[j][q])/(2*length_scale);
	 diff_vy_dx = (vy_field[m][k] - vy_field[n][k])/(2*length_scale);

	 //simply take the difference of the spatial derivatives
	 //\omega = curl V = \nabla \times \vec{v} = dv_x/dy - dv_y/dx
	 curl_field[j][k] = diff_vx_dy - diff_vy_dx;

	 //\[enstrophy] = \Epsilon = (1/2)*\int_{S=Area} (\omega^2 dS)

	 double dS = (2*length_scale)*(2*length_scale);
	 enstrophy_field[j][k] = 0.5*(curl_field[j][k])*(curl_field[j][k])*dS;

	 if( enstrophy_field[j][k] > 1000.0){
	   printf("Enstrophy j,k=%d,%d, %f\n",j,k,enstrophy_field[j][k]);
	   printf("System Enstrophy = %f\n",system_enstrophy);
	   exit(0);

	 }
	 system_enstrophy += enstrophy_field[j][k];
       }
     }
		    
     /*------------------------------------------------------------------*/
     /*------------------Open and print to curl/field file---------------*/
     /*------------------------------------------------------------------*/
     

     //print the data to a file

     if(time == (parameters.maxtime - parameters.writemovietime)){

       if ((coarse_grain = fopen(coarse_grain_file, "w")) == NULL){
	 printf("Can't open %s.\n",coarse_grain_file);
	 exit(-1);
       }

       if(DEBUG == 1){
	 printf("file number %d\n", fileno(coarse_grain));
	 fflush(stdout);
       }
       fprintf(coarse_grain,"#%d \n #%d \n #%f \n",
	       num_vor,time,system_enstrophy);   

       for(j=0;j<x_bins;j++){
	 for(k=0;k<y_bins;k++){
	 
	   double x = ((double)j + 0.5)*length_scale;
	   double y = ((double)k + 0.5)*length_scale;

	   fprintf(coarse_grain,"%f\t %f\t %3.6e\t %3.6e\t %3.6e\t %3.6e\n", \
		   x, y,						\
		   vx_field[j][k],  vy_field[j][k],			\
		   curl_field[j][k], enstrophy_field[j][k]);
	 
	 }
       }
       fflush(coarse_grain);
     
       fclose(coarse_grain);
     }
       
     /*------------------------------------------------------------------*/
     /*------------- Print to the enstrophy frame file ------------------*/
     /*------------------------------------------------------------------*/
    

     fprintf(enstrophy_file,"%d \t %f \n", time, system_enstrophy);
     //printf("%d \t %f \n", time, system_enstrophy);
     //fclose(enstrophy_file);

   return;
}//end subroutine coarse grain
