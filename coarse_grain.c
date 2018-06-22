#include "Movie_Conv.h"

void coarse_grain_field(int num_vor, int time,
			float *x, float *y,
			float *vx, float *vy, int print_flag){


  //Make a grid that simply divides up the space:
  float length_scale = 1.2;  //large bins for testing
  int i, j, k;               //for loop control variable

  int m, n, p, q;            //for grid counting
  
  //file to print to, hardwired name
  char coarse_grain_file[50] = "coarse_grain_field.txt";

  //calculate how many bins in x and y based on the length_scale
  int x_bins = (int)(xmax / length_scale);
  int y_bins = (int)(ymax / length_scale);

  //printf("%d \t %d \n",x_bins,y_bins);
  //exit(0);

  //define two arrays to hold the net vx and vy in each bin
  double vx_field[x_bins][y_bins];
  double vy_field[x_bins][y_bins];
  double curl_field[x_bins][y_bins];
  double enstrophy_field[x_bins][y_bins];
  double system_enstrophy = 0;

  //iterate over total number of bins, zeroing the values
   for(j=0;j<x_bins;j++){
     for(k=0;k<y_bins;k++){
       vx_field[j][k]=0.0;
       vy_field[j][k]=0.0;
       curl_field[j][k]=0.0;
       enstrophy_field[j][k]=0.0;

     }
   }

   //now iterate over the individual values of vx and vy, 
   //writing them to the fields value by value

   for(i=0;i<num_vor;i++){

     //identify x-bin
     j = (x[i]/xmax)*x_bins;
     //identify y-bin
     k = (y[i]/ymax)*y_bins;     

     //make sure the bin is within the array vx_field
     if (j >= 0 && j <= x_bins && k >= 0 && k <= y_bins) {
       vx_field[j][k] += vx[i];
       vy_field[j][k] += vy[i];
     }
     else{
       printf("j=%d,xbins=%d,k=%d,ybins=%d,x=%f,y=%f",j,x_bins,k,y_bins,j*length_scale,k*length_scale);
       printf("\n Something buggy in your coarse graining\n");
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
	 double diff_vx_dy = (vx_field[j][p] - vx_field[j][q])/(2*length_scale);
	 double diff_vy_dx = (vy_field[m][k] - vy_field[n][k])/(2*length_scale);

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
     /*------------------Open and print to enstrophy file----------------*/
     /*------------------------------------------------------------------*/
     if(print_flag){

       //print the data to a file
       if ((coarse_grain = fopen(coarse_grain_file, "w")) == NULL){
	 printf("Can't open %s.\n",coarse_grain_file);
	 exit(-1);
       }

       fprintf(coarse_grain,"#%d \n #%d \n #%f \n",num_vor,time,system_enstrophy);   

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
   

       fclose(coarse_grain);
     }
     /*------------------------------------------------------------------*/
     /*------------- Close the enstrophy frame file ---------------------*/
     /*------------------------------------------------------------------*/
    

     fprintf(enstrophy_file,"%d \t %f \n", time, system_enstrophy);

   return;
}//end subroutine read_single_ascii




