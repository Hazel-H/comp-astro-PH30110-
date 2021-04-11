#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define res 1000

long long random_number(long long seed, int a, int m, int c)  //random number function
{
  return (seed*a + c) % m;
}

double probdist(double x)  //normalized probability distribution function 
{
  return (double) (3./8.)*(1.+(x*x));
}

double f(double y)  //inverse cumulative function where y is random number: 0 to 1.  y_in is scaled: 1/e to 1 so that corresponding tau is between 0 and 1 
{
  double y_min, y_length, y_in;
  double e = 2.7182818284590;
  y_min = (1./e);
  y_length = 1 - y_min; 
  y_in = (double)  y*y_length;
  return (double) -1.*log(1-y_in);
}

//INITIALIZE VARIABLES 
double x_reject, y_reject, fx;  //rejection method variables
double x_im, y_im, w; //importance sampling variables
int i, j; //loop counter
long long seed = 1; long long seed2 = 107;  long long seed3 = 313; long long m = 2147483647;  //variables for random number generation 


int part_1a() 
{ 
  FILE *fptr_reject;
  fptr_reject = fopen("rejection_method.text", "w");

  int N_accept = 0; //accepted x values counter
  int target = 5000; //target number of x values 
  int count = 0;  //count iterations for rejection method

  //***********************************************REJECTION METHOD*********************************************************
  while (N_accept<target)
  {  
    count++;  //counter for number of interations
    seed = random_number(seed, 16807, m, 0);      //generates random numbers from 0 to 1 
    seed2 = random_number(seed2, 48271, m, 0);       //generates random numbers from 0 to 1 

    x_reject = (double) seed*2./m -1.;  //x range of random numbers from -1 to +1
    y_reject = (double) seed2*0.75/m;  //y range of random numbers from 0 to peak (+0.75)
    fx = (double) probdist(x_reject);   //p(x)

    double x_accept[N_accept];   double y_accept[N_accept];   //create two arrays 

    if(y_reject<=fx)  
    {
      N_accept++;
      x_accept[N_accept] = x_reject;
      y_accept[N_accept] = y_reject; 
      fprintf(fptr_reject, "%.3f, %.3f\n", x_accept[N_accept], y_accept[N_accept]);   //print accepted (x,y) values to text file 
    }
  }

  //***********************************************IMPORTANCE SAMPLING*********************************************************
  int count_im = 0;   //count iterations for importance sampling method
  FILE *fptr_importance;
  fptr_importance = fopen("importance_sampling.text", "w");
  
  for(i=0; i<target; i++)
  { 
    count_im++;   //counter for number of interations
    seed = random_number(seed, 16807, m, 0);    
    seed2 = random_number(seed2, 48271, m, 0);
    
    double x_values[target];    double y_values[target];   //create two arrays for importance sampling method

    x_im = (double) seed*2./m -1.;  //x range of random numbers from -1 to +1
    y_im = (double) seed2*0.75/m;  //y range of random numbers from 0 to peak +0.75
    fx = (double) probdist(x_im);

    w = (double) fx/0.75;   //weighting factor 

    x_values[i] = x_im; 
    y_values[i] = (double) (y_im*w);

    fprintf(fptr_importance, "%.3f,%.3f\n", x_values[i], y_values[i]);    //print (x,y) values to text file 
  }
  printf("\nPart 1a\n");
  printf("number of iterations (rejection) = %d, number of iterations (importance) = %d\n", count, count_im);
  return 0;
}

//INITIALIZE VARIABLES 
int z_min = 0; //bottom of atmosphere
int z_max = 1; //top of atmosphere
double a = 1; //albedo = 1 meaning all scattering and no absorption
double x, y, z, dx, dy, dz, tau, rand1,rand2, theta, phi, mu; //variables
double tau_max = 10;  //total optical depth 
int N_photons = 1000000;  //required number of photons 
int N_bins = 10; //number of bins 

int part_1b()
{ 
  FILE *fptr;
  fptr = fopen("part_1b.text", "w");

  //create array of bins and set values to zero 
  int array[N_bins];
  for (j=0; j<N_bins; j++) 
  {
    array[j] = 0;
  }

  double d_theta = (double) (M_PI/2/N_bins);  //bin width 

  int count_below = 0;
  int count_above = 0;     //counters for number of photons above/below atmosphere
  while (count_above<N_photons) 
  {
    z = 0; x = 0; y = 0;  //each new photon starts at origin 
    theta = 0;  phi = 0;  mu = 1; //each new photon originally fired upwards 

      
    while (z>=z_min && z<=z_max) 
    {  //the photon has NOT escaped and is within atmosphere boundaries
      seed = random_number(seed, 16807, m, 0);
      rand1 = (double) seed/m;  //random number between 0 and 1 
      tau = (double) f(rand1);  //random distance tau between 0 to 1 from exponential func using cumulative method
      z = z + tau*1/tau_max;  //update positions 
      x = x + tau*1/tau_max; 
      y = y + tau*1/tau_max; 

      //should we scatter or absorb? 
      if(rand1 < a) 
      { 
        seed2 = random_number(seed2, 48271, m, 0);  //second random number for random new photon direction between 0 and pi 
        rand2 = (double) seed2/m ; //random number between 0 and 1 
        
        phi = 2*M_PI*rand2; 
        theta =  M_PI*rand2; 

        z = z + cos(theta)*1/tau_max;     //update positions 
        x = x + sin(theta)*cos(phi)*1/tau_max; 
        y = y + sin(theta)*sin(phi)*1/tau_max; 
      }

      else 
      {
        break;     //absorb and remove photon
      }  

    }

    if(z<z_min)
    {     
      count_below++; //photon escapes from below atmosphere - restart, new photon 
    }
    
    if(z>z_max)
    {   
      //photon escapes from above atmosphere - bin photon 
      int index = (int) (theta/d_theta);  //index of bin 
      array[index] = array[index]+1;  //add one for every photon escaping with angle within that bin 
      count_above++;
    }
  }
  printf("\nPart 1b\n");
  printf("No. photons escaped above: %d, No. photons escaped below: %d\n", count_above, count_below);

  for (j=0; j<N_bins; j++) 
    {
    double mu_mid = (double) (d_theta*j) + d_theta/2;    //angle from midpoint of bin
    double norm_int = (double) array[j]/N_photons;  //normalized intensity 
    fprintf(fptr, "%d, %e, %e\n", j, norm_int, mu_mid);   //print normalized intensity against midpoint of bin 
    }
    return 0; 
}


int part_1c()
{ 
  FILE *fptr;
  fptr = fopen("part_1c.text", "w");

  //create array of bins and set to zero 
  int array[N_bins];
  for (j=0; j<N_bins; j++) {
    array[j] = 0;
  }

  double d_theta = (double) (M_PI/2/N_bins);  //bin width 

  int count_below = 0;
  int count_above = 0;     

  while (count_above<N_photons) 
  {

    z = 0; x = 0; y = 0;  //each new photon starts at origin 
    theta = 0;  phi = 0;  mu = 1; //each new photon originally fired upwards 

      
    while (z>=z_min && z<=z_max) 
    {   //the photon has NOT escaped  
      seed = random_number(seed, 16807, m, 0);
      rand1 = (double) seed/m;  //random number between 0 and 1 
      tau = (double) f(rand1);  //random distance tau between 0 to 1 from exponential func using cumulative method
      z = z + tau*1/tau_max;   //update positions 
      x = x + tau*1/tau_max; 
      y = y + tau*1/tau_max; 

      //should we scatter or absorb? 
      if(rand1 < a) 
      {   
        int accept = 0; //counter- want one accepted value for mu (from rejection method) for each loop 
        while(accept<1)
        {
          seed2 = random_number(seed2, 48271, m, 0);
          seed3 = random_number(seed3, 69621, m, 0);  

          double x = (double) seed2*2./m -1.;  //x range of random numbers from -1 to +1
          double y = (double) seed3*0.75/m;   //y range of random numbers from 0 to 0.75 
          double fx = (double) probdist(x);   //p(x)

          if (y<fx)
          {
            mu = x;
            accept++;
          }

        }
        theta = acos(mu);
        z = z + cos(theta)*1/tau_max;    //update positions 
        x = x + sin(theta)*cos(phi)*1/tau_max; 
        y = y + sin(theta)*sin(phi)*1/tau_max; 
      }

      else 
      {
        break;     //absorb and remove photon
      }  
    }

    if(z<z_min)
    {     
      count_below++; //photon escapes from below atmosphere - restart, new photon 
    }
    
    if(z>z_max)
    {   
      //photon escapes from above atmosphere - bin photon 
      int index = (int) (theta/d_theta);  //index of bin 
      array[index] = array[index]+1;  //add one for every photon escaping with angle within that bin 
      count_above++;
    }

  }
  printf("\nPart 1c\n");
  printf("No. of photons escaped above: %d, No. of photons escaped below: %d\n\n", count_above, count_below);

  for (j=0; j<N_bins; j++)   //calculating normalized intensity and mid of bin of angle for plotting
  {
    double mu_mid = (double) (d_theta*j) + d_theta/2;    //angle from midpoint of bin
    double norm_int = (double) array[j]/N_photons;  //normalized intensity 
    fprintf(fptr, "%d, %e, %e\n", j, norm_int, mu_mid);   //print normalized intensity against midpoint of bin 
  }
  return 0; 
}


int main() 
{
  part_1a();
  part_1b();
  part_1c();
  return 0;
}

