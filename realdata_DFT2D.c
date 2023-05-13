#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define PI 3.141592

// Creating the 3D DFT function
void dft2d(double complex *freq, double complex *image, int N, double uvw[N][3]){
   #pragma omp parallel for
    for (int u = 0; u < N; u++)
        {
        for (int x = 0; x < N; x++)
            {
            for (int y = 0; y < N; y++)
                {
                    double coef = ((uvw[u][0] * (double)x)/ (double)N) + ((uvw[u][1] * (double)y) / (double)N);
                    double complex test=-2.0 * M_PI * I * coef;
                    freq[u] += (image[x * N + y] * cexp(-2.0 * M_PI * I * coef ));
                    }
                }
            }
        }



// Creating the 3D iDFT function
void idft2d(double complex *reconstructed_image, double complex *freq, int N, int Cols, double uvw[N][Cols]){
  #pragma omp parallel for
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
                for (int u = 0; u < N; u++)
                    {
                    double coef =((uvw[u][0] * (double)x)/ (double)N) + ((uvw[u][1] *(double)y) / (double)N);
                    double complex test=2.0 * M_PI * I * coef;
                    reconstructed_image[x*N+y]+=(freq[u] * cexp(2.0 * M_PI * I * coef ));
                    }
                    reconstructed_image[x*N+ y]= reconstructed_image[x*N +y]/pow(N,2.0);

                }
        }
    }





int main() {
    //Initialising Variables
    int N=351;
    int Cols=3;
    double complex reconstructed_image[N * N];

    //Loading in the UVW samples
    FILE *fp;
    double uvw[N][Cols];
    int row, col;

    fp = fopen("uvw.txt", "r");


    for (row = 0; row < N; row++) {
        for (col = 0; col < Cols; col++) {
            fscanf(fp, "%lf,", &uvw[row][col]);
        }
    }

    fclose(fp);

    //Loading in the visibilities
    double real, imag;
    double complex visibilities[N];

    fp = fopen("vis1.txt", "r");
    for (int i = 0; i < N && fscanf(fp, "%lf %lf", &real, &imag) != EOF; i++) {
        visibilities[i] = real + imag * I;
    }

    fclose(fp);

    //Initialising my reconstructed image array of zeros
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
                reconstructed_image[i * N+j] = 0;
            }
        }

   //Finding the initial dirty image
    idft2d(reconstructed_image, visibilities, N, Cols, uvw);
   
   // To find the maximum point in the dirty image
    double complex max_value;

    int maxindex=0;
    for (int i = 1; i < N*N; ++i)
    {
        if (creal(reconstructed_image[maxindex]) < creal(reconstructed_image[i]))
            maxindex = i;
    }


    //Initialising an array of zeros to input our maximum value into
    double complex image_max[N * N];
        for (int i=0; i<N*N; i++){
        image_max[i]=0;
    }

    //Placing the maximum value into the array of zeros at the right index
    image_max[maxindex]=reconstructed_image[maxindex];

    //Initialising an array of zeros for our CLEAN'ed points to get put into
    double complex image_clean[N * N];
        for (int i=0; i<N*N; i++){
        image_clean[i]=0;
    }

    //Placing our first value into the CLEAN image
    image_clean[maxindex]=reconstructed_image[maxindex];

    //Finding the DFT of our zeros array with the first max value in it
    double complex frequency_max[N];
    dft2d(frequency_max, image_max, N, uvw);

    // Creating a new array where we are taking the original sampled freq and subtracting the sampled position DFT of the maximum value
    double complex frequency_new[N];
    for(int i=0; i<N; i++){
        frequency_new[i]=visibilities[i]-0.2*frequency_max[i];
    }


    //Finding the iDFT of the new frequency to find the next maximum value
    double complex image_new[N*N];
    idft2d(image_new,frequency_new, N, Cols, uvw);

    // Calling the maximum value function
    double complex max_valuenew;
    int maxindex_new=0;

    for (int i = 1; i < N*N; ++i)
    {
        if (creal(image_new[maxindex_new]) < creal(image_new[i]))
            maxindex_new = i;
    }

    //Placing our second value into the CLEAN image
    image_clean[maxindex_new]=image_new[maxindex_new];
    int p=0;
    while (p<800){
  //while (creal(image_new[maxindex_new])>0.28*creal(reconstructed_image[maxindex])){

        //Initialising an array of zeros to input our maximum value into
        for (int i=0; i<N*N; i++){
            image_max[i]=0;
        }

        //Placing the maximum value into the array of zeros at the right index
        image_max[maxindex_new]=image_new[maxindex_new];
   
       //getting the max frequency array
        dft2d(frequency_max, image_max, N, uvw);
     
       //Subtracting to get the new residual frequency array
        for(int i=0; i<N; i++){
            frequency_new[i]=frequency_new[i]-0.2*frequency_max[i];
        }
       
       //New residual image
        idft2d(image_new,frequency_new, N, Cols, uvw);

       //Finding maximum for new residual image
        maxindex_new=0;

        for (int i = 1; i < N*N; ++i)
        {
            if (creal(image_new[maxindex_new]) < creal(image_new[i]))
            maxindex_new = i;
        }
        image_clean[maxindex_new]+=image_new[maxindex_new];
        p=p+1;
        printf("The iteration number is %d \n", p);
    }
   
   
    double complex snapshot[N][N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            int index = i * N+j;
            snapshot[i][j]=reconstructed_image[index];
            }
        }

    fp = fopen("snapshot2DDFT.txt", "w");
 // Write the reconstructed image array to the file to be plotted in MATLAB
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fp, "%.5f ", creal(snapshot[i][j]));
        }
        fprintf(fp, "\n");
    }

    // Close the file
    fclose(fp);

    double complex snapshotCLEAN[N][N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            int index = i * N+j;
            snapshotCLEAN[i][j]=image_clean[index];
            }
        }

    fp = fopen("snapshot2DCLEANDFT.txt", "w");
 // Write the CLEAN image array to the file to be plotted in MATLAB
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fp, "%.5f ", creal(snapshotCLEAN[i][j]));
        }
        fprintf(fp, "\n");
    }

    // Close the file
    fclose(fp);

}
