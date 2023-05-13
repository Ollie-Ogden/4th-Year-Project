#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define PI 3.141592

// Creating the 3D DFT function
void nudft3d(double complex *freq, double complex *image, int N, double uvw[N][3], double kx[N], double ky[N], double kz[N]){
   #pragma omp parallel for
    for (int u = 0; u < N; u++)
        {
        for (int x = 0; x < N; x++)
            {
            for (int y = 0; y < N; y++)
                {
                    for (int z=0; z<N; z++)
                        {
                        double coef = ((uvw[u][0] * kx[x])/ (double)N) + ((uvw[u][1] * ky[y]) / (double)N) + ((uvw[u][2] * kz[z]) / (double)N);
                        freq[u] += (image[x * N*N + y*N +z] * cexp(-2.0 * M_PI * I * coef ));
                        }
                    }
                }
            }
        }



// Creating the 3D iDFT function
void inudft3d(double complex *reconstructed_image, double complex *freq, int N, double uvw[N][3], double kx[N], double ky[N], double kz[N]){
  #pragma omp parallel for
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            for (int z=0; z<N; z++)
            {
                for (int u = 0; u < N; u++)
                {
                    double coef =((uvw[u][0] * kx[x])/ (double)N) + ((uvw[u][1] *ky[y]) / (double)N)+ ((uvw[u][2] * kz[z]) / (double)N);
                    reconstructed_image[x*N*N+y*N +z]+=(freq[u] * cexp(2.0 * M_PI * I * coef ));
                }
                    reconstructed_image[x*N*N+ y*N +z]= reconstructed_image[x*N*N +y*N+z]/pow(N,3.0);

            }
        }
    }
}

int main() {
    //Initialising Variables
    int N=100;
    int Cols=3;
    double complex reconstructed_image[N * N* N];
    double kx[N];
    double ky[N];
    double kz[N];
    double start = -175.5;
    double end = 174.5;
    double step = (end - start) / (N - 1);
   
    //nuDFT indexing for x y and z
    for (int i = 0; i < N; i++) {
        kx[i] = start + i * step;
    }

    for (int i = 0; i < N; i++) {
        ky[i] = start + i * step;
    }

    for (int i = 0; i < N; i++) {
        kz[i] = start + i * step;
    }
   
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
            for(int k=0; k<N; k++){
                reconstructed_image[i * N*N+j*N+k] = 0;
            }
        }
    }
   
   //Getting the initial dirty image
    inudft3d(reconstructed_image, visibilities, N, uvw, kx,ky,kz);

    double complex max_val = reconstructed_image[0];
    int max_index = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int index = i*N*N + j*N + k;
                if (creal(reconstructed_image[index]) > creal(max_val)) {
                    max_val = reconstructed_image[index];
                    max_index = index;
                }
            }
        }
    }

    //Initialising an array of zeros to input our maximum value into
    double complex image_max[N* N* N];
        for (int i=0; i<N*N*N; i++){
        image_max[i]=0;
    }

    //Placing the maximum value into the array of zeros at the right index
    image_max[max_index]=max_val;

    //Initialising an array of zeros for our CLEAN'ed points to get put into
    double complex image_clean[N * N* N];
        for (int i=0; i<N*N*N; i++){
        image_clean[i]=0;
    }

    //Placing our first value into the CLEAN image
    image_clean[max_index]=max_val;

    //Finding the DFT of our zeros array with the first max value in it
    double complex frequency_max[N ];
    for (int i=0; i<N; i++){
        frequency_max[i]=0;
    }
    nudft3d(frequency_max, image_max, N, uvw, kx, ky ,kz);

    // Creating a new array where we are taking the original sampled freq and subtracting the sampled position DFT of the maximum value
    double complex frequency_new[N];
        for (int i=0; i<N; i++){
        frequency_new[i]=0;
    }
    for(int i=0; i<N; i++){
        frequency_new[i]=visibilities[i]-frequency_max[i];
    }

    //Finding the iDFT of the new frequency to find the next maximum value
    double complex image_new[N*N*N];
    inudft3d(image_new,frequency_new, N, uvw, kx,ky,kx);

    double complex max_val1 = image_new[0];
    int max_index1 = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int index = i*N*N + j*N + k;
                if (creal(image_new[index]) > creal(max_val1)) {
                    max_val1 = image_new[index];
                    max_index1 = index;
                }
            }
        }
    }

    //Placing our second value into the CLEAN image
    image_clean[max_index1]=max_val1;

    int p=0;

     while (p<1000){
    //while (creal(image_new[maxindex_new])>0.85*creal(reconstructed_image[maxindex])){
        //Initialising an array of zeros to input our maximum value into

        for (int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                        for(int k=0; k<N;k++){
                        image_max[i*N*N +j*N+k]=0;
                        }
                }
        }

        //Placing the maximum value into the array of zeros at the right index
        image_max[max_index1]=max_val1;
        
        //Getting the new max frequency array
        nudft3d(frequency_max, image_max,N, uvw,kx,ky,kz);
        
        //Subtracting the maximum frequency array from the residual frequency from the previous iteration
        for (int i=0; i<N; i++){
            frequency_new[i]=frequency_new[i]-frequency_max[i];
                        }
        
        //Getting the new residual image
        inudft3d(image_new,frequency_new, N, uvw,kx,ky,kz);

        //Finding the maximum of the residual image
        max_val1 = image_new[0];
        max_index1 = 0;

        int max_indices[N*N*N];
        int num_max = 0;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    int index = i*N*N + j*N + k;
                    if (creal(image_new[index]) > creal(max_val1)) {
                    max_val1 = image_new[index];
                    max_index1 = index;
                    }
                }
            }
        }
        
        //Placing the maxmimum in the clean image array
        image_clean[max_index1]+=max_val1;
        p=p+1;
        printf("The iteration number is %d \n", p);

    }


    printf("\n");
    printf("%0.5f, %d", creal(max_val), max_index);
        printf("\n");
    printf("%0.5f, %d", creal(max_val1), max_index1);


    double complex snapshot[N][N][N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for(int k=0;k<N; k++){
                int index = i * N*N+j*N+k ;
                snapshot[i][j][k]=image_clean[index];
            }
        }
    }

    fp = fopen("snapshot3DCLEAN.txt", "w");
    // Write the CLEAN image array to the file to be plotted in gnuplot or MATLAB
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                fprintf(fp, "%.9f ", creal(snapshot[i][j][k]));
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }

    // Close the file
    fclose(fp);

}
