#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#define PI 3.141592

// Creating the 3D DFT function
void dft3d(double complex *freq, double complex *image, int N){
    #pragma omp parallel for
    for (int u = 0; u < N; u++)
    {
        for (int v = 0; v < N; v++)
        {
            for (int w=0; w<N; w++)
                {
                for (int x = 0; x < N; x++)
                    {
                    for (int y = 0; y < N; y++)
                        {
                        for (int z=0; z< N; z++)
                            {
                            freq[u*N*N + v*N + w] += image[x * N* N + y*N + z] * cexp(-2 * PI * I * (u * x/ (double)N + v * y / (double)N + w * z / (double)N));
                            }
                        }
                    }
                }
            }
        }
    }


// Creating the 3D iDFT function
void idft3d(double complex *reconstructed_image, double complex *freq, int N){
    #pragma omp parallel for
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            for (int z=0; z<N; z++)
                {
                for (int u = 0; u < N; u++)
                    {
                    for (int v = 0; v < N; v++)
                        {
                        for (int w=0; w < N; w++)
                            {
                            reconstructed_image[x*N*N +y*N +z]+=freq[u*N*N + v*N + w] * cexp(2 * PI * I * (u * x/ (double)N + v * y / (double)N + w * z / (double)N));
                            }
                        }
                    }
                    reconstructed_image[x*N*N +y*N +z]=reconstructed_image[x*N*N +y*N +z]/pow(N,3);
                }
            }
        }
    }

#define N 20
#define RADIUS1 4
#define RADIUS2 2
#define RADIUS3 3

typedef struct {
    double x;
    double y;
    double z;
} point;

double distance(point a, point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}


int main() {
    // Initialising Variables
    int i, j, k, idx;
    double x, y, z;
    point center1 = {N/2 - 2, N/2-4, N/2};
    point center2 = {N/2 + 5, N/2 + 4, N/2};
    point center3 = {N/2, N/2 + 5, N/2 + 6};
    double complex matrix[N * N * N];
    double complex sphere[N*N*N];
    double complex simulated_data[N * N * N];
    double complex reconstructed_image[N * N * N];
    double complex random_samples[N * N * N];
    double complex freq[N * N * N];

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                idx = i * N * N + j * N + k;
                x = i;
                y = j;
                z = k;
                double d1 = distance(center1, (point) {x, y, z});
                double d2 = distance(center2, (point) {x, y, z});
                double d3 = distance(center3, (point) {x, y, z});
                if (d1 <= RADIUS1) {
                    matrix[idx] = 1;
                } else if (d2 <= RADIUS2) {
                    matrix[idx] = 1;
                } else if (d3 <= RADIUS3) {
                    matrix[idx] = 1;
                } else {
                    matrix[idx] = 0;
                }
            }
        }
    }

    //Initialising my Frequency array of zeros
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k=0; k<N; k++){
                freq[i * N*N + j*N+k] = 0;
            }
        }
    }

    // DFT of the simulated data to get my frequency data
    dft3d(freq, matrix,N);

    //Initialising my reconstructed image array of zeros
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k=0; k<N; k++){
                random_samples[i * N*N + j*N+k] = 0;
            }
        }
    }

    // Generating a random number of samples from the frequency domain
    int lower = 0, upper = pow(N,3), count = (0.25*(pow(N,3)));


    srand(time(0));
    i = 0;
    int samples[count];
        for (int i=0; i<count; i++){
            samples[i]=0;
    }

    while(i < count) {
        int num = (rand() % (upper - lower + 1)) + lower;
        if(random_samples[num] == 0) {
            random_samples[num]=freq[num];
            samples[i]=num;
            i++;
        }
    }

    //Initialising my reconstructed image array of zeros
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k=0; k<N; k++){
                reconstructed_image[i * N*N + j*N+k] = 0;
            }
        }
    }
    // Finding the inverse 3d DFT to get back to my original array
    idft3d(reconstructed_image, random_samples, N);

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
    double complex frequency_max[N * N* N];
    for (int i=0; i<N*N*N; i++){
        frequency_max[i]=0;
    }
    dft3d(frequency_max, image_max, N);

    //Creates an array of zeros and then inputs only the randomly sampled frequency points back into the array
    double complex Fmax[N*N*N];
    for (int i=0; i<N*N*N; i++){
        Fmax[i]=0;
    }
    int sample_pos;
    for(int i=0; i<count; i++){
        Fmax[samples[i]]=frequency_max[samples[i]];
    }

    // Creating a new array where we are taking the original sampled freq and subtracting the sampled position DFT of the maximum value
    double complex frequency_new[N*N*N];
        for (int i=0; i<N*N*N; i++){
        frequency_new[i]=0;
    }
    for(int i=0; i<N*N*N; i++){
        frequency_new[i]=random_samples[i]-Fmax[i];
    }

    //Finding the iDFT of the new frequency to find the next maximum value
    double complex image_new[N*N*N];
    idft3d(image_new,frequency_new, N);

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
    int maxindexarray[N*N*N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k=0; k<N; k++){
                maxindexarray[i * N*N + j*N+k] = 0;
            }
        }
    }
    maxindexarray[max_index]=max_index;
    maxindexarray[max_index1]=max_index1;

     while (p<150){
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

        dft3d(frequency_max, image_max,N);
        for (int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                        for(int k=0; k<N;k++){
                        Fmax[i*N*N +j*N+k]=0;
                        }
                }
        }
        for(int i=0; i<count; i++){
            Fmax[samples[i]]=frequency_max[samples[i]];
        }

        for (int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                        for(int k=0; k<N;k++){
                        frequency_new[i*N*N +j*N+k]=frequency_new[i*N*N+j*N+k]-0.05*Fmax[i*N*N+j*N+k];
                        }
                }
        }
        idft3d(image_new,frequency_new, N);

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
        image_clean[max_index1]=max_val1;
        p=p+1;
        printf("The iteration number is %d \n", p);

    }


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int index = i*N*N + j*N + k;
                printf("%0.5f ", creal(image_clean[index]));
            }
        }
    }
    printf("\n");
    printf("%0.5f, %d", creal(max_val), max_index);
        printf("\n");
    printf("%0.5f, %d", creal(max_val1), max_index1);

    double complex sphereinput[N][N][N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k=0; k<N; k++){
                int index = i * N*N + j*N+k;
                sphereinput[i][j][k]=reconstructed_image[index];
            }
        }
    }
    //Code for opening up a text file called data and writing the new matrix to it that can be opened and plotted in Matlab.
    FILE *fp;

    fp = fopen("sphereinput40percent.txt", "w");

    // Write the array to the file
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                fprintf(fp, "%.5f ", creal(sphereinput[i][j][k]));
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
        // Close the file
    fclose(fp);

        double complex sphereresidual[N][N][N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k=0; k<N; k++){
                int index = i * N*N + j*N+k;
                sphereresidual[i][j][k]=image_new[index];
            }
        }
    }
    //Code for opening up a text file called data and writing the new matrix to it that can be opened and plotted in Matlab

    fp = fopen("sphereresidual.txt", "w");

    // Write the array to the file
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                fprintf(fp, "%.5f ", creal(sphereresidual[i][j][k]));
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
        // Close the file
    fclose(fp);

        double complex sphereCLEAN[N][N][N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k=0; k<N; k++){
                int index = i * N*N + j*N+k;
                sphereCLEAN[i][j][k]=image_clean[index];
            }
        }
    }
    //Code for opening up a text file called data and writing the new matrix to it that can be opened and plotted in Matlab.


    fp = fopen("3DCLEANTEST.txt", "w");

    // Write the array to the file
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                fprintf(fp, "%.5f ", creal(sphereCLEAN[i][j][k]));
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
        // Close the file
    fclose(fp);
}
