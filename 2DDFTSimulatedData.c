#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define PI 3.141592
// Creating the 2D DFT function
void dft2d(double complex *freq, double complex *f, int M, int N){
    #pragma omp parallel for
    for (int u = 0; u < M; u++)
    {
        for (int v = 0; v < N; v++)
        {
            double complex freq_point = 0;
            for (int x = 0; x < M; x++)
            {
                for (int y = 0; y < N; y++)
                {
                    freq_point += f[x * N + y] * cexp(-2 * PI * I * (((double)u * x / M)+((double)v * y/ N)));
                }
            }
            freq[u * N + v] = freq_point;
        }
    }
}

// function for finding the inverse DFT
void idft2d(double complex *image, double complex *freq, int M, int N){
    #pragma omp parallel for
    for (int x = 0; x < M; x++)
    {
        for (int y = 0; y < N; y++)
        {
            double complex image_points = 0;
            for (int u = 0; u < M; u++)
            {
                for (int v = 0; v < N; v++)
                {
                    image_points += freq[u * N + v] * cexp(2 * PI * I * (((double)u * x / M)+((double)v * y/ N)));
                }
            }
            image[x * N + y] = image_points/(M*N);
        }
    }
}

int main() {
    //Initialising variables
    int M=100;
    int N=100;
    double complex reconstructed_image[M * N];
    double complex random_freq_samples[M*N];
    double complex freq[M * N];

    //Code for opening a text file that has a 2D matrix from MATLAB of Simulated data
    int size;
    double complex *simulated_data;
    FILE *fp;
    fp = fopen("Image.txt", "r");
    int i = 0;
    double real;
    while (fscanf(fp, "%lf", &real) != EOF) {
        i++;
    }
    size = i;
    simulated_data = (double complex *) malloc(size * sizeof(double complex));
    rewind(fp);  // move the file pointer back to the beginning of the file
    i = 0;
    while (fscanf(fp, "%lf", &real) != EOF) {
        simulated_data[i] = real + 0.0 * I;
        i++;
    }
    fclose(fp);

    // Finding the DFT of the simulated data
    dft2d(freq, simulated_data, M, N);

    // Generating a random number of samples from the frequency domain
    int lower = 0, upper = M*N, count = ((M*N)/4);
    int output[count];

    srand(time(0));
    i = 0;
    int samples[count];
        for (int i=0; i<count; i++){
            samples[i]=0;
    }

    while(i < count) {
        int num = (rand() % (upper - lower + 1)) + lower;
        if(random_freq_samples[num] == 0) {
            random_freq_samples[num]=freq[num];
            samples[i]=num;
            i++;
        }
    }

    //Converting the random samples back into the frequency domain
    idft2d(reconstructed_image, random_freq_samples, M, N);

    // Calling the maximum value function
    double complex max_value;

    int maxindex=0;
    for (int i = 1; i < M*N; ++i)
    {
        if (creal(reconstructed_image[maxindex]) < creal(reconstructed_image[i]))
            maxindex = i;
    }


    //Initialising an array of zeros to input our maximum value into
    double complex image_max[M * N];
        for (int i=0; i<M*N; i++){
        image_max[i]=0;
    }

    //Placing the maximum value into the array of zeros at the right index
    image_max[maxindex]=reconstructed_image[maxindex];

    //Initialising an array of zeros for our CLEAN'ed points to get put into
    double complex image_clean[M * N];
        for (int i=0; i<M*N; i++){
        image_clean[i]=0;
    }

    //Placing our first value into the CLEAN image
    image_clean[maxindex]=reconstructed_image[maxindex];

    //Finding the DFT of our zeros array with the first max value in it
    double complex frequency_max[M * N];
    dft2d(frequency_max, image_max, M, N);

    //Creates an array of zeros and then inputs only the randomly sampled frequency points back into the array
    double complex Fmax[M*N];
    for (int i=0; i<M*N; i++){
        Fmax[i]=0;
    }
    int sample_pos;
    for(int i=0; i<count; i++){
        Fmax[samples[i]]=frequency_max[samples[i]];
    }

    // Creating a new array where we are taking the original sampled freq and subtracting the sampled position DFT of the maximum value
    double complex frequency_new[M*N];
    for(int i=0; i<M*N; i++){
        frequency_new[i]=random_freq_samples[i]-Fmax[i];
    }


    //Finding the iDFT of the new frequency to find the next maximum value
    double complex image_new[M*N];
    idft2d(image_new,frequency_new, M, N);

    // Calling the maximum value function
    double complex max_valuenew;
    int maxindex_new=0;

    for (int i = 1; i < M*N; ++i)
    {
        if (creal(image_new[maxindex_new]) < creal(image_new[i]))
            maxindex_new = i;
    }

    //Placing our second value into the CLEAN image
    image_clean[maxindex_new]=image_new[maxindex_new];
    int p=0;
  while (p<351){
  //while (creal(image_new[maxindex_new])>0.28*creal(reconstructed_image[maxindex])){

        //Initialising an array of zeros to input our maximum value into
        for (int i=0; i<M*N; i++){
            image_max[i]=0;
        }

        //Placing the maximum value into the array of zeros at the right index
        image_max[maxindex_new]=image_new[maxindex_new];
      
        // Same process as earlier
        dft2d(frequency_max, image_max, M, N);

        for (int i=0; i<M*N; i++){
            Fmax[i]=0;
        }
        int sample_pos;
        for(int i=0; i<count; i++){
            Fmax[samples[i]]=frequency_max[samples[i]];
        }

        for(int i=0; i<M*N; i++){
            frequency_new[i]=frequency_new[i]-Fmax[i];
        }
        idft2d(image_new,frequency_new, M, N);

        maxindex_new=0;

        for (int i = 1; i < M*N; ++i)
        {
            if (creal(image_new[maxindex_new]) < creal(image_new[i]))
            maxindex_new = i;
        }
        image_clean[maxindex_new]+=image_new[maxindex_new];
        p=p+1;
        printf("The iteration number is %d /n", p);
    }



    //Code for opening up a text file and writing the CLEAN image to it that can be plotted in Matlab.
    fp = fopen("CleanImage1.txt", "w");

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int index = i * N + j;
            fprintf(fp, "%g", creal(image_clean[index]));
            if (j != N - 1) {
                fprintf(fp, " ");
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    //Code for opening up a text file and writing the CL resultant to it that can be plotted in Matlab.
    fp = fopen("ResultantNoise1.txt", "w");

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int index = i * N + j;
            fprintf(fp, "%g", creal(image_new[index]));
            if (j != N - 1) {
                fprintf(fp, " ");
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    //Code for opening up a text file and writing the dirty image to it that can be plotted in Matlab.
    fp = fopen("ReconstructedImage1.txt", "w");

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int index = i * N + j;
            fprintf(fp, "%g", creal(reconstructed_image[index]));
            if (j != N - 1) {
                fprintf(fp, " ");
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);


    return 0;
}

