
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "portaudio.h"


/*////////////////////////////////////////////////////////////*/

#define DTW_PASS_VALUE 4
#define VERY_BIG  (1e30)
#define MY_FRAMES (20000)
#define FRAMENUM 512
#define MYPI 3.14159265359
#define PERSONNUMBER 2

/*////////////////////////////////////////////////////////////*/



/* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
#define SAMPLE_RATE  (44100)
#define FRAMES_PER_BUFFER (512)
#define NUM_SECONDS     (3)
#define NUM_CHANNELS    (2)
/* #define DITHER_FLAG     (paDitherOff) */
#define DITHER_FLAG     (0) /**/
/** Set to 1 if you want to capture the recording to a file. */
#define WRITE_TO_FILE   (1)

/* Select sample format. */
#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif


/*////////////////////////////////////////////////////////////////////////////////////////////*/


float GetCoefficient(float* spectralData, unsigned int samplingRate, unsigned int NumFilters, unsigned int binSize, unsigned int m);
float NormalizationFactor(int NumFilters, int m);
float GetFilterParameter(unsigned int samplingRate, unsigned int binSize, unsigned int frequencyBand, unsigned int filterBand);
float GetMagnitudeFactor(unsigned int filterBand);
float GetCenterFrequency(unsigned int filterBand);



float dtw(char file1Name[], char file2Name[], char outputfileName[], int xsize, int ysize, int params){

double **globdist;
double **Dist;

double top, mid, bot, cheapest, total;
unsigned short int **move;
unsigned short int **warp;
unsigned short int **temp;

unsigned int I, X, Y, n, i, j, k;

unsigned int debug; /* debug flag */

float **x, **y; /*now 2 dimensional*/

FILE *file1, *file2, *glob, *debug_file, *output_file;

 /* open x-parameter file */

file1=fopen(file1Name, "rb");

/* open y-parameter file */

file2=fopen(file2Name,"rb");

/* allocate memory for x and y matrices */

x = malloc(xsize * sizeof(float *));

for (i=0; i < xsize; i++)
     x[i] = malloc(params * sizeof(float));

y = malloc(ysize * sizeof(float *));

for (i=0; i < ysize; i++)
     y[i] = malloc(params * sizeof(float));

/* allocate memory for Dist */

Dist = malloc(xsize * sizeof(double *));

for (i=0; i < xsize; i++)
	Dist[i] = malloc(ysize * sizeof(double));

     /* allocate memory for globdist */

globdist = malloc(xsize * sizeof(double *));

for (i=0; i < xsize; i++)
	globdist[i] = malloc(ysize * sizeof(double));

     /* allocate memory for move */

move = malloc(xsize * sizeof(short *));

for (i=0; i < xsize; i++)
	move[i] = malloc(ysize * sizeof(short));

     /* allocate memory for temp */

temp = malloc(xsize * 2 * sizeof(short *));

for (i=0; i < xsize*2; i++)
	temp[i] = malloc(2 * sizeof(short));

     /* allocate memory for warp */

warp = malloc(xsize * 2 * sizeof(short *));

for (i=0; i < xsize*2; i++)
	warp[i] = malloc(2 * sizeof(short));


/*read x parameter in x[]*/

for (i=0; i < xsize; i++)
  for (k=0; k < params; k++)
    fscanf(file1,"%f ",&x[i][k]);

/*read y parameter in y[]*/

for (i=0; i < ysize; i++)
  for (k=0; k < params; k++)
  	fscanf(file2,"%f ",&y[i][k]);


/*Compute distance matrix*/

for(i=0;i<xsize;i++) {
  for(j=0;j<ysize;j++) {
    total = 0;
    for (k=0;k<params;k++) {
      total = total + ((x[i][k] - y[j][k]) * (x[i][k] - y[j][k]));
    }
    Dist[i][j] = total;
  }
}


/*% for first frame, only possible match is at (0,0)*/

globdist[0][0] = Dist[0][0];
for (j=1; j<xsize; j++)
	globdist[j][0] = VERY_BIG;

globdist[0][1] = VERY_BIG;
globdist[1][1] = globdist[0][0] + Dist[1][1];
move[1][1] = 2;

for(j=2;j<xsize;j++)
	globdist[j][1] = VERY_BIG;

for(i=2;i<ysize;i++) {
	globdist[0][i] = VERY_BIG;
	globdist[1][i] = globdist[0][i-1] + Dist[1][i];

	for(j=2;j<xsize;j++) {
		top = globdist[j-1][i-2] + Dist[j][i-1] + Dist[j][i];
		mid = globdist[j-1][i-1] + Dist[j][i];
		bot = globdist[j-2][i-1] + Dist[j-1][i] + Dist[j][i];
		if( (top < mid) && (top < bot))
		{
		cheapest = top;
		I = 1;
		}
	else if (mid < bot)
		{
		cheapest = mid;
		I = 2;
		}
	else {cheapest = bot;
		I = 3;
		}

/*if all costs are equal, pick middle path*/
       if( ( top == mid) && (mid == bot))
	 I = 2;

	globdist[j][i] = cheapest;
	move[j][i] = I;
      }
}



X = ysize-1; Y = xsize-1; n=0;
warp[n][0] = X; warp[n][1] = Y;


while (X > 0 && Y > 0) {
n=n+1;


if (n>ysize *2) {fprintf (stderr,"Warning: warp matrix too large!");
exit(1);
}

if (move[Y] [X] == 1 )
	{
	warp[n][0] = X-1; warp[n][1] = Y;
	n=n+1;
	X=X-2; Y = Y-1;
	}
else if (move[Y] [X] == 2)
	{
	X=X-1; Y = Y-1;
	}
else if (move[Y] [X] == 3 )
	{
	warp[n] [0] = X;
	warp[n] [1] = Y-1; 
	n=n+1;
	X=X-1; Y = Y-2;
      }
else {fprintf(stderr,"Error: move not defined for X = %d Y = %d\n",X,Y); 
}
warp[n] [0] =X;
warp[n] [1] =Y;

}


/*flip warp*/
for (i=0;i<=n;i++) {
  temp[i][0] = warp[n-i][0];
  temp[i][1] = warp[n-i][1];

}

for (i=0;i<=n;i++) {
  warp[i][0] = temp[i][0];
  warp[i][1] = temp[i][1];

}


/* open output file */
output_file=fopen(outputfileName,"wb");

/*print warped trajectory to stdout*/
for (i=0;i<=n;i++)
     fprintf(output_file,"%d %d\n",warp[i][0]+1,warp[i][1]+1);

     fclose(output_file);

/* print global distance to globfile*/     

if ((glob=fopen("glob","w"))==NULL)
     fprintf(stderr,"Cannot open file glob\n");

fprintf(glob,"%f\n",globdist[xsize-1][ysize-1]);
fclose(glob);

	return globdist[xsize-1][ysize-1];

}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

void hamming(int len,float window[]){

    int i;

    float arg;

    arg = 2*MYPI / (len - 1);

    for (i = 0; i < len; i++){

        window[i] = window[i] *  (0.54 - 0.46 * cos(i * arg));

    }

}






void dft(float inreal[],int n) {

	int k;

	int t;

	float angle;

	float sumreal;

    	for (k = 0; k < n; k++) { 

         	sumreal = 0;

        	for (int t = 0; t < n; t++) {  
             
	    		angle = 2 * MYPI * t * k / n;

            		sumreal = sumreal+ inreal[t] * cos(angle);

        }

        inreal[k] = sumreal;

    }


}






void frameblocking(float pre[299],float  frame[299][1764]){
		
	int len,wid,s=0;
		
	for(len=0;len<299;len++){

		for(wid=0;wid<1764;wid++){

			frame[len][wid]=pre[s];

			s++;

		}

		s-=881;
	}


}







/*/////////////////////////////////////////////////////////////////////////////////////////*/




float GetCoefficient(float spectralData[], unsigned int samplingRate, unsigned int NumFilters, unsigned int binSize, unsigned int m)
{
  float result = 0.0f;
  float outerSum = 0.0f;
  float innerSum = 0.0f;
  unsigned int k, l;

  // 0 <= m < L
  if(m >= NumFilters)
  {
    // This represents an error condition - the specified coefficient is greater than or equal to the number of filters. The behavior in this case is undefined.
    return 0.0f;
  }

  result = NormalizationFactor(NumFilters, m);


  for(l = 1; l <= NumFilters; l++)
  {
    // Compute inner sum
    innerSum = 0.0f;
    for(k = 0; k < binSize - 1; k++)
    {
      innerSum += fabs(spectralData[k] * GetFilterParameter(samplingRate, binSize, k, l));

    }

    if(innerSum > 0.0f)
    {
      innerSum = log(innerSum); // The log of 0 is undefined, so don't use it
    }

    innerSum = innerSum * cos(((m * MYPI) / NumFilters) * (l - 0.5f));

    outerSum += innerSum;
  }

  result *= outerSum;

  return result;
}

/* 
 * Computes the Normalization Factor (Equation 6)
 * Used for internal computation only - not to be called directly
 */
float NormalizationFactor(int NumFilters, int m)
{
  float normalizationFactor = 0.0f;

  if(m == 0)
  {
    normalizationFactor = sqrt(1.0f / NumFilters);
  }
  else 
  {
    normalizationFactor = sqrt(2.0f / NumFilters);
  }
  
  return normalizationFactor;
}

/* 
 * Compute the filter parameter for the specified frequency and filter bands (Eq. 2)
 * Used for internal computation only - not the be called directly
 */
float GetFilterParameter(unsigned int samplingRate, unsigned int binSize, unsigned int frequencyBand, unsigned int filterBand)
{
  float filterParameter = 0.0f;

  float boundary = (frequencyBand * samplingRate) / binSize;    // k * Fs / N
  float prevCenterFrequency = GetCenterFrequency(filterBand - 1);    // fc(l - 1) etc.
  float thisCenterFrequency = GetCenterFrequency(filterBand);
  float nextCenterFrequency = GetCenterFrequency(filterBand + 1);

  if(boundary >= 0 && boundary < prevCenterFrequency)
  {
    filterParameter = 0.0f;
  }
  else if(boundary >= prevCenterFrequency && boundary < thisCenterFrequency)
  {
    filterParameter = (boundary - prevCenterFrequency) / (thisCenterFrequency - prevCenterFrequency);
    filterParameter *= GetMagnitudeFactor(filterBand);
  }
  else if(boundary >= thisCenterFrequency && boundary < nextCenterFrequency)
  {
    filterParameter = (boundary - nextCenterFrequency) / (thisCenterFrequency - nextCenterFrequency);
    filterParameter *= GetMagnitudeFactor(filterBand);
  }
  else if(boundary >= nextCenterFrequency && boundary < samplingRate)
  {
    filterParameter = 0.0f;
  }

  return filterParameter;
}

/* 
 * Compute the band-dependent magnitude factor for the given filter band (Eq. 3)
 * Used for internal computation only - not the be called directly
 */
float GetMagnitudeFactor(unsigned int filterBand)
{
  float magnitudeFactor = 0.0f;
  
  if(filterBand >= 1 && filterBand <= 14)
  {
    magnitudeFactor = 0.015;
  }
  else if(filterBand >= 15 && filterBand <= 48)
  {
    magnitudeFactor = 2.0f / (GetCenterFrequency(filterBand + 1) - GetCenterFrequency(filterBand -1));
  }

  return magnitudeFactor;
}

/*
 * Compute the center frequency (fc) of the specified filter band (l) (Eq. 4)
 * This where the mel-frequency scaling occurs. Filters are specified so that their
 * center frequencies are equally spaced on the mel scale
 * Used for internal computation only - not the be called directly
 */
float GetCenterFrequency(unsigned int filterBand)
{
  float centerFrequency = 0.0f;
  float exponent;

  if(filterBand == 0)
  {
    centerFrequency = 0;
  }
  else if(filterBand >= 1 && filterBand <= 14)
  {
    centerFrequency = (200.0f * filterBand) / 3.0f;
  }
  else
  {
    exponent = filterBand - 14.0f;
    centerFrequency = pow(1.0711703, exponent);
    centerFrequency *= 1073.4;
  }
  
  return centerFrequency;
}


/*************************************************************************************************************************/


typedef struct
{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
}
paTestData;

/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int recordCallback( const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData )
{
    paTestData *data = (paTestData*)userData;
    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) outputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;

    if( framesLeft < framesPerBuffer )
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }

    if( inputBuffer == NULL )
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = SAMPLE_SILENCE;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
        }
    }
    else
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}

/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int playCallback( const void *inputBuffer, void *outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo* timeInfo,
                         PaStreamCallbackFlags statusFlags,
                         void *userData )
{
    paTestData *data = (paTestData*)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE*)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) inputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;

    if( framesLeft < framesPerBuffer )
    {
        /* final buffer... */
        for( i=0; i<framesLeft; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
        for( ; i<framesPerBuffer; i++ )
        {
            *wptr++ = 0;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = 0;  /* right */
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else
    {
        for( i=0; i<framesPerBuffer; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}

/*******************************************************************/

int portaudio(int per)
{
    PaStreamParameters  inputParameters,
                        outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    char str[30];
    int                 i;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;
    SAMPLE              max, val;
    double              average;

    printf("patest_record.c\n"); fflush(stdout);

    data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; /* Record for a few seconds. */
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;
    numBytes = numSamples * sizeof(SAMPLE);
    data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
    if( data.recordedSamples == NULL )
    {
        printf("Could not allocate record array.\n");
        goto done;
    }
    for( i=0; i<numSamples; i++ ) data.recordedSamples[i] = 0;

    err = Pa_Initialize();
    if( err != paNoError ) goto done;

    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto done;
    }
    inputParameters.channelCount = 2;                    /* stereo input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    /* Record some audio. -------------------------------------------- */
    err = Pa_OpenStream(
              &stream,
              &inputParameters,
              NULL,                  /* &outputParameters, */
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              recordCallback,
              &data );
    if( err != paNoError ) goto done;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto done;
    printf("\n=== Now recording!! Please speak into the microphone. ===\n"); fflush(stdout);

    while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
    {
        Pa_Sleep(1000);
        printf("index = %d\n", data.frameIndex ); fflush(stdout);
    }
    if( err < 0 ) goto done;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto done;

    /* Measure maximum peak amplitude. */
    max = 0;
    average = 0.0;
    for( i=0; i<numSamples; i++ )
    {
        val = data.recordedSamples[i];
        if( val < 0 ) val = -val; /* ABS */
        if( val > max )
        {
            max = val;
        }
        average += val;
    }

    average = average / (double)numSamples;

    printf("sample max amplitude = "PRINTF_S_FORMAT"\n", max );
    printf("sample average = %lf\n", average );

    /* Write recorded data to a file. */
#if WRITE_TO_FILE
    {
        FILE  *fid;
        sprintf(str, "recorded%d.raw",per);
        fid = fopen(str, "wb");
        if( fid == NULL )
        {
            printf("Could not open file.");
        }
        else
        {
            fwrite( data.recordedSamples, NUM_CHANNELS * sizeof(SAMPLE), totalFrames, fid );
            fclose( fid );
            printf("Wrote data to 'recorded.raw'\n");
        }
    }
#endif



/*////////////////////////////////////////////////////////////////////////////////////////*/

FILE *fpointer ;
FILE *fpointer1 ;
FILE *fpointer2 ;
FILE *fpointer3 ;
FILE *fpointer4 ;
FILE *fpointer5 ;
FILE *fpointer6 ;
FILE *fpointer7 ;
FILE *fpointer8 ;


int l,len,wid,s=0,k=0,coeff;

double mfcc_result[13];

float pre_emph[264600];
float frame_block[299][1764];
float frame_temp[527436];


sprintf(str, "sample_pure %d.txt",per);
fpointer =fopen(str,"w+");
sprintf(str, "pre_emph %d.txt",per);
fpointer1 =fopen(str,"w+");
sprintf(str, "frame_blocking2d %d.txt",per);
fpointer2 =fopen(str,"w+");
sprintf(str, "frame_blocking1d %d.txt",per);
fpointer3 =fopen(str,"w+");
sprintf(str, "hamming2d %d.txt",per);
fpointer4 =fopen(str,"w+");
sprintf(str, "hamming1d %d.txt",per);
fpointer5 =fopen(str,"w+");
sprintf(str, "dft2d %d.txt",per);
fpointer6 =fopen(str,"w+");
sprintf(str, "dft1d %d.txt",per);
fpointer7 =fopen(str,"w+");
sprintf(str, "dtw%d.txt",per);
fpointer8 =fopen(str,"w+");




	printf("\n person %d pre_emphesis in progress\n",per);

	for(l = 0 ; l< numSamples ; l++){

		fprintf(fpointer ,"%d : %f\n",l,data.recordedSamples[l]);

		if(l!=0)
			pre_emph[l]=(data.recordedSamples[l]-(0.95*data.recordedSamples[l-1]));
		if(l==0)
			pre_emph[l]=(0.05*data.recordedSamples[l]);
	
		fprintf(fpointer1 ,"%d : %f\n",l,pre_emph[l]);	


	}

	printf(" person %d samples and pre_emph saved to file\n",per);
	
	printf("\n person %d frameblocking in progress\n",per);
	
	frameblocking(pre_emph,frame_block);
	

	for(len=0;len<299;len++){

		for(wid=0;wid<1764;wid++){

			fprintf(fpointer2 , "%d : %f\n",k,frame_block[len][wid]);

			frame_temp[k]=frame_block[len][wid];

			k++;

			s++;

		}

		s-=881;
	}

	printf(" person %d frameblock_2d saved to file\n",per);
	
	
	
	for(k=0;k<527436;k++){


		fprintf(fpointer3 , "%d : %f\n",k,frame_temp[k]);


	}

	printf(" person %d frameblock_1d saved to file\n",per);


	printf("\n person %d hamming_window in progress\n",per);

	for(k=0;k<299;k++){
	
		hamming(1764,frame_block[k]);
	
	}	


	printf(" person %d hamming done\n",per);




	k=0;

	for(len=0;len<299;len++){

		for(wid=0;wid<1764;wid++){

			fprintf(fpointer4 , "%d : %f\n",s,frame_block[len][wid]);

			frame_temp[k]=frame_block[len][wid];

			k++;

		}

	}

		
	printf(" person %d hamming_2d saved to file\n",per);


	
	for(k=0;k<527436;k++){


		fprintf(fpointer5 , "%d : %f\n",k,frame_temp[k]);


	}
	
	printf(" person %d hamming_1d saved to file\n",per);

	printf("\n person %d dft in progress\n",per);

	for(k=0;k<299;k++){

		
	
		dft(frame_block[k],1764);
	
	}

	
	printf(" person %d dft done\n",per);


	k=0;

	for(len=0;len<299;len++){

		for(wid=0;wid<1764;wid++){

			fprintf(fpointer6 , "%d : %f\n",k,frame_block[len][wid]);

			frame_temp[k]=frame_block[len][wid];

			k++;

		}

	}

	printf(" person %d dft_2d saved to file\n",per);



	for(k=0;k<527436;k++){


		fprintf(fpointer7 , "%d : %f\n",k,frame_temp[k]);


	}

	printf(" person %d dft_1d saved to file\n",per);


	printf("\n person %d getcoefficient in progress\n",per);


	for(coeff = 0; coeff < 13; coeff++){

			mfcc_result[coeff] = GetCoefficient(frame_temp, 44100, 20, 128, coeff);

	}

	printf(" person %d getcoefficient done\n",per);
	
		
	for(coeff = 1; coeff < 13; coeff++){
	fprintf(fpointer8 , "%f\n",mfcc_result[coeff]);
	printf("%lf\n",mfcc_result[coeff]);
	}
	

fclose(fpointer);
fclose(fpointer1);
fclose(fpointer2);
fclose(fpointer3);
fclose(fpointer4);
fclose(fpointer5);
fclose(fpointer6);
fclose(fpointer7);
fclose(fpointer8);


/*////////////////////////////////////////////////////////////////////////////////////////*/




    /* Playback recorded data.  -------------------------------------------- */
    data.frameIndex = 0;

    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto done;
    }
    outputParameters.channelCount = 2;                     /* stereo output */
    outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    printf("\n=== Now playing back. ===\n"); fflush(stdout);
    err = Pa_OpenStream(
              &stream,
              NULL, /* no input */
              &outputParameters,
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              playCallback,
              &data );
    if( err != paNoError ) goto done;

    if( stream )
    {
        err = Pa_StartStream( stream );
        if( err != paNoError ) goto done;
        
        printf("Waiting for playback to finish.\n"); fflush(stdout);

        while( ( err = Pa_IsStreamActive( stream ) ) == 1 ) Pa_Sleep(100);
        if( err < 0 ) goto done;
        
        err = Pa_CloseStream( stream );
        if( err != paNoError ) goto done;
        
        printf("Done.\n"); fflush(stdout);
    }

done:
    Pa_Terminate();
    if( data.recordedSamples )       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    if( err != paNoError )
    {
        fprintf( stderr, "An error occured while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          /* Always return 0 or 1, but no other return codes. */
    }
    return err;
}







int main(){

	int i;
	
	int personnumber=2;

	for(i=0;i<personnumber;i++){

		printf("Press Enter:");
		
		getchar();

		portaudio(i);

	}


	float dtwanswer = dtw("dtw0.txt","dtw1.txt", "output.txt", 12, 12, 1);
	
	printf("%f\n", dtwanswer);
	
	if(dtwanswer < DTW_PASS_VALUE)
		puts("PASS");
	else
		puts("FAIL");

}

