// ConsoleApplication2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

#include <thread>

#include <filesystem>

#include <math.h>     /* gcc mf.c -o wa -lm  */
#include <string.h>   /* gcc -O3 mf.c -o wa -lm  */


#pragma comment(lib, "libsndfile-1.lib")

#include "sndfile.h"  /* gcc -O3 mf.c -o wa -lm -lsndfile */ 
#include "const4mfsq.h"  /* gcc -O3 mfsq.c const4mfsq.h -o wa -lm -lsndfile  */
/* thresholdsq, fs, dt, N, N2, fp, TP, PBmin, TPsq, PBsqmin, first_ind, last_ind, HH
*/

#define PI      3.14159265358979323846    /* pi to machine precision, defined in math.h */
#define TWOPI   (2.0*PI)
# define numpts 30720000   /* ulimit -s unlimited */


int fileProc(const std::wstring& fileW, const std::wstring& pathDir,
             const bool bUsePBsqmi, const bool bScrPrint);
void four1(float* data, int nn, int isign);
void realft(float* data, int n, int isign);
void convolution(float* y, float* c);
void calc_PB(float* y, float& PB);
void calc_PBsq(float* y, float& PBsq);

int main()
{
    std::wstring pathDir;

    std::wifstream fsin;
    fsin.open("mf.ini");
    if (fsin.is_open())
    {
        std::getline(fsin, pathDir);
        fsin.close();
        std::wcout << std::endl << "Folder: " << pathDir.c_str() << std::endl << std::endl;
    }
    else
    {
        std::cout << "open file mf.ini failed.";
    }

    if (!std::filesystem::exists(pathDir))
    {
        std::wcout << std::endl << "Folder: " << pathDir.c_str() << " doesn't exist." << std::endl;
        return -1;
    }

    std::vector<std::wstring> files;
    for (const auto& entry : std::filesystem::directory_iterator(pathDir))
    {
        std::wstring fileW{ entry.path() };
        if (std::size_t found = fileW.find_last_of(L"/\\"); found != std::string::npos)
        {
            fileW = fileW.substr(found + 1);
            if (std::size_t found = fileW.find_last_of(L'.'); found != std::string::npos)
            {
                std::wstring ext = fileW.substr(found + 1);
                std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
                if (L"WAV" == ext)
                {
                    files.push_back(fileW);
                }
            }
        }
        //std::cout << entry.path() << std::endl;
    }

    std::wstring dataDir{L"DATA"};
    if (std::size_t found = pathDir.find_last_of(L"/\\"); found != std::string::npos)
    {
        dataDir = pathDir.substr(found + 1);
        pathDir = pathDir.substr(0, found + 1);
    }

    for (const auto& fileW : files)
    {
        std::wcout << fileW.c_str() << std::endl;
    }
    std::cout << std::endl << "------------------------------------------" << std::endl;

    printf("thresholdsq = %12.6e wav units\n", thresholdsq);
    printf("fs = %i samples per second\n", fs);
    printf("dt = %12.8e time step (seconds)\n", dt);
    printf("N  = %i sample length for matched filter \n", N);
    printf("N2 = %i sample overlap for matched filter\n", N2);
    printf("fp = %i frequency of porpoise click (Hz) \n", fp);

    printf("\n TP[65] = ");
    for (int m = 0; m < 65; m++)
    {
        printf("%12.8e, ", TP[m]);
    }
    printf("\n");
    std::cout << std::endl << "------------------------------------------" << std::endl;

    bool bUsePBsqmi{ false };
    if (bUsePBsqmi)
    {
        printf("PBsqmin = %12.8f \n", PBsqmin);
    }
    else
    {
        printf("PBmin = %12.8f \n", PBmin);
    }

    /*for (const auto& fileW : files)
    {
        fileProc(fileW, pathDir, bUsePBsqmi, true);
    }*/

    const int step{3};
    std::vector<std::thread> vthreads(step);
    for(int i = 0; i < files.size(); i += step)
    {
        for(int j = 0; j < step && i + j < files.size(); j++)
        {
            vthreads[j] = std::thread(fileProc, files[i + j], pathDir, bUsePBsqmi, false);
        }
        for (int j = 0; j < step && i + j < files.size(); j++)
        {
            if(vthreads[j].joinable())
                vthreads[j].join();
        }
    }

    return 0;
}

int fileProc(const std::wstring &fileW, const std::wstring& pathDir,
             const bool bUsePBsqmi, const bool bScrPrint)
{
    SNDFILE* sf{ 0 };
    SF_INFO info{ 0 };
    int num_channels{ 0 };
    int num{ 0 }, num_items{ 0 };
    int* buf{ 0 };
    int f{ 0 }, sr{ 0 }, c{ 0 };
    int nr{ 0 };

    int num_reads = numpts / 512 - 1;

    FILE* fpout{ nullptr };

    int   raw_data[1024]{ 0 };
    float y[1024]{ 0 };       /* data (normalized) */
    float yt[1024]{ 0 };      /* fourier transform */
    float ymean{ 0 };
    float cenvsq[1024]{ 0 };    /* envelope of the matched filter */
    float cmaxsq{ 0 }, cmax{ 0 }, cmaxdB{ 0 }, cms{ 0 }, PB{ 0 }, PBsq{ 0 };
    int k{ 0 };
    int end_number = 0;    /* last position that has been read in the wav file */

    if (bScrPrint)
    {
        /* printf("\n       HH[2048] = \n");
            printf("2*m      real       imaginary \n");
            for (m=0; m<1024; m++){
            printf(" %i  %12.8e  %12.8e  \n",2*m,HH[2*m],HH[2*m+1]);
            }*/

        std::cout << std::endl << "------------------------------------------" << std::endl;

        printf("\nSearch %ls for porpoise clicks.\n\n", fileW.c_str());
    }

    std::wstring fileN{ fileW };
    std::wstring directory_file{ pathDir + L"DATA/" + fileN };
    std::wstring fout{ fileN };
    if (std::size_t found = fileN.find_last_of('.'); found != std::string::npos)
    {
        fout = fileN.substr(0, found);
    }
    std::wstring directory_out{ pathDir + L"OUTPUT/" + fout + L".matches" };

    if (bScrPrint)
    {
        std::wcout << directory_file << std::endl;
        std::wcout << directory_out << std::endl << std::endl;
    }

    /* Open the WAV file. */
    info.format = 0;
    /*  sf = sf_open("file.wav",SFM_READ,&info); */
    /* SFM_READ  for read only mode, SFM_WRITE for write only, SFM_RDWR */
    std::string directory_fileS(directory_file.begin(), directory_file.end());
    sf = sf_open(directory_fileS.c_str(), SFM_READ, &info);
    if (sf == NULL)
    {
        printf("Failed to open the file.\n"); exit(-1);
    }
    /* Print some of the info, and figure out how much data to read. */
    f = info.frames;
    sr = info.samplerate;
    c = info.channels;
    if (bScrPrint)
    {
        printf("frames=%d\n", f);
        printf("samplerate=%d\n", sr);
        printf("channels=%d\n", c);
    }
    num_items = f * c;
    if (bScrPrint)
    {
        printf("num_items=%d\n", num_items);
    }

    /* Allocate space for the data to be read, then read it. */
    /*    buf = (int *) malloc(num_items*sizeof(int));
    num = sf_read_int(sf,buf,num_items); */
    num = sf_read_int(sf, raw_data, 512);
    end_number = 512;
    if (bScrPrint)
    {
        printf("------------------------------------------------- \n");
        printf("Read %d items, should equal %d\n", num, end_number);
        printf("raw_dat 0: %d \n", raw_data[0]);
        printf("raw_dat 1: %d \n", raw_data[1]);
        printf("raw_dat 510: %d \n", raw_data[510]);
        printf("raw_dat 511: %d \n", raw_data[511]);
        printf("------------------------------------------------- \n");
    }

    errno_t err;
    if ((err = _wfopen_s(&fpout, directory_out.c_str(), L"w")) == 0)
    {
        if (bScrPrint)
        {
            printf("\nOpened %ls in mode w\n", directory_out.c_str());
        }
    }
    else
    {
        printf("\nFailed to open %ls in mode w\n", directory_out.c_str());
        return -1;
    }

    if (nullptr != fpout)
    {
        fprintf(fpout, "%% mf.c run on %ls \n", directory_file.c_str());
        fprintf(fpout, "%% point-position  time-seconds  convolution  PB \n");
    }

    num_reads = num_items / 512;
    if (bScrPrint)
    {
        printf("num_reads: %d\n", num_reads);
    }

    for (nr = 1; nr <= num_reads; nr++) {
        num = sf_read_int(sf, raw_data + 512, 512);
        end_number = end_number + 512;
        /*      printf("Read another %d items, end_number %d\n",num,end_number); */
        for (int m = 0; m < N; m++) { y[m] = raw_data[m] / 2147483648.; }
        /*      for (m=0; m<256; m++) { y[m]=raw_data[m]/32768.; } */
        /*      printf("N = %d \n",N);
                printf("y[0]=%e\n",y[0]);
                printf("y[255]=%e\n",y[255]); */


        ymean = 0.;
        for (int m = 0; m < N; m++) { yt[m] = y[m]; ymean = ymean + y[m]; }
        ymean = ymean / N;
        for (int m = 0; m < N; m++) { yt[m] = yt[m] - ymean; }   /* extract the mean */
        realft(yt - 1, N2, 1);  /* yt-1 moves pointer 1 position to the left */
        convolution(yt, cenvsq);

        cmaxsq = cenvsq[256]; k = 256;
        for (int m = 256; m < N - 256; m++) {
            if (cenvsq[m] > cmaxsq) {
                cmaxsq = cenvsq[m];
                k = m;
            }
        }
        cms = 0.;
        /*      mm=0;
                for (m=0; m<k-34; m++){ cms=cms+cenvsq[m]; mm=mm+1; }
                for (m=k+33; m<N-100; m++){ cms=cms+cenvsq[m]; mm=mm+1; }
                printf("mm=%d should be 857\n",mm);                        */
        for (int m = 0; m < k - 34; m++) { cms = cms + cenvsq[m]; }
        for (int m = k + 33; m < N - 100; m++) { cms = cms + cenvsq[m]; }
        cms = cms / 857;                  /* mean square */

        if (bUsePBsqmi)
        {
            if (cmaxsq > 100 * cms && cmaxsq > thresholdsq) {
                calc_PBsq(y + k - 64, PBsq);
                if (PBsq > PBsqmin) {
                    int m = end_number - 1023 + k;
                    cmax = sqrt(cmaxsq);
                    cmaxdB = 20 * log10(3 * cmax) + 169; /* convert wav to dB */
                    fprintf(fpout, "%12d %15.8f %8.1f %8.2f \n", m, m * dt, cmaxdB, PBsq);
                    if (bScrPrint)
                    {
                        printf("%10d %12d %15.8f %8.1f %8.2f \n", nr, m, m * dt, cmaxdB, PBsq);
                    }
                }
            }
        }
        else
        {
            if (cmaxsq > 100 * cms && cmaxsq > thresholdsq) {
                calc_PB(y + k - 64, PB);
                if (PB > PBmin) {
                    int m = end_number - 1023 + k;
                    cmax = sqrt(cmaxsq);
                    cmaxdB = 20 * log10(3 * cmax) + 169;
                    fprintf(fpout, "%12d %15.8f %8.1f %8.2f \n", m, m * dt, cmaxdB, PB);
                    if (bScrPrint)
                    {
                        printf("%10d %12d %15.8f %8.1f %8.2f \n", nr, m, m * dt, cmaxdB, PB);
                    }
                }
            }
        }

        /*  shift raw data to accomodate the next 512 points */
        for (int m = 0; m < 512; m++) { raw_data[m] = raw_data[512 + m]; }
    }         /* end of nr loop that reads the data file */

    if (bScrPrint)
    {
        printf("end_number %d\n", end_number);
    }

    if (nullptr != fpout)
    {
        fclose(fpout);
    }

    printf("\nClosed %ls \n", directory_out.c_str());

    if (bScrPrint)
    {
        printf("------------------------------------------------- \n");
        printf("raw_dat 0: %d \n", raw_data[0]);
        printf("raw_dat 1: %d \n", raw_data[1]);
        printf("raw_dat 510: %d \n", raw_data[510]);
        printf("raw_dat 511: %d \n", raw_data[511]);
        printf("------------------------------------------------- \n");
    }
}

void calc_PB(float *y, float &PB)
{
    float yt[128], absyt, P, B;
    int m, mm, mp;

    for (m = 0; m < 128; m++) { yt[m] = y[m]; }
    realft(yt - 1, 64, 1);  /* data-1 moves pointer 1 position to the left */

    P = 0.; B = 0.;
    for (m = first_ind; m < last_ind + 1; m++) {
        mm = 2 * m;  mp = mm + 1;
        absyt = sqrt(pow(yt[mm], 2) + pow(yt[mp], 2));
        P = P + TP[m] * absyt;
        B = B + (1 - TP[m]) * absyt;
    }
    PB = P / B;
}

void calc_PBsq(float *y, float &PBsq)
{
    float yt[128], ytsq, P, B;
    int m, mm, mp;

    for (m = 0; m < 128; m++) { yt[m] = y[m]; }
    realft(yt - 1, 64, 1);  /* data-1 moves pointer 1 position to the left */

    P = 0.; B = 0.;
    for (m = first_ind; m < last_ind + 1; m++) {
        mm = 2 * m;  mp = mm + 1;
        ytsq = pow(yt[mm], 2) + pow(yt[mp], 2);
        P = P + TPsq[m] * ytsq;
        B = B + (1 - TPsq[m]) * ytsq;
    }
    PBsq = P / B;
}


void convolution(float *y, float *c)
{
    float yt[N], yr, yi;
    int m, mm, mp;

    /* Apply the filter -product of fourier transforms- (complex arithmetic) */
    for (m = 0; m < N2; m++) {
        mm = 2 * m; mp = mm + 1;
        yr = y[mm] * HH[mm] - y[mp] * HH[mp]; /* real part of the product */
        yi = y[mm] * HH[mp] + y[mp] * HH[mm]; /* imaginary part of the product */
        yt[mm] = yr;
        yt[mp] = yi;
    }
    /* printf("From convolution: HH times Fourier transform \n");
       m=200;
       printf("HH*yt[%d]=%e  %e\n",m,yt[2*m],yt[2*m+1]);

       printf("From convolution: filtered yt \n");
       for (m=0;m<N2;m++) {printf("yt[%d]=%e  %e\n",m,yt[2*m],yt[2*m+1]);}  */

       /* to reconstruct the signal with a pi/2 phase shift */
       /* Matlab code: X(2:N/2)=-i*X(2:N/2);  */
    c[0] = yt[0]; c[1] = yt[1];
    for (m = 1; m < N2; m++) {
        mm = 2 * m; mp = mm + 1;
        c[mm] = yt[mp];
        c[mp] = -yt[mm];
    }

    /* printf("From convolution: filtered c = yt with phase shift \n");
       for (m=0;m<N2;m++) {printf("c[%d]=%e  %e\n",m,c[2*m],c[2*m+1]);} */

    realft(yt - 1, N2, -1);   /* inverse transform */
    realft(c - 1, N2, -1);   /* inverse transform */

    for (m = 0; m < N; m++) { c[m] = c[m] / N2; yt[m] = yt[m] / N2; }

    /* calculate the square of the envelope */
    for (m = 0; m < N; m++) { c[m] = pow(c[m], 2) + pow(yt[m], 2); }

    /* printf("From convolution: envelope function is \n");     */
    /* for (m=0;m<N;m++) {printf("cenv[%d]=%e \n",m,c[m]);}     */
    /* printf("cnv=[ ");
       for (m=0;m<N;m++) {printf(" %8.3e ",c[m]);}
       printf("]; ");                                           */

}

/*-------------------------------------------------------------------------*
 * realft()
 *
 * from "numerical recipes in C".
 * Calculates the Fourier Transform of a set of 2*n real-valued data points.
 * Replaces this data (which is stored in the array data[1..2n]) by the
 * positive frequency half of its complex Fourier Transform. The real-valued
 * first and last components of the complex transform are returned as elements
 * data[1] and data[2] respectively. n must be a power of 2. This routine
 * also calculates the inverse transform of a complex data array if it is the
 * transform of real data. (Results in this case must be multiplied by 1/n.)
 *--------------------------------------------------------------------------*/

void realft(float *data, int n, int isign)
{
    int i, i1, i2, i3, i4, n2p3;
    float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    double wr, wi, wpr, wpi, wtemp, theta;

    theta = PI / (double)n;
    if (isign == 1) {      /* the forward transform here */
        c2 = -0.5;
        four1(data, n, 1);
    }
    else {                /* otherwise set up for the inverse transform */
        c2 = 0.5;
        theta = -theta;
    }
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    n2p3 = 2 * n + 3;
    /* case i=1 done separately below */
    for (i = 2; i <= n / 2; i++) {
        i4 = 1 + (i3 = n2p3 - (i2 = 1 + (i1 = i + i - 1)));
        /* the two separate transforms are separated out of data */
        h1r = c1 * (data[i1] + data[i3]);
        h1i = c1 * (data[i2] - data[i4]);
        h2r = -c2 * (data[i2] + data[i4]);
        h2i = c2 * (data[i1] - data[i3]);
        /* here they are recombined to form the true transform
        of the original real data */
        data[i1] = h1r + wr * h2r - wi * h2i;
        data[i2] = h1i + wr * h2i + wi * h2r;
        data[i3] = h1r - wr * h2r + wi * h2i;
        data[i4] = -h1i + wr * h2i + wi * h2r;
        /* the recurrence */
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    if (isign == 1) {
        /* squeeze the first and the last data together to get them
        all within the original data */
        data[1] = (h1r = data[1]) + data[2];
        data[2] = h1r - data[2];
    }
    else {
        /* this is the inverse transform for the case isign=-1 */
        data[1] = c1 * ((h1r = data[1]) + data[2]);
        data[2] = c1 * (h1r - data[2]);
        four1(data, n, -1);
    }
}


/*-------------------------------------------------------------------*
 * four1()
 *
 *  From "numerical recipes in C".
 *  Replace data by its DFT, if isign is input as 1; or replace data
 *  by nn times its inverse-DFT, if isign is input as -1.
 *  data is a complex array of length nn, input as a real
 *  array data[1...2nn]. nn must be an integer power of 2
 *-------------------------------------------------------------------*/

void four1(float *data, int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;
    /* this is the bit-reversal section of the routine */
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            /* exchange the two complex numbers */
            tempr = data[j];     data[j] = data[i];     data[i] = tempr;
            tempr = data[j + 1]; data[j + 1] = data[i + 1]; data[i + 1] = tempr;
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    /* here begins the Danielson-Lanczos section of the routine */
    /* Outer loop executed log2(nn) times */
    while (n > mmax) {
        istep = 2 * mmax;
        /* initialization for the trigonometric recurrence */
        theta = TWOPI / (isign * mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        /* here are the two nested loops */
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                /* this is Danielson-Lanczos formula */
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            /* trigonometric recurrence */
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
