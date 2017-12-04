/*
 * DSP Example is part of the DSP Lecture at TEC-Costa Rica
 * Copyright (C) 2017  Pablo Alvarado
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file   controlVolume.cpp
 *         Implements 10 bands eq with IIR filters, fourth order
 * \authors Rolando Coto, Daniel Leon, Daniel Sandoval
 * \date   2017.19.06
 *
 * $Id: controlVolume.cpp $
 */

#include "controlvolume.h"
#include <cmath>
#include <cstring>
#include <fftw3.h>

// se agregaron para la DFT
#include <iostream>
#include <fstream>




#undef _DSP_DEBUG
#define _DSP_DEBUG

#ifdef _DSP_DEBUG
#define _debug(x) std::cerr << x
#include <iostream>
#else
#define _debug(x)
#endif




#define REAL 0
#define IMAG 1

void fft(double ent[2048][2], double sal[2048][2]){
    fftw_plan planfft = fftw_plan_dft_1d(2048, ent, sal, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(planfft);
    fftw_destroy_plan(planfft);
    fftw_cleanup();
}

void idft(double ent[2048][2],double sal[2048][2]){
    fftw_plan planidft = fftw_plan_dft_1d(2048, ent ,sal, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(planidft);
    fftw_destroy_plan(planidft);
    fftw_cleanup();

    for (int i=0; i<2048; i++){ //Normalización
        sal[i][REAL] /= 2048;
        sal[i][IMAG] /= 2048;
    }
}


void hn(double K,double a1,double a2,double a3,double a4,double a5,double a6,double b1,double b2,double b3,double b4,double b5,double b6, double H[2048][2]){


    float y_1=0, y_2=0, y_3=0, y_4=0, y_5=0, y_6=0; //Init Condiciones iniciales

    for(int i=0; i<2048; i++){
        H[i][IMAG] = 0;
        if(i<300){
            switch(i){
            case 0: //Caso n=0. Solo la entrada no retrasada
                H[i][REAL]= K;
                break;
            case 1:
                H[i][REAL]= K*b1 - a1*y_1; //Caso n=1. Solo entrada retrada con k=1 y salidas enteriores corrrespondietes
                break;
            case 2:
                H[i][REAL]= K*b2 - a1*y_1 - a2*y_2; //Caso n=2, x(n-k) k=2, salidas correspondientes
                break;
            case 3:
                H[i][REAL]= K*b3 - a1*y_1 - a2*y_2 - a3*y_3; //Caso n=3, k=3.
                break;
            case 4:
                H[i][REAL]= K*b4 - a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4; //n = k = 4
                break;
            case 5:
                H[i][REAL]= K*b5 - a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4 - a5*y_5; // Idem
                break;
            case 6:
                H[i][REAL]= K*b6 - a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4 - a5*y_5 - a6*y_6; //Idem
                break;
            default:
                H[i][REAL]= -a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4 - a5*y_5 - a6*y_6; // n > 6
                break;
            }
            y_6 = y_5; //Reasignación de salidas anteriores
            y_5 = y_4;
            y_4 = y_3;
            y_3 = y_2;
            y_2 = y_1;
            y_1 = H[i][REAL];
        }
        else{
            H[i][REAL]=0;
        }
    }
}

/**
 * Constructor
 */
controlVolume::controlVolume(){
    // se inicializan variables para asignar condiciones iniciales
    a=0;b=0;c=0;d=0;
    e=0;f=0;g=0;h=0;
    i=0,j=0,k=0,l=0;
    m=0; q=0;o=0;p=0;
    r=0;s=0; t=0;u=0;
    aa=0;ab=0;ac=0;ad=0;
    ba=0;bb=0;bc=0;bd=0;
    ca=0;cb=0;cc=0;cd=0;
    da=0;db=0;dc=0;dd=0;
    ea=0;eb=0;ec=0;ed=0;



    // variables de energia
    energia31=0;   // variable donde se guarda la cantidad de energia de la senal
    energia64=0;energia125=0;energia250=0;energia500=0;energia1000=0;energia2000=0;
    energia4000=0;energia8000=0;energia16000=0;

    // variables de temporizacion
    y1=0;y2=0;y3=0;y4=0;y5=0;     //Variable para circuito RC
    y6=0;y7=0;y8=0;y9=0;y10=0;

}


/*
 * Destructor
 */

controlVolume::~controlVolume(){



}


/**
 * Filtro de 16 kHz de orden 2
 * Parametros de entrada:
 * tama;o de bloque
 * ganancia
 * condicion de primer frame
 * puntero entrada
 * puntero salida
 */
//-------------------------------------------------FILTRO DE 16KHz de orden 2----------------------------------------------------------------///
void controlVolume::filter(int blockSize, int volumeGain, bool inicial, float *in, float *out){

        int N = 2048;

        fftw_complex *x_n;
        x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *X;
        X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *y;
        y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *Y;
        Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *h10;
        h10 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *H10;
        H10 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        double *out_10 = new double[blockSize];
        double *in_10  = new double[blockSize];

        if (haux){
            sol[300] = {0};
            haux=false;
        }

    /*

        double K = 0.0322030478247811;
        double b1 = 4.08388652813182;
        double b2 = 8.10881739310181;
        double b3 = 10.6835815861358;
        double b4 = 10.3961708909533;
        double b5 = 7.5146551670359;
        double b6 = 3.88925305960443;
        double b7 = 1.32749953645413;
        double b8 = 0.234132025906004;
        double a1 = 0.258281798956747;
        double a2 = -3.45362835066451;
        double a3 = -0.221624362721121;
        double a4 = 4.98062979981792;
        double a5 = -0.221624362721121;
        double a6 = -3.45362835066451;
        double a7 = 0.258281798956747;
        double a8 = 1;
        hn8(K,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,h10);
*/

        double K=0.07875912;
        double a1=3.44387319;
        double a2=5.29844359;
        double a3=5.13528975;
        double a4=3.47957963;
        double a5=1.47740514;
        double a6=0.28411879;
        double b1=0.04060383;
        double b2=-2.91809914;
        double b3=0;
        double b4=2.91809914;
        double b5=-0.04060383;
        double b6=-1;


        hn(K,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,h10);
        fft(h10,H10); //calcula H[k]

        for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
            if (i<300){
                    x_n[i][REAL]= sol[i];
                    x_n[i][IMAG]= 0;
                }

            else {
                if(i<1324){
                        in_10[i-300]=static_cast<double>(in[i-300]);
                        x_n[i][REAL] = in_10[i-300];
                        x_n[i][IMAG] = 0;
                   }

                else {
                    x_n[i][REAL] = 0;
                    x_n[i][IMAG] = 0;
                }
            }
        }

        for(int i=0;i<300;i++){  // rellena con la entrada actual
            sol[i] = in_10[724 + i];
        }

        fft(x_n,X); //calcula X[k]

        for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
            Y[i][REAL]=H10[i][REAL]*X[i][REAL]-H10[i][IMAG]*X[i][IMAG];
            Y[i][IMAG]=H10[i][REAL]*X[i][IMAG]+H10[i][IMAG]*X[i][REAL];
        }

        idft(Y,y);        //calcula y[n]

        for(int n=0;n<1024;++n){
            out_10[n]= y[n + 299][REAL];
            out_10[n]=(0.02)*(volumeGain)*(out_10[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
            out[n]=static_cast<float>(out_10[n]);// se hace conversion de double a float
        }

       // energia16000=FFT(blockSize,out_10);//se determina la energia de la banda

        delete out_10;
        delete in_10;
        fftw_free(x_n);
        fftw_free(X);
        fftw_free(Y);
        fftw_free(h10);
        fftw_free(H10);

        /*
        for(int i =0; i<blockSize; i++){
            out[i] = 0;
        }
        */
}


 //-------------------------------------------------FILTRO DE 8KHz------------------------------------------------------------------------------//

void controlVolume::filter_8k(int blockSize, int volumeGain, bool inicial, float *in, float *out){

     int N = 2048;

     fftw_complex *x_n;
     x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

     fftw_complex *X;
     X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

     fftw_complex *y;
     y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

     fftw_complex *Y;
     Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

     fftw_complex *h8;
     h8 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

     fftw_complex *H8;
     H8 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

     double *out_8 = new double[blockSize];
     double *in_8  = new double[blockSize];

     if (haux8){
         sol8[14] = {0};
         haux8=false;
     }

     for(int i=0;i<2048;i++){
         switch (i) {
         case 0:
             h8[i][REAL] =  -0.004131972221000227016;
             h8[i][IMAG] =  0;
             break;
         case 1:
             h8[i][REAL] =  0.018860988212350455;
             h8[i][IMAG] =  0;
             break;
         case 2:
             h8[i][REAL] = 0.07724184701232166;
             h8[i][IMAG] =  0;
             break;
         case 3:
             h8[i][REAL] =  0.043525981075392392;
             h8[i][IMAG] =  0;
             break;
         case 4:
             h8[i][REAL] =  -0.10368955092873815;
             h8[i][IMAG] =  0;
             break;
         case 5:
             h8[i][REAL] =  -0.17462974943888687;
             h8[i][IMAG] =  0;
             break;
         case 6:
             h8[i][REAL] = -0.026317712854047076;
             h8[i][IMAG] =  0;
             break;
         case 7:
             h8[i][REAL] =  0.18362657716751676;
             h8[i][IMAG] =  0;
             break;
         case 8:
             h8[i][REAL] =  0.18362657716751676;
             h8[i][IMAG] =  0;
             break;
         case 9:
             h8[i][REAL] =  -0.026317712854047076;
             h8[i][IMAG] =  0;
             break;
         case 10:
             h8[i][REAL] =  -0.17462974943888687;
             h8[i][IMAG] =  0;
             break;

         case 11:
             h8[i][REAL] =  -0.10368955092873815;
             h8[i][IMAG] =  0;
             break;
         case 12:
             h8[i][REAL] =  0.043525981075392392;
             h8[i][IMAG] =  0;
             break;
         case 13:
             h8[i][REAL] =  0.07724184701232166;
             h8[i][IMAG] =  0;
             break;
         case 14:
             h8[i][REAL] =  0.018860988212350455;
             h8[i][IMAG] =  0;
             break;
         case 15:
             h8[i][REAL] =  -0.004131972221000227016;
             h8[i][IMAG] =  0;
             break;

         default:
             h8[i][REAL] =  0;
             h8[i][IMAG] =  0;
             break;
         }
     }

     fft(h8,H8); //calcula H[k]

     for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
         if (i<14){
                 x_n[i][REAL]= sol8[i];
                 x_n[i][IMAG]= 0;
             }

         else {
             if(i<(1024 +14)){
                     in_8[i-14]=static_cast<double>(in[i-14]);
                     x_n[i][REAL] = in_8[i-14];
                     x_n[i][IMAG] = 0;
                }

             else {
                 x_n[i][REAL] = 0;
                 x_n[i][IMAG] = 0;
             }
         }
     }

     for(int i=0;i < 14;i++){  // rellena con la entrada actual
         sol8[i] = in_8[(1024-14) + i];
     }

     fft(x_n,X); //calcula X[k]

     for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
         Y[i][REAL]=H8[i][REAL]*X[i][REAL]-H8[i][IMAG]*X[i][IMAG];
         Y[i][IMAG]=H8[i][REAL]*X[i][IMAG]+H8[i][IMAG]*X[i][REAL];
     }

     idft(Y,y);        //calcula y[n]

     for(int n=0;n<1024;++n){
         out_8[n]= y[n + 14][REAL];
         out_8[n]=(0.02)*(volumeGain)*(out_8[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
         out[n]=static_cast<float>(out_8[n]);// se hace conversion de double a float
     }

    // energia8000=FFT(blockSize,out_8);//se determina la energia de la banda

     delete out_8;
     delete in_8;
     fftw_free(x_n);
     fftw_free(X);
     fftw_free(Y);
     fftw_free(h8);
     fftw_free(H8);

     /*
     for(int i =0; i<blockSize; i++){
         out[i] = 0;
     }
     */

 }


//-------------------------------------------------FILTRO DE 4KHz------------------------------------------------------------------------------//
void controlVolume::filter_4k(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 2kHz
    int N = 2048;

    fftw_complex *x_n;
    x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *h4;
    h4 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *H4;
    H4 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    double *out_4 = new double[blockSize];
    double *in_4  = new double[blockSize];

    if (haux4){
        sol4[30] = {0};
        haux4=false;
    }

    for(int i=0;i<2048;i++){
        switch (i) {
        case 0:
            h4[i][REAL] =  -0.004006833477669;
            h4[i][IMAG] =  0;
            break;
        case 1:
            h4[i][REAL] =  -0.006966907948778;
            h4[i][IMAG] =  0;
            break;
        case 2:
            h4[i][REAL] = -0.000081447279846;
            h4[i][IMAG] =  0;
            break;
        case 3:
            h4[i][REAL] =  0.016899788357053;
            h4[i][IMAG] =  0;
            break;
        case 4:
            h4[i][REAL] =   0.036758700596718;
            h4[i][IMAG] =  0;
            break;
        case 5:
            h4[i][REAL] =   0.047570718829542;
            h4[i][IMAG] =  0;
            break;
        case 6:
            h4[i][REAL] =  0.038900469575622;
            h4[i][IMAG] =  0;
            break;
        case 7:
            h4[i][REAL] =  0.008325708765166;
            h4[i][IMAG] =  0;
            break;
        case 8:
            h4[i][REAL] =  -0.035398535410993;
            h4[i][IMAG] =  0;
            break;
        case 9:
            h4[i][REAL] =  -0.074781635551810;
            h4[i][IMAG] =  0;
            break;
        case 10:
            h4[i][REAL] =  -0.091057430123658;
            h4[i][IMAG] =  0;
            break;

        case 11:
            h4[i][REAL] =  -0.073304794560737;
            h4[i][IMAG] =  0;
            break;
        case 12:
            h4[i][REAL] =  -0.024679936892059;
            h4[i][IMAG] =  0;
            break;
        case 13:
            h4[i][REAL] =  0.037653168801664;
            h4[i][IMAG] =  0;
            break;
        case 14:
            h4[i][REAL] =  0.089315840237288;
            h4[i][IMAG] =  0;
            break;
        case 15:
            h4[i][REAL] =  0.109287718669758;
            h4[i][IMAG] =  0;
            break;
        case 16:
                h4[i][REAL] =  0.089315840237288;
                h4[i][IMAG] =  0;
                break;
        case 17:
                h4[i][REAL] =  0.037653168801664;
                h4[i][IMAG] =  0;
                break;
        case 18:
                h4[i][REAL] =  -0.024679936892059;
                h4[i][IMAG] =  0;
                break;
        case 19:
                h4[i][REAL] =  -0.073304794560737;
                h4[i][IMAG] =  0;
                break;
        case 20:
                h4[i][REAL] =  -0.091057430123658;
                h4[i][IMAG] =  0;
                break;
        case 21:
                h4[i][REAL] =  -0.074781635551810;
                h4[i][IMAG] =  0;
                break;
        case 22:
                h4[i][REAL] =  -0.035398535410993;
                h4[i][IMAG] =  0;
                break;
        case 23:
                h4[i][REAL] =  0.008325708765166;
                h4[i][IMAG] =  0;
                break;
        case 24:
                h4[i][REAL] =  0.038900469575622;
                h4[i][IMAG] =  0;
                break;
        case 25:
                h4[i][REAL] =  0.047570718829542;
                h4[i][IMAG] =  0;
                  break;
        case 26:
                h4[i][REAL] =  0.036758700596718;
                h4[i][IMAG] =  0;
                break;
        case 27:
                h4[i][REAL] =  0.016899788357053;
                h4[i][IMAG] =  0;
                break;
        case 28:
                h4[i][REAL] =  -0.000081447279846;
                h4[i][IMAG] =  0;
                break;
        case 29:
                h4[i][REAL] =  -0.006966907948778;
                h4[i][IMAG] =  0;
                break;
        case 30:
                h4[i][REAL] =  -0.004006833477669;
                h4[i][IMAG] =  0;
                break;
        default:
            h4[i][REAL] =  0;
            h4[i][IMAG] =  0;
            break;
        }
    }

    fft(h4,H4); //calcula H[k]

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<29){
                x_n[i][REAL]= sol4[i];
                x_n[i][IMAG]= 0;
            }

        else {
            if(i<(1024 +29)){
                    in_4[i-29]=static_cast<double>(in[i-29]);
                    x_n[i][REAL] = in_4[i-29];
                    x_n[i][IMAG] = 0;
               }

            else {
                x_n[i][REAL] = 0;
                x_n[i][IMAG] = 0;
            }
        }
    }

    for(int i=0;i < 29;i++){  // rellena con la entrada actual
        sol4[i] = in_4[(1024-29) + i];
    }

    fft(x_n,X); //calcula X[k]

    for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
        Y[i][REAL]=H4[i][REAL]*X[i][REAL]-H4[i][IMAG]*X[i][IMAG];
        Y[i][IMAG]=H4[i][REAL]*X[i][IMAG]+H4[i][IMAG]*X[i][REAL];
    }

    idft(Y,y);        //calcula y[n]

    for(int n=0;n<1024;++n){
        out_4[n]= y[n + 29][REAL];
        out_4[n]=(0.02)*(volumeGain)*(out_4[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
        out[n]=static_cast<float>(out_4[n]);// se hace conversion de double a float
    }

    //energia4000=FFT(blockSize,out_4);//se determina la energia de la banda

    delete out_4;
    delete in_4;
    fftw_free(x_n);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h4);
    fftw_free(H4);

    /*
    for(int i =0; i<blockSize; i++){
        out[i] = 0;
    }
    */

}


//-------------------------------------------------FILTRO DE 2KHz------------------------------------------------------------------------------//
void controlVolume::filter_2k(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 2kHz
    int N = 2048;

    fftw_complex *x_n;
    x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *h2;
    h2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *H2;
    H2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    double *out_2 = new double[blockSize];
    double *in_2  = new double[blockSize];

    if (haux2){
        sol2[35] = {0};
        haux2=false;
    }

    for(int i=0;i<2048;i++){
        switch (i) {
        case 0:
            h2[i][REAL] =   0.021451373614944;
            h2[i][IMAG] =   0;
            break;
        case 1:
            h2[i][REAL] =  0.011457401363579;
            h2[i][IMAG] =   0;
            break;
        case 2:
            h2[i][REAL] = -0.001287396491463;
            h2[i][IMAG] =  0;
            break;
        case 3:
            h2[i][REAL] =  -0.015624607935512;
            h2[i][IMAG] =   0;
            break;
        case 4:
            h2[i][REAL] =  -0.030078670032301;
            h2[i][IMAG] =   0;
            break;
        case 5:
            h2[i][REAL] =  -0.043019296183234;
            h2[i][IMAG] =  0;
            break;
        case 6:
            h2[i][REAL] =  -0.052854689840613;
            h2[i][IMAG] =  0;
            break;
        case 7:
            h2[i][REAL] =  -0.058231625868760;
            h2[i][IMAG] =  0;
            break;
        case 8:
            h2[i][REAL] =  -0.058216547443086;
            h2[i][IMAG] =  0;
            break;
        case 9:
            h2[i][REAL] =  -0.052433322472396;
            h2[i][IMAG] =  0;
            break;
        case 10:
            h2[i][REAL] =  -0.041138122202598;
            h2[i][IMAG] =  0;
            break;
        case 11:
            h2[i][REAL] =  -0.025219458682030;
            h2[i][IMAG] =  0;
            break;
        case 12:
            h2[i][REAL] =  -0.006120789280146;
            h2[i][IMAG] =  0;
            break;
        case 13:
            h2[i][REAL] =  0.014306957414070;
            h2[i][IMAG] =  0;
            break;
        case 14:
            h2[i][REAL] =  0.034006382569578;
            h2[i][IMAG] =  0;
            break;
        case 15:
            h2[i][REAL] =  0.050944387549309;
            h2[i][IMAG] =  0;
            break;
        case 16:
            h2[i][REAL] =   0.063344733346520;
            h2[i][IMAG] =   0;
            break;
        case 17:
            h2[i][REAL] =   0.069893768198607;
            h2[i][IMAG] =   0;
            break;
        case 18:
            h2[i][REAL] =   0.069893768198607;
            h2[i][IMAG] =   0;
            break;
        case 19:
            h2[i][REAL] =   0.063344733346520;
            h2[i][IMAG] =   0;
            break;
        case 20:
            h2[i][REAL] =   0.050944387549309;
            h2[i][IMAG] =   0;
            break;
        case 21:
            h2[i][REAL] =   0.034006382569578;
            h2[i][IMAG] =   0;
            break;
        case 22:
            h2[i][REAL] =   0.014306957414070;
            h2[i][IMAG] =   0;
            break;
        case 23:
            h2[i][REAL] =   -0.006120789280146;
            h2[i][IMAG] =   0;
            break;
        case 24:
            h2[i][REAL] =   -0.025219458682030;
            h2[i][IMAG] =   0;
            break;
        case 25:
            h2[i][REAL] =   -0.041138122202598;
            h2[i][IMAG] =   0;
            break;
        case 26:
            h2[i][REAL] =   -0.052433322472396;
            h2[i][IMAG] =   0;
            break;
        case 27:
            h2[i][REAL] =   -0.058216547443086;
            h2[i][IMAG] =   0;
            break;
        case 28:
            h2[i][REAL] =   -0.058231625868760;
            h2[i][IMAG] =   0;
            break;
        case 29:
            h2[i][REAL] =   -0.052854689840613;
            h2[i][IMAG] =   0;
            break;
        case 30:
            h2[i][REAL] =   -0.043019296183234;
            h2[i][IMAG] =   0;
            break;
        case 31:
            h2[i][REAL] =   -0.030078670032301;
            h2[i][IMAG] =   0;
            break;
        case 32:
            h2[i][REAL] =   -0.015624607935512;
            h2[i][IMAG] =   0;
            break;
        case 33:
            h2[i][REAL] =   -0.001287396491463;
            h2[i][IMAG] =   0;
            break;
        case 34:
            h2[i][REAL] =   0.011457401363579;
            h2[i][IMAG] =   0;
            break;
        case 35:
            h2[i][REAL] =   0.021451373614944;
            h2[i][IMAG] =   0;
            break;

        default:
            h2[i][REAL] =  0;
            h2[i][IMAG] =  0;
            break;
        }
    }

    fft(h2,H2); //calcula H[k]

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<34){
                x_n[i][REAL]= sol2[i];
                x_n[i][IMAG]= 0;
            }

        else {
            if(i<(1024 +34)){
                    in_2[i-34]=static_cast<double>(in[i-34]);
                    x_n[i][REAL] = in_2[i-34];
                    x_n[i][IMAG] = 0;
               }

            else {
                x_n[i][REAL] = 0;
                x_n[i][IMAG] = 0;
            }
        }
    }

    for(int i=0;i < 34;i++){  // rellena con la entrada actual
        sol2[i] = in_2[(1024-34) + i];
    }

    fft(x_n,X); //calcula X[k]

    for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
        Y[i][REAL]=H2[i][REAL]*X[i][REAL]-H2[i][IMAG]*X[i][IMAG];
        Y[i][IMAG]=H2[i][REAL]*X[i][IMAG]+H2[i][IMAG]*X[i][REAL];
    }

    idft(Y,y);        //calcula y[n]

    for(int n=0;n<1024;++n){
        out_2[n]= y[n + 34][REAL];
        out_2[n]=(0.02)*(volumeGain)*(out_2[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
        out[n]=static_cast<float>(out_2[n]);// se hace conversion de double a float
    }

    // energia2000=FFT(blockSize,out_2);//se determina la energia de la banda

    delete out_2;
    delete in_2;
    fftw_free(x_n);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h2);
    fftw_free(H2);

    /*
    for(int i =0; i<blockSize; i++){
        out[i] = 0;
    }
    */
}

    //-------------------------------------------------FILTRO DE 1KHz------------------------------------------------------------------------------//
void controlVolume::filter_1k(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 1kHz

        int N = 2048;

        fftw_complex *x_n;
        x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *X;
        X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *y;
        y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *Y;
        Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *h1;
        h1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        fftw_complex *H1;
        H1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

        double *out_1 = new double[blockSize];
        double *in_1  = new double[blockSize];

        if (haux1){    //calcula el h[n] una sola vez a partir de los coeficientes de la ecuacion de diferencias
            sol1[49] = {0};
            haux1=false;
        }

        for(int i=0;i<2048;i++){
            switch (i) {
            case 0:
                h1[i][REAL] =  -0.023650272042231;
                h1[i][IMAG] =   0;
                break;
            case 1:
                h1[i][REAL] =  -0.026618037667980;
                h1[i][IMAG] =   0;
                break;
            case 2:
                h1[i][REAL] = -0.029057374349268;
                h1[i][IMAG] =  0;
                break;
            case 3:
                h1[i][REAL] =  -0.030880288723826;
                h1[i][IMAG] =   0;
                break;
            case 4:
                h1[i][REAL] =  -0.032013396671849;
                h1[i][IMAG] =   0;
                break;
            case 5:
                h1[i][REAL] =  -0.032400569572086;
                h1[i][IMAG] =  0;
                break;
            case 6:
                h1[i][REAL] =  -0.032005107162238;
                h1[i][IMAG] =  0;
                break;
            case 7:
                h1[i][REAL] =  -0.030811357783519;
                h1[i][IMAG] =  0;
                break;
            case 8:
                h1[i][REAL] =  -0.028825723583863;
                h1[i][IMAG] =  0;
                break;
            case 9:
                h1[i][REAL] =  -0.026077007362987;
                h1[i][IMAG] =  0;
                break;
            case 10:
                h1[i][REAL] =  -0.022616078544242;
                h1[i][IMAG] =  0;
                break;
            case 11:
                h1[i][REAL] =  -0.018514857567865;
                h1[i][IMAG] =  0;
                break;
            case 12:
                h1[i][REAL] =  -0.013864640091910;
                h1[i][IMAG] =  0;
                break;
            case 13:
                h1[i][REAL] =  -0.008773804018898;
                h1[i][IMAG] =  0;
                break;
            case 14:
                h1[i][REAL] =  -0.003364962805808;
                h1[i][IMAG] =  0;
                break;
            case 15:
                h1[i][REAL] =  0.002228352934717;
                h1[i][IMAG] =  0;
                break;
            case 16:
                h1[i][REAL] =   0.007865387511180;
                h1[i][IMAG] =   0;
                break;
            case 17:
                h1[i][REAL] =   0.013402114937636;
                h1[i][IMAG] =   0;
                break;
            case 18:
                h1[i][REAL] =   0.01869533605378;
                h1[i][IMAG] =   0;
                break;
            case 19:
                h1[i][REAL] =   0.023606792082434;
                h1[i][IMAG] =   0;
                break;
            case 20:
                h1[i][REAL] =   0.028007167282193;
                h1[i][IMAG] =   0;
                break;
            case 21:
                h1[i][REAL] =   0.031779856061465;
                h1[i][IMAG] =   0;
                break;
            case 22:
                h1[i][REAL] =   0.034824376710089;
                h1[i][IMAG] =   0;
                break;
            case 23:
                h1[i][REAL] =   0.037059324567409;
                h1[i][IMAG] =   0;
                break;
            case 24:
                h1[i][REAL] =   0.038424771639987;
                h1[i][IMAG] =   0;
                break;
            case 25:
                h1[i][REAL] =   0.038884036946932;
                h1[i][IMAG] =   0;
                break;
            case 50:
                h1[i][REAL] =  -0.023650272042231;
                h1[i][IMAG] =   0;
                break;
            case 49:
                h1[i][REAL] =  -0.026618037667980;
                h1[i][IMAG] =   0;
                break;
            case 48:
                h1[i][REAL] = -0.029057374349268;
                h1[i][IMAG] =  0;
                break;
            case 47:
                h1[i][REAL] =  -0.030880288723826;
                h1[i][IMAG] =   0;
                break;
            case 46:
                h1[i][REAL] =  -0.032013396671849;
                h1[i][IMAG] =   0;
                break;
            case 45:
                h1[i][REAL] =  -0.032400569572086;
                h1[i][IMAG] =  0;
                break;
            case 44:
                h1[i][REAL] =  -0.032005107162238;
                h1[i][IMAG] =  0;
                break;
            case 43:
                h1[i][REAL] =  -0.030811357783519;
                h1[i][IMAG] =  0;
                break;
            case 42:
                h1[i][REAL] =  -0.028825723583863;
                h1[i][IMAG] =  0;
                break;
            case 41:
                h1[i][REAL] =  -0.026077007362987;
                h1[i][IMAG] =  0;
                break;
            case 40:
                h1[i][REAL] =  -0.022616078544242;
                h1[i][IMAG] =  0;
                break;
            case 39:
                h1[i][REAL] =  -0.018514857567865;
                h1[i][IMAG] =  0;
                break;
            case 38:
                h1[i][REAL] =  -0.013864640091910;
                h1[i][IMAG] =  0;
                break;
            case 37:
                h1[i][REAL] =  -0.008773804018898;
                h1[i][IMAG] =  0;
                break;
            case 36:
                h1[i][REAL] =  -0.003364962805808;
                h1[i][IMAG] =  0;
                break;
            case 35:
                h1[i][REAL] =  0.002228352934717;
                h1[i][IMAG] =  0;
                break;
            case 34:
                h1[i][REAL] =   0.007865387511180;
                h1[i][IMAG] =   0;
                break;
            case 33:
                h1[i][REAL] =   0.013402114937636;
                h1[i][IMAG] =   0;
                break;
            case 32:
                h1[i][REAL] =   0.01869533605378;
                h1[i][IMAG] =   0;
                break;
            case 31:
                h1[i][REAL] =   0.023606792082434;
                h1[i][IMAG] =   0;
                break;
            case 30:
                h1[i][REAL] =   0.028007167282193;
                h1[i][IMAG] =   0;
                break;
            case 29:
                h1[i][REAL] =   0.031779856061465;
                h1[i][IMAG] =   0;
                break;
            case 28:
                h1[i][REAL] =   0.034824376710089;
                h1[i][IMAG] =   0;
                break;
            case 27:
                h1[i][REAL] =   0.037059324567409;
                h1[i][IMAG] =   0;
                break;
            case 26:
                h1[i][REAL] =   0.038424771639987;
                h1[i][IMAG] =   0;
                break;
            default:
                h1[i][REAL] =  0;
                h1[i][IMAG] =  0;
                break;
            }
        }

        fft(h1,H1); //calcula H[k]

        for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
            if (i<49){
                    x_n[i][REAL]= sol1[i];
                    x_n[i][IMAG]= 0;
                }

            else {
                if(i<(1024 +49)){
                        in_1[i-49]=static_cast<double>(in[i-49]);
                        x_n[i][REAL] = in_1[i-49];
                        x_n[i][IMAG] = 0;
                   }

                else {
                    x_n[i][REAL] = 0;
                    x_n[i][IMAG] = 0;
                }
            }
        }

        for(int i=0;i < 49;i++){  // rellena con la entrada actual
            sol1[i] = in_1[(1024-49) + i];
        }

        fft(x_n,X); //calcula X[k]

        for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
            Y[i][REAL]=H1[i][REAL]*X[i][REAL]-H1[i][IMAG]*X[i][IMAG];
            Y[i][IMAG]=H1[i][REAL]*X[i][IMAG]+H1[i][IMAG]*X[i][REAL];
        }

        idft(Y,y);        //calcula y[n]

        for(int n=0;n<1024;++n){
            out_1[n]= y[n + 49][REAL];
            out_1[n]=(0.02)*(volumeGain)*(out_1[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
            out[n]=static_cast<float>(out_1[n]);// se hace conversion de double a float
        }

        // energia1000=FFT(blockSize,out_1);//se determina la energia de la banda

        delete out_1;
        delete in_1;
        fftw_free(x_n);
        fftw_free(X);
        fftw_free(Y);
        fftw_free(h1);
        fftw_free(H1);

        /*
        for(int i =0; i<blockSize; i++){
            out[i] = 0;
        }
        */

}


void controlVolume::filter_500(int blockSize, int volumeGain, bool inicial, float *in, float *out){

    int N = 2048;

    fftw_complex *x_n;
    x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *h500;
    h500 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *H500;
    H500 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    double *out_500 = new double[blockSize];
    double *in_500  = new double[blockSize];

    if (haux500){
        sol500[100] = {0};
        haux500=false;
    }


    for(int i=0;i<2048;i++){
        switch (i) {
        case 0:
            h500[i][REAL] =   -0.011941756482142;
            h500[i][IMAG] =   0;
            break;
        case 1:
            h500[i][REAL] =  -0.012719031321949;
            h500[i][IMAG] =   0;
            break;
        case 2:
            h500[i][REAL] = -0.013435709228244  ;
            h500[i][IMAG] =  0;
            break;
        case 3:
            h500[i][REAL] =  -0.014085663651104;
            h500[i][IMAG] =   0;
            break;
        case 4:
            h500[i][REAL] =  -0.014663134109047;
            h500[i][IMAG] =   0;
            break;
        case 5:
            h500[i][REAL] =  -0.015162775357471;
            h500[i][IMAG] =  0;
            break;
        case 6:
            h500[i][REAL] =  -0.015579703949401;
            h500[i][IMAG] =  0;
            break;
        case 7:
            h500[i][REAL] =  -0.015909541785888;
            h500[i][IMAG] =  0;
            break;
        case 8:
            h500[i][REAL] = -0.016148456274583;
            h500[i][IMAG] =  0;
            break;
        case 9:
            h500[i][REAL] = -0.016293196739476;
            h500[i][IMAG] =  0;
            break;
        case 10:
            h500[i][REAL] = -0.016341126752486;
            h500[i][IMAG] =  0;
            break;
        case 11:
            h500[i][REAL] = -0.016290252088183;
            h500[i][IMAG] =  0;
            break;
        case 12:
            h500[i][REAL] =  -0.016139244036277;
            h500[i][IMAG] =  0;
            break;
        case 13:
            h500[i][REAL] =  -0.015887457842252;
            h500[i][IMAG] =  0;
            break;
        case 14:
            h500[i][REAL] = -0.015534946084492;
            h500[i][IMAG] =  0;
            break;
        case 15:
            h500[i][REAL] = -0.015082466835963;
            h500[i][IMAG] =  0;
            break;
        case 16:
            h500[i][REAL] = -0.014531486499872 ;
            h500[i][IMAG] = 0;
            break;
        case 17:
            h500[i][REAL] = -0.013884177251126;
            h500[i][IMAG] = 0;
            break;
        case 18:
            h500[i][REAL] = -0.013143409058780;
            h500[i][IMAG] = 0;
            break;
        case 19:
            h500[i][REAL] = -0.012312736308353;
            h500[i][IMAG] = 0;
            break;
        case 20:
            h500[i][REAL] = -0.011396379086820;
            h500[i][IMAG] = 0;
            break;
        case 21:
            h500[i][REAL] = -0.010399199236621;
            h500[i][IMAG] = 0;
            break;
        case 22:
            h500[i][REAL] = -0.009326671328000;
            h500[i][IMAG] = 0;
            break;
        case 23:
            h500[i][REAL] = -0.008184848740949;
            h500[i][IMAG] = 0;
            break;
        case 24:
            h500[i][REAL] = -0.006980325088576;
            h500[i][IMAG] = 0;
            break;
        case 25:
            h500[i][REAL] = -0.005720191252632;
            h500[i][IMAG] = 0;
            break;
        case 26:
            h500[i][REAL] = -0.004411988338748;
            h500[i][IMAG] = 0;
            break;
        case 27:
            h500[i][REAL] = -0.003063656893404;
            h500[i][IMAG] = 0;
            break;
        case 28:
            h500[i][REAL] = -0.001683482756497;
            h500[i][IMAG] = 0;
            break;
        case 29:
            h500[i][REAL] = -0.000280039952229;
            h500[i][IMAG] = 0;
            break;
        case 30:
            h500[i][REAL] = 0.001137868953242;
            h500[i][IMAG] = 0;
            break;
        case 31:
            h500[i][REAL] = 0.002561274576663;
            h500[i][IMAG] = 0;
            break;
        case 32:
            h500[i][REAL] = 0.003981104055680;
            h500[i][IMAG] = 0;
            break;
        case 33:
            h500[i][REAL] = 0.005388245515643;
            h500[i][IMAG] = 0;
            break;
        case 34:
            h500[i][REAL] = 0.006773613169977;
            h500[i][IMAG] = 0;
            break;
        case 35:
            h500[i][REAL] = 0.008128212551858;
            h500[i][IMAG] = 0;
            break;
        case 36:
            h500[i][REAL] = 0.009443205372090;
            h500[i][IMAG] = 0;
            break;
        case 37:
            h500[i][REAL] = 0.010709973499306;
            h500[i][IMAG] = 0;
            break;
        case 38:
            h500[i][REAL] = 0.011920181564002;
            h500[i][IMAG] = 0;
            break;
        case 39:
            h500[i][REAL] = 0.013065837697338;
            h500[i][IMAG] = 0;
            break;
        case 40:
            h500[i][REAL] = 0.014139351929149;
            h500[i][IMAG] = 0;
            break;
        case 41:
            h500[i][REAL] = 0.015133591786964;
            h500[i][IMAG] = 0;
            break;
        case 42:
            h500[i][REAL] = 0.016041934659007;
            h500[i][IMAG] = 0;
            break;
        case 43:
            h500[i][REAL] = 0.016858316508920;
            h500[i][IMAG] = 0;
            break;
        case 44:
            h500[i][REAL] = 0.017577276558116;
            h500[i][IMAG] = 0;
            break;
        case 45:
            h500[i][REAL] = 0.018193997583028;
            h500[i][IMAG] = 0;
            break;
        case 46:
            h500[i][REAL] = 0.018704341508784;
            h500[i][IMAG] = 0;
            break;
        case 47:
            h500[i][REAL] = 0.019104880017753;
            h500[i][IMAG] = 0;
            break;
        case 48:
            h500[i][REAL] = 0.019392919930647;
            h500[i][IMAG] = 0;
            break;
        case 49:
            h500[i][REAL] = 0.019566523159108;
            h500[i][IMAG] = 0;
            break;
        case 50:
            h500[i][REAL] = 0.019624521071662;
            h500[i][IMAG] = 0;
            break;

        case 100:
            h500[i][REAL] =   -0.011941756482142;
            h500[i][IMAG] =   0;
            break;
        case 99:
            h500[i][REAL] =  -0.012719031321949;
            h500[i][IMAG] =   0;
            break;
        case 98:
            h500[i][REAL] = -0.013435709228244  ;
            h500[i][IMAG] =  0;
            break;
        case 97:
            h500[i][REAL] =  -0.014085663651104;
            h500[i][IMAG] =   0;
            break;
        case 96:
            h500[i][REAL] =  -0.014663134109047;
            h500[i][IMAG] =   0;
            break;
        case 95:
            h500[i][REAL] =  -0.015162775357471;
            h500[i][IMAG] =  0;
            break;
        case 94:
            h500[i][REAL] =  -0.015579703949401;
            h500[i][IMAG] =  0;
            break;
        case 93:
            h500[i][REAL] =  -0.015909541785888;
            h500[i][IMAG] =  0;
            break;
        case 92:
            h500[i][REAL] = -0.016148456274583;
            h500[i][IMAG] =  0;
            break;
        case 91:
            h500[i][REAL] = -0.016293196739476;
            h500[i][IMAG] =  0;
            break;
        case 90:
            h500[i][REAL] = -0.016341126752486;
            h500[i][IMAG] =  0;
            break;
        case 89:
            h500[i][REAL] = -0.016290252088183;
            h500[i][IMAG] =  0;
            break;
        case 88:
            h500[i][REAL] =  -0.016139244036277;
            h500[i][IMAG] =  0;
            break;
        case 87:
            h500[i][REAL] =  -0.015887457842252;
            h500[i][IMAG] =  0;
            break;
        case 86:
            h500[i][REAL] = -0.015534946084492;
            h500[i][IMAG] =  0;
            break;
        case 85:
            h500[i][REAL] = -0.015082466835963;
            h500[i][IMAG] =  0;
            break;
        case 84:
            h500[i][REAL] = -0.014531486499872 ;
            h500[i][IMAG] = 0;
            break;
        case 83:
            h500[i][REAL] = -0.013884177251126;
            h500[i][IMAG] = 0;
            break;
        case 82:
            h500[i][REAL] = -0.013143409058780;
            h500[i][IMAG] = 0;
            break;
        case 81:
            h500[i][REAL] = -0.012312736308353;
            h500[i][IMAG] = 0;
            break;
        case 80:
            h500[i][REAL] = -0.011396379086820;
            h500[i][IMAG] = 0;
            break;
        case 79:
            h500[i][REAL] = -0.010399199236621;
            h500[i][IMAG] = 0;
            break;
        case 78:
            h500[i][REAL] = -0.009326671328000;
            h500[i][IMAG] = 0;
            break;
        case 77:
            h500[i][REAL] = -0.008184848740949;
            h500[i][IMAG] = 0;
            break;
        case 76:
            h500[i][REAL] = -0.006980325088576;
            h500[i][IMAG] = 0;
            break;
        case 75:
            h500[i][REAL] = -0.005720191252632;
            h500[i][IMAG] = 0;
            break;
        case 74:
            h500[i][REAL] = -0.004411988338748;
            h500[i][IMAG] = 0;
            break;
        case 73:
            h500[i][REAL] = -0.003063656893404;
            h500[i][IMAG] = 0;
            break;
        case 72:
            h500[i][REAL] = -0.001683482756497;
            h500[i][IMAG] = 0;
            break;
        case 71:
            h500[i][REAL] = -0.000280039952229;
            h500[i][IMAG] = 0;
            break;
        case 70:
            h500[i][REAL] = 0.001137868953242;
            h500[i][IMAG] = 0;
            break;
        case 69:
            h500[i][REAL] = 0.002561274576663;
            h500[i][IMAG] = 0;
            break;
        case 68:
            h500[i][REAL] = 0.003981104055680;
            h500[i][IMAG] = 0;
            break;
        case 67:
            h500[i][REAL] = 0.005388245515643;
            h500[i][IMAG] = 0;
            break;
        case 66:
            h500[i][REAL] = 0.006773613169977;
            h500[i][IMAG] = 0;
            break;
        case 65:
            h500[i][REAL] = 0.008128212551858;
            h500[i][IMAG] = 0;
            break;
        case 64:
            h500[i][REAL] = 0.009443205372090;
            h500[i][IMAG] = 0;
            break;
        case 63:
            h500[i][REAL] = 0.010709973499306;
            h500[i][IMAG] = 0;
            break;
        case 62:
            h500[i][REAL] = 0.011920181564002;
            h500[i][IMAG] = 0;
            break;
        case 61:
            h500[i][REAL] = 0.013065837697338;
            h500[i][IMAG] = 0;
            break;
        case 60:
            h500[i][REAL] = 0.014139351929149;
            h500[i][IMAG] = 0;
            break;
        case 59:
            h500[i][REAL] = 0.015133591786964;
            h500[i][IMAG] = 0;
            break;
        case 58:
            h500[i][REAL] = 0.016041934659007;
            h500[i][IMAG] = 0;
            break;
        case 57:
            h500[i][REAL] = 0.016858316508920;
            h500[i][IMAG] = 0;
            break;
        case 56:
            h500[i][REAL] = 0.017577276558116;
            h500[i][IMAG] = 0;
            break;
        case 55:
            h500[i][REAL] = 0.018193997583028;
            h500[i][IMAG] = 0;
            break;
        case 54:
            h500[i][REAL] = 0.018704341508784;
            h500[i][IMAG] = 0;
            break;
        case 53:
            h500[i][REAL] = 0.019104880017753;
            h500[i][IMAG] = 0;
            break;
        case 52:
            h500[i][REAL] = 0.019392919930647;
            h500[i][IMAG] = 0;
            break;
        case 51:
            h500[i][REAL] = 0.019566523159108;
            h500[i][IMAG] = 0;
            break;
        default:
            h500[i][REAL] =  0;
            h500[i][IMAG] =  0;
            break;
        }
    }


    fft(h500,H500); //calcula H[k]

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<99){
                x_n[i][REAL]= sol500[i];
                x_n[i][IMAG]= 0;
            }

        else {
            if(i<(1024 + 99)){
                    in_500[i-99]=static_cast<double>(in[i-99]);
                    x_n[i][REAL] = in_500[i-99];
                    x_n[i][IMAG] = 0;
               }
            else {
                x_n[i][REAL] = 0;
                x_n[i][IMAG] = 0;
            }
        }
    }

    for(int i=0;i < 99;i++){  // rellena con la entrada actual
        sol500[i] = in_500[(1024-99) + i];
    }

    fft(x_n,X); //calcula X[k]

    for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
        Y[i][REAL]=H500[i][REAL]*X[i][REAL]-H500[i][IMAG]*X[i][IMAG];
        Y[i][IMAG]=H500[i][REAL]*X[i][IMAG]+H500[i][IMAG]*X[i][REAL];
    }

    idft(Y,y);        //calcula y[n]

    for(int n=0;n<1024;++n){
        out_500[n]= y[n + 99][REAL];
        out_500[n]=(0.02)*(volumeGain)*(out_500[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
        out[n]=static_cast<float>(out_500[n]);// se hace conversion de double a float
    }

    //energia500=FFT(blockSize,out_6);

    delete out_500;
    delete in_500;
    fftw_free(x_n);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h500);
    fftw_free(H500);

    /*
    for(int i =0; i<blockSize; i++){
        out[i] = 0;
    }
    */
}


/*
//-------------------------------------------------FILTRO DE 250Hz------------------------------------------------------------------------------//
void controlVolume::filter_250(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 250Hz
    double s_250=0.015896786349464623;
    double a_0_250=-1.9797045895876129;
    double a_1_250=0.98206529318049718;
    double b_0_250=0.44000811887250152;
    double a_2_250=-1.9897991375567881;
    double a_3_250=0.99046350451803944;
    double b_1_250=-1.9999997456454339;

    double *tmpout_9=new double[blockSize];
    double *tmpout_0=new double[blockSize];
    double *z=new double[blockSize];
    double *out_4=new double[blockSize];



    if(inicial){
        _debug("filtro de 250Hz---------------------------------------------------------------------------------------------------" << std::endl);
        for(int n=0;n<blockSize;++n){

            if(n==0){
                tmpout_9[n]=(s_250)*in[n];
                z[n]=tmpout_9[n];
                tmpout_0[n]=(s_250)*z[n];
                out_4[n]=tmpout_0[n];
            }
            else if(n==1){
                tmpout_9[n]=(s_250)*in[n]-(a_0_250)*tmpout_9[n-1];
                z[n]=tmpout_9[n]+(b_0_250)*tmpout_9[n-1];
                tmpout_0[n]=(s_250)*z[n]-(a_2_250)*tmpout_0[n-1];
                out_4[n]=tmpout_0[n]+(b_1_250)*tmpout_0[n-1];
            }
            else{
                tmpout_9[n]=(s_250)*in[n]-(a_0_250)*tmpout_9[n-1]-(a_1_250)*tmpout_9[n-2];
                z[n]=tmpout_9[n]+(b_0_250)*tmpout_9[n-1]+tmpout_9[n-2];
                tmpout_0[n]=(s_250)*z[n]-(a_2_250)*tmpout_0[n-1]-(a_3_250)*tmpout_0[n-2];
                out_4[n]=tmpout_0[n]+(b_1_250)*tmpout_0[n-1]+tmpout_0[n-2];
            }
            out_4[n]=(0.02)*(volumeGain)*out_4[n];
            out[n]=static_cast<float>(out_4[n]);
        }

     }
    else{

            for(int n=0;n<blockSize;++n){

                if(n==0){
                    tmpout_9[n]=(s_250)*in[n]-(a_0_250)*j-(a_1_250)*i;
                    z[n]=tmpout_9[n]+(b_0_250)*j+i;
                    tmpout_0[n]=(s_250)*z[n]-(a_2_250)*l-(a_3_250)*k;
                    out_4[n]=tmpout_0[n]+(b_1_250)*l+k;
                }
                else if(n==1){
                    tmpout_9[n]=(s_250)*in[n]-(a_0_250)*tmpout_9[n-1]-(a_1_250)*j;
                    z[n]=tmpout_9[n]+(b_0_250)*tmpout_9[n-1]+j;
                    tmpout_0[n]=(s_250)*z[n]-(a_2_250)*tmpout_0[n-1]-(a_3_250)*l;
                    out_4[n]=tmpout_0[n]+(b_1_250)*tmpout_0[n-1]+l;
                }
                else{
                    tmpout_9[n]=(s_250)*in[n]-(a_0_250)*tmpout_9[n-1]-(a_1_250)*tmpout_9[n-2];
                    z[n]=tmpout_9[n]+(b_0_250)*tmpout_9[n-1]+tmpout_9[n-2];
                    tmpout_0[n]=(s_250)*z[n]-(a_2_250)*tmpout_0[n-1]-(a_3_250)*tmpout_0[n-2];
                    out_4[n]=tmpout_0[n]+(b_1_250)*tmpout_0[n-1]+tmpout_0[n-2];
                }
                out_4[n]=(0.02)*(volumeGain)*out_4[n];
                out[n]=static_cast<float>(out_4[n]);
            }
       }

      energia250=FFT(blockSize,out_4);
      i=tmpout_9[1022];//w1(-2)
      j=tmpout_9[1023];//w1(-1)
      k=tmpout_0[1022];//w2(-2)
      l=tmpout_0[1023];//w2(-1)

      delete tmpout_9;
      delete tmpout_0;
      delete z;
      delete out_4;

} */

//-------------------------------------------------FILTRO DE 250Hz------------------------------------------------------------------------------//
void controlVolume::filter_250(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 250Hz
    /*
    int N = 2048;

    fftw_complex *x_n;
    x_n = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *h250;
    h250 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex *H250;
    H250 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    double *out_250 = new double[blockSize];
    double *in_250  = new double[blockSize];

    if (haux250){    //calcula el h[n] una sola vez a partir de los coeficientes de la ecuacion de diferencias
        sol250[199] = {0};
        haux250=false;
    }

    for(int i=0;i<2048;i++){
        switch (i) {
        case 0:
            h250[i][REAL] = -0.003492796533952;
            h250[i][IMAG] =   0;
            break;
        case 1:
            h250[i][REAL] = -0.003725541579766;
            h250[i][IMAG] =   0;
            break;
        case 2:
            h250[i][REAL] = -0.003955562564568;
            h250[i][IMAG] =  0;
            break;
        case 3:
            h250[i][REAL] = -0.004182435240641;
            h250[i][IMAG] =   0;
            break;
        case 4:
            h250[i][REAL] = -0.004405736663093;
            h250[i][IMAG] =   0;
            break;
        case 5:
            h250[i][REAL] = -0.004625046007301;
            h250[i][IMAG] =  0;
            break;
        case 6:
            h250[i][REAL] = -0.004839945388123;
            h250[i][IMAG] =  0;
            break;
        case 7:
            h250[i][REAL] = -0.005050020679241;
            h250[i][IMAG] =  0;
            break;
        case 8:
            h250[i][REAL] = -0.005254862330963;
            h250[i][IMAG] =  0;
            break;
        case 9:
            h250[i][REAL] = -0.005454066184830;
            h250[i][IMAG] =  0;
            break;
        case 10:
            h250[i][REAL] = -0.005647234283361;
            h250[i][IMAG] =  0;
            break;
        case 11:
            h250[i][REAL] = -0.005833975673277;
            h250[i][IMAG] =  0;
            break;
        case 12:
            h250[i][REAL] = -0.006013907200536;
            h250[i][IMAG] =  0;
            break;
        case 13:
            h250[i][REAL] = -0.006186654295549;
            h250[i][IMAG] =  0;
            break;
        case 14:
            h250[i][REAL] = -0.006351851746919;
            h250[i][IMAG] =  0;
            break;
        case 15:
            h250[i][REAL] = -0.006509144462102;
            h250[i][IMAG] =  0;
            break;
        case 16:
            h250[i][REAL] = -0.006658188213376;
            h250[i][IMAG] =   0;
            break;
        case 17:
            h250[i][REAL] = -0.006798650367548;
            h250[i][IMAG] =   0;
            break;
        case 18:
            h250[i][REAL] = -0.006930210597838;
            h250[i][IMAG] =   0;
            break;
        case 19:
            h250[i][REAL] = -0.007052561576431;
            h250[i][IMAG] =   0;
            break;
        case 20:
            h250[i][REAL] = -0.007165409646189;
            h250[i][IMAG] =   0;
            break;
        case 21:
            h250[i][REAL] = -0.007268475470083;
            h250[i][IMAG] =   0;
            break;
        case 22:
            h250[i][REAL] = -0.007361494656927;
            h250[i][IMAG] =   0;
            break;
        case 23:
            h250[i][REAL] = -0.007444218362035;
            h250[i][IMAG] =   0;
            break;
        case 24:
            h250[i][REAL] = -0.007516413861479;
            h250[i][IMAG] =   0;
            break;
        case 25:
            h250[i][REAL] = -0.007577865098673;
            h250[i][IMAG] =   0;
            break;
        case 26:
            h250[i][REAL] = -0.007628373202039;
            h250[i][IMAG] =   0;
            break;
        case 27:
            h250[i][REAL] = -0.007667756972603;
            h250[i][IMAG] =   0;
            break;
        case 28:
            h250[i][REAL] = -0.007695853340380;
            h250[i][IMAG] =   0;
            break;
        case 29:
            h250[i][REAL] = -0.007712517788506;
            h250[i][IMAG] =   0;
            break;
        case 30:
            h250[i][REAL] = -0.007717624744110;
            h250[i][IMAG] =   0;
            break;
        case 31:
            h250[i][REAL] = -0.007711067934994;
            h250[i][IMAG] =   0;
            break;
        case 32:
            h250[i][REAL] = -0.007692760711255;
            h250[i][IMAG] =   0;
            break;
        case 33:
            h250[i][REAL] = -0.007662636331045;
            h250[i][IMAG] =   0;
            break;
        case 34:
            h250[i][REAL] = -0.007620648209735;
            h250[i][IMAG] =   0;
            break;
        case 35:
            h250[i][REAL] = -0.007566770131835;
            h250[i][IMAG] =   0;
            break;
        case 36:
            h250[i][REAL] = -0.007500996425070;
            h250[i][IMAG] =   0;
            break;
        case 37:
            h250[i][REAL] = -0.007423342096109;
            h250[i][IMAG] =   0;
            break;
        case 38:
            h250[i][REAL] = -0.007333842927503;
            h250[i][IMAG] =   0;
            break;
        case 39:
            h250[i][REAL] = -0.007232555535483;
            h250[i][IMAG] =   0;
            break;
        case 40:
            h250[i][REAL] = -0.007119557388325;
            h250[i][IMAG] =   0;
            break;
        case 41:
            h250[i][REAL] = -0.006994946785089;
            h250[i][IMAG] =   0;
            break;
        case 42:
            h250[i][REAL] = -0.006858842794601;
            h250[i][IMAG] =   0;
            break;
        case 43:
            h250[i][REAL] = -0.006711385154638;
            h250[i][IMAG] =   0;
            break;
        case 44:
            h250[i][REAL] = -0.006552734131354;
            h250[i][IMAG] =   0;
            break;
        case 45:
            h250[i][REAL] = -0.006383070339060;
            h250[i][IMAG] =   0;
            break;
        case 46:
            h250[i][REAL] = -0.006202594520561;
            h250[i][IMAG] =   0;
            break;
        case 47:
            h250[i][REAL] = -0.006011527288317;
            h250[i][IMAG] =   0;
            break;
        case 48:
            h250[i][REAL] = -0.005810108826802;
            h250[i][IMAG] =   0;
            break;
        case 49:
            h250[i][REAL] = -0.005598598556479;
            h250[i][IMAG] =   0;
            break;
        case 50:
            h250[i][REAL] = -0.005377274759917;
            h250[i][IMAG] =   0;
            break;
        case 51:
            h250[i][REAL] = -0.005146434170647;
            h250[i][IMAG] =   0;
            break;
        case 52:
            h250[i][REAL] = -0.004906391525414;
            h250[i][IMAG] =   0;
            break;
        case 53:
            h250[i][REAL] = -0.004657479080577;
            h250[i][IMAG] =   0;
            break;
        case 54:
            h250[i][REAL] = -0.004400046093480;
            h250[i][IMAG] =   0;
            break;
        case 55:
            h250[i][REAL] = -0.004134458269679;
            h250[i][IMAG] =   0;
            break;
        case 56:
            h250[i][REAL] = -0.003861097176999;
            h250[i][IMAG] =   0;
            break;
        case 57:
            h250[i][REAL] = -0.003580359627434;
            h250[i][IMAG] =   0;
            break;
        case 58:
            h250[i][REAL] = -0.003292657028015;
            h250[i][IMAG] =   0;
            break;
        case 59:
            h250[i][REAL] = -0.002998414701791;
            h250[i][IMAG] =   0;
            break;
        case 60:
            h250[i][REAL] = -0.002698071180161;
            h250[i][IMAG] =   0;
            break;
        case 61:
            h250[i][REAL] = -0.002392077467845;
            h250[i][IMAG] =   0;
            break;
        case 62:
            h250[i][REAL] = -0.002080896281838;
            h250[i][IMAG] =   0;
            break;
        case 63:
            h250[i][REAL] = -0.001765001265747;
            h250[i][IMAG] =   0;
            break;
        case 64:
            h250[i][REAL] = -0.001444876180976;
            h250[i][IMAG] =   0;
            break;
        case 65:
            h250[i][REAL] = -0.001121014076257;
            h250[i][IMAG] =   0;
            break;
        case 66:
            h250[i][REAL] = -0.000793916437071;
            h250[i][IMAG] =   0;
            break;
        case 67:
            h250[i][REAL] = -0.000464092316577;
            h250[i][IMAG] =   0;
            break;
        case 68:
            h250[i][REAL] = -0.000132057449656;
            h250[i][IMAG] =   0;
            break;
        case 69:
            h250[i][REAL] = 0.000201666648229;
            h250[i][IMAG] =   0;
            break;
        case 70:
            h250[i][REAL] = 0.000536553595665;
            h250[i][IMAG] =   0;
            break;
        case 71:
            h250[i][REAL] = 0.000872073071679;
            h250[i][IMAG] =   0;
            break;
        case 72:
            h250[i][REAL] = 0.001207691752098;
            h250[i][IMAG] =   0;
            break;
        case 73:
            h250[i][REAL] = 0.001542874254124;
            h250[i][IMAG] =   0;
            break;
        case 74:
            h250[i][REAL] = 0.001877084087277;
            h250[i][IMAG] =   0;
            break;
        case 75:
            h250[i][REAL] = 0.002209784608918;
            h250[i][IMAG] =   0;
            break;
        case 76:
            h250[i][REAL] = 0.002540439982485;
            h250[i][IMAG] =   0;
            break;
        case 77:
            h250[i][REAL] = 0.002868516136617;
            h250[i][IMAG] =   0;
            break;
        case 78:
            h250[i][REAL] = 0.003193481723293;
            h250[i][IMAG] =   0;
            break;
        case 79:
            h250[i][REAL] = 0.003514809073163;
            h250[i][IMAG] =   0;
            break;
        case 80:
            h250[i][REAL] = 0.003831975146184;
            h250[i][IMAG] =   0;
            break;
        case 81:
            h250[i][REAL] = 0.004144462475753;
            h250[i][IMAG] =   0;
            break;
        case 82:
            h250[i][REAL] = 0.004451760104474;
            h250[i][IMAG] =   0;
            break;
        case 83:
            h250[i][REAL] = 0.004753364509754;
            h250[i][IMAG] =   0;
            break;
        case 84:
            h250[i][REAL] = 0.005048780517428;
            h250[i][IMAG] =   0;
            break;
       case 85:
            h250[i][REAL] = 0.005337522201614;
            h250[i][IMAG] =   0;
            break;
        case 86:
            h250[i][REAL] = 0.005619113769061;
            h250[i][IMAG] =  0;
            break;
        case 87:
            h250[i][REAL] = 0.005893090426249;
            h250[i][IMAG] =  0;
            break;
        case 88:
            h250[i][REAL] = 0.006158999227546;
            h250[i][IMAG] =  0;
            break;
        case 89:
            h250[i][REAL] = 0.006416399902769;
            h250[i][IMAG] =  0;
            break;
        case 90:
            h250[i][REAL] = 0.006664865662513;
            h250[i][IMAG] =  0;
            break;
        case 91:
            h250[i][REAL] = 0.006903983979682;
            h250[i][IMAG] =  0;
            break;
        case 92:
            h250[i][REAL] = 0.007133357345681;
            h250[i][IMAG] =  0;
            break;
        case 93:
            h250[i][REAL] = 0.007352603999776;
            h250[i][IMAG] =  0;
            break;
        case 94:
            h250[i][REAL] = 0.007561358630206;
            h250[i][IMAG] =  0;
            break;
        case 95:
            h250[i][REAL] = 0.007759273045649;
            h250[i][IMAG] =  0;
            break;
        case 96:
            h250[i][REAL] = 0.007946016815726;
            h250[i][IMAG] =  0;
            break;
        case 97:
            h250[i][REAL] = 0.008121277879292;
            h250[i][IMAG] =   0;
            break;
        case 98:
            h250[i][REAL] = 0.008284763119299;
            h250[i][IMAG] =   0;
            break;
        case 99:
            h250[i][REAL] = 0.008436198903098;
            h250[i][IMAG] =  0;
            break;
        case 100:
            h250[i][REAL] = 0.008575331587137;
            h250[i][IMAG] =   0;
            break;
        case 101:
            h250[i][REAL] = 0.008701927985023;
            h250[i][IMAG] =   0;
            break;
        case 102:
            h250[i][REAL] = 0.008815775798054;
            h250[i][IMAG] =  0;
            break;
        case 103:
            h250[i][REAL] = 0.008916684007356;
            h250[i][IMAG] =  0;
            break;
        case 104:
            h250[i][REAL] = 0.009004483226852;
            h250[i][IMAG] =  0;
            break;
        case 105:
            h250[i][REAL] = 0.009079026016366;
            h250[i][IMAG] =  0;
            break;
        case 106:
            h250[i][REAL] = 0.009140187154248;
            h250[i][IMAG] =   0;
            break;
        case 107:
            h250[i][REAL] = 0.009187863868967;
            h250[i][IMAG] =   0;
            break;
        case 108:
            h250[i][REAL] = 0.009221976029235;
            h250[i][IMAG] =  0;
            break;
        case 109:
            h250[i][REAL] = 0.009242466292264;
            h250[i][IMAG] =   0;
            break;
        case 110:
            h250[i][REAL] = 0.009249300209878;
            h250[i][IMAG] =   0;
            break;



        case 220:
            h250[i][REAL] = -0.003492796533952;
            h250[i][IMAG] =   0;
            break;
        case 219:
            h250[i][REAL] = -0.003725541579766;
            h250[i][IMAG] =   0;
            break;
        case 218:
            h250[i][REAL] = -0.003955562564568;
            h250[i][IMAG] =  0;
            break;
        case 217:
            h250[i][REAL] = -0.004182435240641;
            h250[i][IMAG] =   0;
            break;
        case 216:
            h250[i][REAL] = -0.004405736663093;
            h250[i][IMAG] =   0;
            break;
        case 215:
            h250[i][REAL] = -0.004625046007301;
            h250[i][IMAG] =  0;
            break;
        case 214:
            h250[i][REAL] = -0.004839945388123;
            h250[i][IMAG] =  0;
            break;
        case 213:
            h250[i][REAL] = -0.005050020679241;
            h250[i][IMAG] =  0;
            break;
        case 212:
            h250[i][REAL] = -0.005254862330963;
            h250[i][IMAG] =  0;
            break;
        case 211:
            h250[i][REAL] = -0.005454066184830;
            h250[i][IMAG] =  0;
            break;
        case 210:
            h250[i][REAL] = -0.005647234283361;
            h250[i][IMAG] =  0;
            break;
        case 209:
            h250[i][REAL] = -0.005833975673277;
            h250[i][IMAG] =  0;
            break;
        case 208:
            h250[i][REAL] = -0.006013907200536;
            h250[i][IMAG] =  0;
            break;
        case 207:
            h250[i][REAL] = -0.006186654295549;
            h250[i][IMAG] =  0;
            break;
        case 206:
            h250[i][REAL] = -0.006351851746919;
            h250[i][IMAG] =  0;
            break;
        case 205:
            h250[i][REAL] = -0.006509144462102;
            h250[i][IMAG] =  0;
            break;
        case 204:
            h250[i][REAL] = -0.006658188213376;
            h250[i][IMAG] =   0;
            break;
        case 203:
            h250[i][REAL] = -0.006798650367548;
            h250[i][IMAG] =   0;
            break;
        case 202:
            h250[i][REAL] = -0.006930210597838;
            h250[i][IMAG] =   0;
            break;
        case 201:
            h250[i][REAL] = -0.007052561576431;
            h250[i][IMAG] =   0;
            break;
        case 200:
            h250[i][REAL] = -0.007165409646189;
            h250[i][IMAG] =   0;
            break;
        case 199:
            h250[i][REAL] = -0.007268475470083;
            h250[i][IMAG] =   0;
            break;
        case 198:
            h250[i][REAL] = -0.007361494656927;
            h250[i][IMAG] =   0;
            break;
        case 197:
            h250[i][REAL] = -0.007444218362035;
            h250[i][IMAG] =   0;
            break;
        case 196:
            h250[i][REAL] = -0.007516413861479;
            h250[i][IMAG] =   0;
            break;
        case 195:
            h250[i][REAL] = -0.007577865098673;
            h250[i][IMAG] =   0;
            break;
        case 194:
            h250[i][REAL] = -0.007628373202039;
            h250[i][IMAG] =   0;
            break;
        case 193:
            h250[i][REAL] = -0.007667756972603;
            h250[i][IMAG] =   0;
            break;
        case 192:
            h250[i][REAL] = -0.007695853340380;
            h250[i][IMAG] =   0;
            break;
        case 191:
            h250[i][REAL] = -0.007712517788506;
            h250[i][IMAG] =   0;
            break;
        case 190:
            h250[i][REAL] = -0.007717624744110;
            h250[i][IMAG] =   0;
            break;
        case 189:
            h250[i][REAL] = -0.007711067934994;
            h250[i][IMAG] =   0;
            break;
        case 188:
            h250[i][REAL] = -0.007692760711255;
            h250[i][IMAG] =   0;
            break;
        case 187:
            h250[i][REAL] = -0.007662636331045;
            h250[i][IMAG] =   0;
            break;
        case 186:
            h250[i][REAL] = -0.007620648209735;
            h250[i][IMAG] =   0;
            break;
        case 185:
            h250[i][REAL] = -0.007566770131835;
            h250[i][IMAG] =   0;
            break;
        case 184:
            h250[i][REAL] = -0.007500996425070;
            h250[i][IMAG] =   0;
            break;
        case 183:
            h250[i][REAL] = -0.007423342096109;
            h250[i][IMAG] =   0;
            break;
        case 182:
            h250[i][REAL] = -0.007333842927503;
            h250[i][IMAG] =   0;
            break;
        case 181:
            h250[i][REAL] = -0.007232555535483;
            h250[i][IMAG] =   0;
            break;
        case 180:
            h250[i][REAL] = -0.007119557388325;
            h250[i][IMAG] =   0;
            break;
        case 179:
            h250[i][REAL] = -0.006994946785089;
            h250[i][IMAG] =   0;
            break;
        case 178:
            h250[i][REAL] = -0.006858842794601;
            h250[i][IMAG] =   0;
            break;
        case 177:
            h250[i][REAL] = -0.006711385154638;
            h250[i][IMAG] =   0;
            break;
        case 176:
            h250[i][REAL] = -0.006552734131354;
            h250[i][IMAG] =   0;
            break;
        case 175:
            h250[i][REAL] = -0.006383070339060;
            h250[i][IMAG] =   0;
            break;
        case 174:
            h250[i][REAL] = -0.006202594520561;
            h250[i][IMAG] =   0;
            break;
        case 173:
            h250[i][REAL] = -0.006011527288317;
            h250[i][IMAG] =   0;
            break;
        case 172:
            h250[i][REAL] = -0.005810108826802;
            h250[i][IMAG] =   0;
            break;
        case 171:
            h250[i][REAL] = -0.005598598556479;
            h250[i][IMAG] =   0;
            break;
        case 170:
            h250[i][REAL] = -0.005377274759917;
            h250[i][IMAG] =   0;
            break;
        case 169:
            h250[i][REAL] = -0.005146434170647;
            h250[i][IMAG] =   0;
            break;
        case 168:
            h250[i][REAL] = -0.004906391525414;
            h250[i][IMAG] =   0;
            break;
        case 167:
            h250[i][REAL] = -0.004657479080577;
            h250[i][IMAG] =   0;
            break;
        case 166:
            h250[i][REAL] = -0.004400046093480;
            h250[i][IMAG] =   0;
            break;
        case 165:
            h250[i][REAL] = -0.004134458269679;
            h250[i][IMAG] =   0;
            break;
        case 164:
            h250[i][REAL] = -0.003861097176999;
            h250[i][IMAG] =   0;
            break;
        case 163:
            h250[i][REAL] = -0.003580359627434;
            h250[i][IMAG] =   0;
            break;
        case 162:
            h250[i][REAL] = -0.003292657028015;
            h250[i][IMAG] =   0;
            break;
        case 161:
            h250[i][REAL] = -0.002998414701791;
            h250[i][IMAG] =   0;
            break;
        case 160:
            h250[i][REAL] = -0.002698071180161;
            h250[i][IMAG] =   0;
            break;
        case 159:
            h250[i][REAL] = -0.002392077467845;
            h250[i][IMAG] =   0;
            break;
        case 158:
            h250[i][REAL] = -0.002080896281838;
            h250[i][IMAG] =   0;
            break;
        case 157:
            h250[i][REAL] = -0.001765001265747;
            h250[i][IMAG] =   0;
            break;
        case 156:
            h250[i][REAL] = -0.001444876180976;
            h250[i][IMAG] =   0;
            break;
        case 155:
            h250[i][REAL] = -0.001121014076257;
            h250[i][IMAG] =   0;
            break;
        case 154:
            h250[i][REAL] = -0.000793916437071;
            h250[i][IMAG] =   0;
            break;
        case 153:
            h250[i][REAL] = -0.000464092316577;
            h250[i][IMAG] =   0;
            break;
        case 152:
            h250[i][REAL] = -0.000132057449656;
            h250[i][IMAG] =   0;
            break;
        case 151:
            h250[i][REAL] = 0.000201666648229;
            h250[i][IMAG] =   0;
            break;
        case 150:
            h250[i][REAL] = 0.000536553595665;
            h250[i][IMAG] =   0;
            break;
        case 149:
            h250[i][REAL] = 0.000872073071679;
            h250[i][IMAG] =   0;
            break;
        case 148:
            h250[i][REAL] = 0.001207691752098;
            h250[i][IMAG] =   0;
            break;
        case 147:
            h250[i][REAL] = 0.001542874254124;
            h250[i][IMAG] =   0;
            break;
        case 146:
            h250[i][REAL] = 0.001877084087277;
            h250[i][IMAG] =   0;
            break;
        case 145:
            h250[i][REAL] = 0.002209784608918;
            h250[i][IMAG] =   0;
            break;
        case 144:
            h250[i][REAL] = 0.002540439982485;
            h250[i][IMAG] =   0;
            break;
        case 143:
            h250[i][REAL] = 0.002868516136617;
            h250[i][IMAG] =   0;
            break;
        case 142:
            h250[i][REAL] = 0.003193481723293;
            h250[i][IMAG] =   0;
            break;
        case 141:
            h250[i][REAL] = 0.003514809073163;
            h250[i][IMAG] =   0;
            break;
        case 140:
            h250[i][REAL] = 0.003831975146184;
            h250[i][IMAG] =   0;
            break;
        case 139:
            h250[i][REAL] = 0.004144462475753;
            h250[i][IMAG] =   0;
            break;
        case 138:
            h250[i][REAL] = 0.004451760104474;
            h250[i][IMAG] =   0;
            break;
        case 137:
            h250[i][REAL] = 0.004753364509754;
            h250[i][IMAG] =   0;
            break;
        case 136:
            h250[i][REAL] = 0.005048780517428;
            h250[i][IMAG] =   0;
            break;
       case 135:
            h250[i][REAL] = 0.005337522201614;
            h250[i][IMAG] =   0;
            break;
        case 134:
            h250[i][REAL] = 0.005619113769061;
            h250[i][IMAG] =  0;
            break;
        case 133:
            h250[i][REAL] = 0.005893090426249;
            h250[i][IMAG] =  0;
            break;
        case 132:
            h250[i][REAL] = 0.006158999227546;
            h250[i][IMAG] =  0;
            break;
        case 131:
            h250[i][REAL] = 0.006416399902769;
            h250[i][IMAG] =  0;
            break;
        case 130:
            h250[i][REAL] = 0.006664865662513;
            h250[i][IMAG] =  0;
            break;
        case 129:
            h250[i][REAL] = 0.006903983979682;
            h250[i][IMAG] =  0;
            break;
        case 128:
            h250[i][REAL] = 0.007133357345681;
            h250[i][IMAG] =  0;
            break;
        case 127:
            h250[i][REAL] = 0.007352603999776;
            h250[i][IMAG] =  0;
            break;
        case 126:
            h250[i][REAL] = 0.007561358630206;
            h250[i][IMAG] =  0;
            break;
        case 125:
            h250[i][REAL] = 0.007759273045649;
            h250[i][IMAG] =  0;
            break;
        case 124:
            h250[i][REAL] = 0.007946016815726;
            h250[i][IMAG] =  0;
            break;
        case 123:
            h250[i][REAL] = 0.008121277879292;
            h250[i][IMAG] =   0;
            break;
        case 122:
            h250[i][REAL] = 0.008284763119299;
            h250[i][IMAG] =   0;
            break;
        case 121:
            h250[i][REAL] = 0.008436198903098;
            h250[i][IMAG] =  0;
            break;
        case 120:
            h250[i][REAL] = 0.008575331587137;
            h250[i][IMAG] =   0;
            break;
        case 119:
            h250[i][REAL] = 0.008701927985023;
            h250[i][IMAG] =   0;
            break;
        case 118:
            h250[i][REAL] = 0.008815775798054;
            h250[i][IMAG] =  0;
            break;
        case 117:
            h250[i][REAL] = 0.008916684007356;
            h250[i][IMAG] =  0;
            break;
        case 116:
            h250[i][REAL] = 0.009004483226852;
            h250[i][IMAG] =  0;
            break;
        case 115:
            h250[i][REAL] = 0.009079026016366;
            h250[i][IMAG] =  0;
            break;
        case 114:
            h250[i][REAL] = 0.009140187154248;
            h250[i][IMAG] =   0;
            break;
        case 113:
            h250[i][REAL] = 0.009187863868967;
            h250[i][IMAG] =   0;
            break;
        case 112:
            h250[i][REAL] = 0.009221976029235;
            h250[i][IMAG] =  0;
            break;
        case 111:
            h250[i][REAL] = 0.009242466292264;
            h250[i][IMAG] =   0;
            break;
        default:
            h250[i][REAL] =  0;
            h250[i][IMAG] =  0;
            break;
        }
    }

    fft(h250,H250); //calcula H[k]

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<(219)){
                x_n[i][REAL]= sol250[i];
                x_n[i][IMAG]= 0;
            }

        else {
            if(i<(1024 + 219)){
                    in_250[i-219]=static_cast<double>(in[i- (199+20)]);
                    x_n[i][REAL] = in_250[i-219];
                    x_n[i][IMAG] = 0;
               }

            else {
                x_n[i][REAL] = 0;
                x_n[i][IMAG] = 0;
            }
        }
    }

    for(int i=0;i < 219;i++){  // rellena con la entrada actual
        sol250[i] = in_250[(1024 - 219) + i];
    }

    fft(x_n,X); //calcula X[k]

    for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
        Y[i][REAL]=H250[i][REAL]*X[i][REAL]-H250[i][IMAG]*X[i][IMAG];
        Y[i][IMAG]=H250[i][REAL]*X[i][IMAG]+H250[i][IMAG]*X[i][REAL];
    }

    idft(Y,y);        //calcula y[n]

    for(int n=0;n<1024;++n){
        out_250[n]= y[n + 219][REAL];
        out_250[n]=(0.02)*(volumeGain)*(out_250[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
        out[n]=static_cast<float>(out_250[n]);// se hace conversion de double a float
    }

//   energia250=FFT(blockSize,out_250);//se determina la energia de la banda

    delete out_250;
    delete in_250;
    fftw_free(x_n);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h250);
    fftw_free(H250);
    */

    for(int i =0; i<blockSize; i++){
        out[i] = 0;
    }

}



//-------------------------------------------------FILTRO DE 125Hz------------------------------------------------------------------------------//
void controlVolume::filter_125(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 125Hz
    double s_125=0.011733114619205466;
    double a_0_125=-1.9904603982304521;
    double a_1_125=0.99105006098634807;
    double b_0_125=-0.88514816518196571;
    double a_2_125=-1.9950728402576143;
    double a_3_125=0.99523931323078052;
    double b_1_125=-1.9999999360407035;

    double *tmpout_5=new double[blockSize];
    double *tmpout_6=new double[blockSize];
    double *s=new double[blockSize];
    double *out_2=new double[blockSize];



    if(inicial){
        _debug("filtro de 125Hz---------------------------------------------------------------------------------------------------" << std::endl);
        for(int n=0;n<blockSize;++n){

            if(n==0){
                tmpout_5[n]=(s_125)*in[n];
                s[n]=tmpout_5[n];
                tmpout_6[n]=(s_125)*s[n];
                out_2[n]=tmpout_5[n];
            }
            else if(n==1){
                tmpout_5[n]=(s_125)*in[n]-(a_0_125)*tmpout_5[n-1];
                s[n]=tmpout_5[n]+(b_0_125)*tmpout_5[n-1];
                tmpout_6[n]=(s_125)*s[n]-(a_2_125)*tmpout_6[n-1];
                out_2[n]=tmpout_6[n]+(b_1_125)*tmpout_6[n-1];
            }
            else{
                tmpout_5[n]=(s_125)*in[n]-(a_0_125)*tmpout_5[n-1]-(a_1_125)*tmpout_5[n-2];
                s[n]=tmpout_5[n]+(b_0_125)*tmpout_5[n-1]+tmpout_5[n-2];
                tmpout_6[n]=(s_125)*s[n]-(a_2_125)*tmpout_6[n-1]-(a_3_125)*tmpout_6[n-2];
                out_2[n]=tmpout_6[n]+(b_1_125)*tmpout_6[n-1]+tmpout_6[n-2];
            }
            out_2[n]=(0.02)*(volumeGain)*out_2[n];
            out[n]=static_cast<float>(out_2[n]);
        }

     }
    else{

            for(int n=0;n<blockSize;++n){

                if(n==0){
                    tmpout_5[n]=(s_125)*in[n]-(a_0_125)*f-(a_1_125)*e;
                    s[n]=tmpout_5[n]+(b_0_125)*f+e;
                    tmpout_6[n]=(s_125)*s[n]-(a_2_125)*h-(a_3_125)*g;
                    out_2[n]=tmpout_6[n]+(b_1_125)*h+g;
                }
                else if(n==1){
                    tmpout_5[n]=(s_125)*in[n]-(a_0_125)*tmpout_5[n-1]-(a_1_125)*f;
                    s[n]=tmpout_5[n]+(b_0_125)*tmpout_5[n-1]+f;
                    tmpout_6[n]=(s_125)*s[n]-(a_2_125)*tmpout_6[n-1]-(a_3_125)*h;
                    out_2[n]=tmpout_6[n]+(b_1_125)*tmpout_6[n-1]+h;
                }
                else{
                    tmpout_5[n]=(s_125)*in[n]-(a_0_125)*tmpout_5[n-1]-(a_1_125)*tmpout_5[n-2];
                    s[n]=tmpout_5[n]+(b_0_125)*tmpout_5[n-1]+tmpout_5[n-2];
                    tmpout_6[n]=(s_125)*s[n]-(a_2_125)*tmpout_6[n-1]-(a_3_125)*tmpout_6[n-2];
                    out_2[n]=tmpout_6[n]+(b_1_125)*tmpout_6[n-1]+tmpout_6[n-2];
                }
                out_2[n]=(0.02)*(volumeGain)*out_2[n];
                out[n]=static_cast<float>(out_2[n]);
            }
       }

      energia125=FFT(blockSize,out_2);
      e=tmpout_5[1022];//w1(-2)
      f=tmpout_5[1023];//w1(-1)
      g=tmpout_6[1022];//w2(-2)
      h=tmpout_6[1023];//w2(-1)

      delete tmpout_5;
      delete tmpout_6;
      delete s;
      delete out_2;

}
//-------------------------------------------------FILTRO DE 63Hz------------------------------------------------------------------------------//
void controlVolume::filter_63(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 63Hz

    double s_63=0.01045358600130269;
    double a_0_63=-1.995367098268765;
    double a_1_63=0.99551484954682889;
    double b_0_63=-1.6476735248656746;
    double a_2_63=-1.9975751976645493;
    double a_3_63=0.99761686592156373;
    double b_1_63=-1.9999999840095546;
    double *tmpout_7=new double[blockSize]; //temporal out de filtro 1KHz
    double *tmpout_8=new double[blockSize];
    double *w=new double[blockSize];
    double *out_3=new double[blockSize];



    if(inicial){
        _debug("filtro de 63Hz---------------------------------------------------------------------------------------------------" << std::endl);
        for(int n=0;n<blockSize;++n){

            if(n==0){
                tmpout_7[n]=(s_63)*in[n];
                w[n]=tmpout_7[n];
                tmpout_8[n]=(s_63)*w[n];
                out_3[n]=tmpout_8[n];
            }
            else if(n==1){
                tmpout_7[n]=(s_63)*in[n]-(a_0_63)*tmpout_7[n-1];
                w[n]=tmpout_7[n]+(b_0_63)*tmpout_7[n-1];
                tmpout_8[n]=(s_63)*w[n]-(a_2_63)*tmpout_8[n-1];
                out_3[n]=tmpout_8[n]+(b_1_63)*tmpout_8[n-1];
            }
            else{
                tmpout_7[n]=(s_63)*in[n]-(a_0_63)*tmpout_7[n-1]-(a_1_63)*tmpout_7[n-2];
                w[n]=tmpout_7[n]+(b_0_63)*tmpout_7[n-1]+tmpout_7[n-2];
                tmpout_8[n]=(s_63)*w[n]-(a_2_63)*tmpout_8[n-1]-(a_3_63)*tmpout_8[n-2];
                out_3[n]=tmpout_8[n]+(b_1_63)*tmpout_8[n-1]+tmpout_8[n-2];
            }
            out_3[n]=(0.02)*(volumeGain)*out_3[n];
            out[n]=static_cast<float>(out_3[n]);
        }

     }
    else{

            for(int n=0;n<blockSize;++n){

                if(n==0){
                    tmpout_7[n]=(s_63)*in[n]-(a_0_63)*s-(a_1_63)*r;
                    w[n]=tmpout_7[n]+(b_0_63)*s+r;
                    tmpout_8[n]=(s_63)*w[n]-(a_2_63)*u-(a_3_63)*t;
                    out_3[n]=tmpout_8[n]+(b_1_63)*u+t;
                }
                else if(n==1){
                    tmpout_7[n]=(s_63)*in[n]-(a_0_63)*tmpout_7[n-1]-(a_1_63)*s;
                    w[n]=tmpout_7[n]+(b_0_63)*tmpout_7[n-1]+s;
                    tmpout_8[n]=(s_63)*w[n]-(a_2_63)*tmpout_8[n-1]-(a_3_63)*u;
                    out_3[n]=tmpout_8[n]+(b_1_63)*tmpout_8[n-1]+u;
                }
                else{
                    tmpout_7[n]=(s_63)*in[n]-(a_0_63)*tmpout_7[n-1]-(a_1_63)*tmpout_7[n-2];
                    w[n]=tmpout_7[n]+(b_0_63)*tmpout_7[n-1]+tmpout_7[n-2];
                    tmpout_8[n]=(s_63)*w[n]-(a_2_63)*tmpout_8[n-1]-(a_3_63)*tmpout_8[n-2];
                    out_3[n]=tmpout_8[n]+(b_1_63)*tmpout_8[n-1]+tmpout_8[n-2];
                }
                out_3[n]=(0.02)*(volumeGain)*out_3[n];
                out[n]=static_cast<float>(out_3[n]);
            }
       }

      energia64=FFT(blockSize,out_3);
      r=tmpout_7[1022];//w1(-2)//e
      s=tmpout_7[1023];//w1(-1)//f
      t=tmpout_8[1022];//w2(-2)//g
      u=tmpout_8[1023];//w2(-1)//h

      delete tmpout_7;
      delete tmpout_8;
      delete w;
      delete out_3;
}
//--------------------------------------------------FILTRO DE 31.5Hz---------------------------------------------------------------------------//
void controlVolume::filter_31_5(int blockSize, int volumeGain, bool inicial, float *in, float *out){
    double s_31=0.01011125550119148;
    double a_0_31=-1.9977179089948256;
    double a_1_31=0.99775488857850614;
    double b_0_31=-1.9056911829229985;
    double a_2_31=-1.9987973055584447;
    double a_3_31=0.99880772887378833;
    double b_1_31=-1.9999999960023487;
    double *tmpout_i=new double[blockSize]; //temporal out de filtro 1KHz
    double *tmpout_j=new double[blockSize];
    double *uu=new double[blockSize];
    double *out_0=new double[blockSize];



    if(inicial){
        _debug("filtro de 31.5Hz---------------------------------------------------------------------------------------------------" << std::endl);
        for(int n=0;n<blockSize;++n){

            if(n==0){
                tmpout_i[n]=(s_31)*in[n];
                uu[n]=tmpout_i[n];
                tmpout_j[n]=(s_31)*uu[n];
                out_0[n]=tmpout_j[n];
            }
            else if(n==1){
                tmpout_i[n]=(s_31)*in[n]-(a_0_31)*tmpout_i[n-1];
                uu[n]=tmpout_i[n]+(b_0_31)*tmpout_i[n-1];
                tmpout_j[n]=(s_31)*uu[n]-(a_2_31)*tmpout_j[n-1];
                out_0[n]=tmpout_j[n]+(b_1_31)*tmpout_j[n-1];
            }
            else{
                tmpout_i[n]=(s_31)*in[n]-(a_0_31)*tmpout_i[n-1]-(a_1_31)*tmpout_i[n-2];
                uu[n]=tmpout_i[n]+(b_0_31)*tmpout_i[n-1]+tmpout_i[n-2];
                tmpout_j[n]=(s_31)*uu[n]-(a_2_31)*tmpout_j[n-1]-(a_3_31)*tmpout_j[n-2];
                out_0[n]=tmpout_j[n]+(b_1_31)*tmpout_j[n-1]+tmpout_j[n-2];
            }
            out_0[n]=(0.02)*(volumeGain)*out_0[n];
            out[n]=static_cast<float>(out_0[n]);
        }

     }
    else{

            for(int n=0;n<blockSize;++n){

                if(n==0){
                    tmpout_i[n]=(s_31)*in[n]-(a_0_31)*eb-(a_1_31)*ea;
                    uu[n]=tmpout_i[n]+(b_0_31)*eb+ea;
                    tmpout_j[n]=(s_31)*uu[n]-(a_2_31)*ed-(a_3_31)*ec;
                    out_0[n]=tmpout_j[n]+(b_1_31)*ed+ec;
                }
                else if(n==1){
                    tmpout_i[n]=(s_31)*in[n]-(a_0_31)*tmpout_i[n-1]-(a_1_31)*eb;
                    uu[n]=tmpout_i[n]+(b_0_31)*tmpout_i[n-1]+eb;
                    tmpout_j[n]=(s_31)*uu[n]-(a_2_31)*tmpout_j[n-1]-(a_3_31)*ed;
                    out_0[n]=tmpout_j[n]+(b_1_31)*tmpout_j[n-1]+ed;
                }
                else{
                    tmpout_i[n]=(s_31)*in[n]-(a_0_31)*tmpout_i[n-1]-(a_1_31)*tmpout_i[n-2];
                    uu[n]=tmpout_i[n]+(b_0_31)*tmpout_i[n-1]+tmpout_i[n-2];
                    tmpout_j[n]=(s_31)*uu[n]-(a_2_31)*tmpout_j[n-1]-(a_3_31)*tmpout_j[n-2];
                    out_0[n]=tmpout_j[n]+(b_1_31)*tmpout_j[n-1]+tmpout_j[n-2];
                }

                out_0[n]=(0.02)*(volumeGain)*out_0[n];
                out[n]=static_cast<float>(out_0[n]);
            }
       }

      energia31=FFT(blockSize,out_0);
      ea=tmpout_i[1022];//w1(-2)//m
      eb=tmpout_i[1023];//w1(-1)//o
      ec=tmpout_j[1022];//w2(-2)//p
      ed=tmpout_j[1023];//w2(-1)//q

      delete tmpout_i;
      delete tmpout_j;
      delete uu;
      delete out_0;
}

/**
 * Ecualizador Digital
 * Parametros de entrada:
 * tama;o de bloque
 * ganancia de volumen
 * ganancia de slider 1:31.5 Hz
 * ganancia de slider 2: 63 Hz
 * ganancia de slider 3: 125 Hz
 * ganancia de slider 4: 250 Hz
 * ganancia de slider 5: 500 Hz
 * ganancia de slider 6: 1kHz
 * ganancia de slider 7: 2kHz
 * ganancia de slider 8: 4kHz
 * ganancia de slider 9: 8kHz
 * ganancia de slider 10: 16 kHz
 * condicion de inicio de primer frame
 * puntero de entrada direccionado a datos de cancion
 * punetro de salida de DSP aplicado
 */
//-------------------------------------------------Ecualizador Final------------------------------------------------------------------------------//

void controlVolume::eq(int blockSize, int volumeGain, int g1, int g2, int g3, int g4, int g5, int g6, int g7, int g8, int g9,int g10, bool inicial, float* in, float* out){
    // se definen punterios temporales de tipo float para la salida de cada filtro
    float *tmp= new float[blockSize];
    float *tmp1=new float[blockSize];
    float *tmp2=new float[blockSize];
    float *tmp3=new float[blockSize];
    float *tmp4=new float[blockSize];
    float *tmp5=new float[blockSize];
    float *tmp6=new float[blockSize];
    float *tmp7=new float[blockSize];
    float *tmp8=new float[blockSize];
    float *tmp9=new float[blockSize];

    // se llama a cada filtro y se pasan los parametroa de ganancia, tama;no de bloque, entrada y salida
    // conexion en paralelo de filtros de ecualizador de 10 bandas
    filter(blockSize,g10,inicial,in,tmp9);

    filter_8k(blockSize,g9,inicial,in,tmp8);

    filter_4k(blockSize,g8,inicial,in,tmp7);

    filter_2k(blockSize,g7,inicial,in,tmp6);

    filter_1k(blockSize,g6,inicial,in,tmp5);

    filter_500(blockSize,g5,inicial,in,tmp4);

    filter_250(blockSize,g4,inicial,in,tmp3);

    filter_125(blockSize,g3,inicial,in,tmp2);

    filter_63(blockSize,g2,inicial,in,tmp1);

    filter_31_5(blockSize,g1,inicial,in,tmp);



    // se suma el resultado de cada una de las salidas de los filtros, nodo principal de suma de contribucion de cada banda

    for(int n=0;n<blockSize;++n){
        out[n]=(0.038)*(volumeGain)*(tmp[n]+tmp1[n]+tmp2[n]+tmp3[n]+tmp4[n]+tmp5[n]+tmp6[n]+tmp7[n]+tmp9[n]+tmp8[n]);

    }


// se eliminan los punteros temporales para evitar conflictos de memoria

    delete tmp;
    delete tmp1;
    delete tmp2;
    delete tmp3;
    delete tmp4;
    delete tmp5;
    delete tmp6;
    delete tmp7;
    delete tmp8;
    delete tmp9;

}

//-------------------------------------------------------------------------------------------------------------Funciones de DFT-------------------------------------------------

/**
 * controlVolume::FFT: funcion que utiliza la biblioteca de FFT para determinar la DFT de cada senal
 * Se utiliza la relacion de parceval para calcular la energia
 * blockSize_: tamano de la muestra y tamano de la DFT
 * input: puntero a la senal de entrada
 * Se retorna el valor de la energia
 */

int controlVolume::FFT(int blockSize_,double *input){

    int N;     // se define el tamano de la DFT
    N=blockSize_;
    double outtotal=0;  // variable para calcular la energia



    fftw_complex *out;
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan my_plan; // se crea variable para almacenar el metodo a ejecutar de la biblbioteca
    my_plan =fftw_plan_dft_r2c_1d(N, input, out,FFTW_ESTIMATE); // se dan a tributos al metodo a utilizar
    fftw_execute(my_plan); // se ejecuta la FFT, se ejecuta el metodo


    for(int i=0;i<(N/2)+1;i++) // Se recorre el arreglo y se calcula la energia
    {
        outtotal=outtotal+((pow(out[i][0],2)+pow(out[i][1],2))/N);    // Relacion  de Parseval para obtener la energia de cada banda.
    }

    outtotal=10*(2*outtotal-((pow(out[N/2][0],2)+pow(out[N/2][1],2))/N)); // se da una ganancia a la energia para que se visualice correctamente en las barras

    fftw_destroy_plan(my_plan); // se elimina el m[etodo utilizado
    fftw_free(out); // se libera puntero de salida de tipo complejo




    //_debug("Llegue aqui-------------------------------------------------------"<<outtotal << std::endl);
    return static_cast<int>(outtotal); // se retorna como numero entero el resultado de Parceval

}

//---------------------------------------------------------------------------------------------------------------------------------------------
/**
 * controlVolume::valorbarraN_: las siguientes funciones devuelven el valor de la energia de cada banda
 * dentro de cada funcion se implenta un circuito RC para retardar el movimiento de las barras
 */

int controlVolume::valorbarra1_(){      // funcion para devolver el valor de la energia como un entero
    int salida=0;
    salida=energia31+0.5*y1;     // implementacion de circuito RC
    y1=salida;
    return salida;
}

int controlVolume::valorbarra2_(){
    int salida=0;
    salida=energia64+0.5*y2;
    y2=salida;
    return salida;

}

int controlVolume::valorbarra3_(){
    int salida=0;
    salida=energia125+0.5*y3;
    y3=salida;
    return salida;
}

int controlVolume::valorbarra4_(){

    int salida=0;
    salida=energia250+0.5*y4;
    y4=salida;
    return salida;
}

int controlVolume::valorbarra5_(){
    int salida=0;
    salida=energia500+0.5*y5;
    y5=salida;
    return salida;
}

int controlVolume::valorbarra6_(){
    int salida=0;
    salida=energia1000+0.5*y6;
    y6=salida;
    return salida;
}
int controlVolume::valorbarra7_(){
    int salida=0;
    salida=energia2000+0.5*y7;
    y7=salida;
    return salida;
}

int controlVolume::valorbarra8_(){
    int salida=0;
    salida=energia4000+0.5*y8;
    y8=salida;
    return salida;
}

int controlVolume::valorbarra9_(){
    int salida=0;
    salida=energia8000+0.5*y9;
    y9=salida;
    return salida;
}

int controlVolume::valorbarra10_(){
    int salida=0;
    salida=energia16000+0.5*y10;
    y10=salida;
    return salida;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
