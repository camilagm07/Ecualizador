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


void hn(double K,double a1,double a2,double a3,double a4,double a5,double a6, double a7, double a8, double b1,double b2,double b3,double b4,double b5,double b6, double b7, double b8, double H[2048][2]){


    float y_1=0, y_2=0, y_3=0, y_4=0, y_5=0, y_6=0, y_7=0, y_8=0; //Init Condiciones iniciales

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
            case 7:
                H[i][REAL]= K*b6 - a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4 - a5*y_5 - a6*y_6 - a7*y_7; //Idem
                break;
            case 8:
                H[i][REAL]= K*b6 - a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4 - a5*y_5 - a6*y_6 - a7*y_7 - a8*y_8; //Idem
                break;
            default:
                H[i][REAL]= -a1*y_1 - a2*y_2 - a3*y_3 - a4*y_4 - a5*y_5 - a6*y_6 - a7*y_7 - a8*y_8; // n > 6
                break;
            }

            y_8 = y_7;
            y_7 = y_6;
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

    //****************************************************************************
}


/*
 * Destructor
 */

controlVolume::~controlVolume(){

}



//-------------------------------------------------FILTRO DE 16KHz de orden 6----------------------------------------------------------------///

/*
 * Filtro de 16 kHz de orden 2
 * Parametros de entrada:
 * tama;o de bloque
 * ganancia
 * condicion de primer frame
 * puntero entrada
 * puntero salida
 */
void controlVolume::filter_16k(int blockSize, int volumeGain, bool inicial, float *in, float *out){

    int tamano = 2048;

    fftw_complex *x;
    x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *h10;
    h10 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *H10;
    H10 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    double *out_10 = new double[blockSize];
    double *in_10  = new double[blockSize];

    if (haux){    //calcula el h[n] una sola vez a partir de los coeficientes de la ecuacion de diferencias
        sol[300] = {0};
        haux=false;
    }

    double K=0.07875912;
    double a1=3.44387319;
    double a2=5.29844359;
    double a3=5.13528975;
    double a4=3.47957963;
    double a5=1.47740514;
    double a6=0.28411879;
    double a7=0;
    double a8=0;
    double b1=0.04060383;
    double b2=-2.91809914;
    double b3=0;
    double b4=2.91809914;
    double b5=-0.04060383;
    double b6=-1;
    double b7=0;
    double b8=0;


    hn(K,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,h10);

    fft(h10,H10); //calcula H[k]



    // SOLAPAMIENTO

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<300){
                x[i][REAL]= sol[i];
                x[i][IMAG]= 0;
                sol[i] = in_10[724 + i];
            }

        else {
            if(i<1324){
                    in_10[i-300]=static_cast<double>(in[i-300]);
                    x[i][REAL] = in_10[i-300];
                    x[i][IMAG] = 0;
               }

            else {
                x[i][REAL] = 0;
                x[i][IMAG] = 0;
            }
        }
    }

    for(int i=0;i<300;i++){  // rellena con la entrada actual
        sol[i] = in_10[724 + i];
    }

    fft(x,X); //calcula X[k]

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

    energia16000=FFT(blockSize,out_10);//se determina la energia de la banda

    delete out_10;
    delete in_10;
    fftw_free(x);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h10);
    fftw_free(H10);
}


//-------------------------------------------------FILTRO DE 8KHz------------------------------------------------------------------------------//

/*
 * Filtro de 8 kHz de orden 4
 * Parametros de entrada:
 * tama;o de bloque
 * ganancia
 * condicion de primer frame
 * puntero entrada
 * puntero salida
 */

void controlVolume::filter_8k(int blockSize, int volumeGain, bool inicial, float *in, float *out){

        int tamano = 2048;
        fftw_complex *x;
        x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

        fftw_complex *X;
        X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

        fftw_complex *y;
        y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

        fftw_complex *Y;
        Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

        fftw_complex *h8;
        h8 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

        fftw_complex *H8;
        H8 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

        double *out_8 = new double[blockSize];
        double *in_8  = new double[blockSize];

        double K = 0.0305879973691657;
        double a1 = -1.7904118543858;
        double a2 = 2.77307718603103;
        double a3 = -2.48561443893466;
        double a4 = 2.04168398689785;
        double a5 = -0.910812362274974;
        double a6 = 0.36902133643413;
        double a7 = 0;
        double a8 = 0;
        double b1 = -0.0418178919003749;
        double b2 = -2.87523422267106;
        double b3 = 0;
        double b4 = 2.87523422267106;
        double b5 = 0.0418178919003749;
        double b6 = -1;
        double b7 = 0;
        double b8 = 0;

        hn(K,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,h8);
        fft(h8,H8); //calcula H[k]

        for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
            if (i<300){
                    x[i][REAL]= sol[i];
                    x[i][IMAG]= 0;
                }

            else {
                if(i<1324){
                        in_8[i-300]=static_cast<double>(in[i-300]);
                        x[i][REAL] = in_8[i-300];
                        x[i][IMAG] = 0;
                   }

                else {
                    x[i][REAL] = 0;
                    x[i][IMAG] = 0;
                }
            }
        }

        fft(x,X); //calcula X[k]

        for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
            Y[i][REAL]=H8[i][REAL]*X[i][REAL]-H8[i][IMAG]*X[i][IMAG];
            Y[i][IMAG]=H8[i][REAL]*X[i][IMAG]+H8[i][IMAG]*X[i][REAL];
        }

        idft(Y,y);        //calcula y[n]

        for(int n=0;n<1024;++n){
            out_8[n]= y[n + 300][REAL];
            out_8[n]=(0.02)*(volumeGain)*(out_8[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
            out[n]=static_cast<float>(out_8[n]);// se hace conversion de double a float
        }

//        energia8000=FFT(blockSize,out_8);//se determina la energia de la banda

        delete out_8;
        delete in_8;
        fftw_free(x);
        fftw_free(X);
        fftw_free(Y);
        fftw_free(h8);
        fftw_free(H8);
}

//-------------------------------------------------FILTRO DE 4KHz------------------------------------------------------------------------------//
void controlVolume::filter_4k(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 2kHz

    int tamano = 2048;
    fftw_complex *x;
    x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *h4;
    h4 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *H4;
    H4 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    double *out_4 = new double[blockSize];
    double *in_4  = new double[blockSize];

    double K = 0.000649931517402985;
    double a1 = -5.11984739057556;
    double a2 = 11.5291867282101;
    double a3 = -14.5094894816439;
    double a4 = 10.7524500972001;
    double a5 = -4.45342292600299;
    double a6 = 0.811437702266877;
    double a7 = 0;
    double a8 = 0;
    double b1 = -1.49111638930048;
    double b2 = 0.00365581508920376;
    double b3 = 0;
    double b4 = -0.00365581508920376;
    double b5 = 1.49111638930048;
    double b6 = -1;
    double b7 = 0;
    double b8 = 0;

    hn(K,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,h4);
    fft(h4,H4); //calcula H[k]

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<300){
                x[i][REAL]= sol[i];
                x[i][IMAG]= 0;
            }

        else {
            if(i<1324){
                    in_4[i-300]=static_cast<double>(in[i-300]);
                    x[i][REAL] = in_4[i-300];
                    x[i][IMAG] = 0;
               }

            else {
                x[i][REAL] = 0;
                x[i][IMAG] = 0;
            }
        }
    }

    fft(x,X); //calcula X[k]

    for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
        Y[i][REAL]=H4[i][REAL]*X[i][REAL]-H4[i][IMAG]*X[i][IMAG];
        Y[i][IMAG]=H4[i][REAL]*X[i][IMAG]+H4[i][IMAG]*X[i][REAL];
    }

    idft(Y,y);        //calcula y[n]

    for(int n=0;n<1024;++n){
        out_4[n]= y[n + 300][REAL];
        out_4[n]=(0.02)*(volumeGain)*(out_4[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
        out[n]=static_cast<float>(out_4[n]);// se hace conversion de double a float
    }

//        energia4000=FFT(blockSize,out_4);//se determina la energia de la banda

    delete out_4;
    delete in_4;
    fftw_free(x);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h4);
    fftw_free(H4);
}

//-------------------------------------------------FILTRO DE 1KHz------------------------------------------------------------------------------//
void controlVolume::filter_1k(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 1kHz

    int tamano = 2048;
    fftw_complex *x;
    x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *X;
    X = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *y;
    y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *Y;
    Y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *h1;
    h1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    fftw_complex *H1;
    H1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);

    double *out_1 = new double[blockSize];
    double *in_1 = new double[blockSize];

    double K = 0.000409262614240481;
    double a1 = -5.73904408499587;
    double a2 = 13.7892758259963;
    double a3 = -17.754602031527;
    double a4 = 12.9204032623488;
    double a5 = -5.0387879593448;
    double a6 = 0.82276253431545;
    double a7 = 0;
    double a8 = 0;
    double b1 = -2.11909798107263;
    double b2 = 1.23841599355237;
    double b3 = 0;
    double b4 = -1.23841599355237;
    double b5 = 2.11909798107263;
    double b6 = -1;
    double b7 = 0;
    double b8 = 0;

    hn(K,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,h1);
    fft(h1,H1); //calcula H[k]

    for(int i=0;i<2048;i++){  // introduce los ultimos 300 del x[n-1]
        if (i<300){
                x[i][REAL]= sol[i];
                x[i][IMAG]= 0;
            }

        else {
            if(i<1324){
                    in_1[i-300]=static_cast<double>(in[i-300]);
                    x[i][REAL] = in_1[i-300];
                    x[i][IMAG] = 0;
               }

            else {
                x[i][REAL] = 0;
                x[i][IMAG] = 0;
            }
        }
    }

    fft(x,X); //calcula X[k]

    for(int i=0;i<2048;i++){  // Y(k)=H[k]X[k]
        Y[i][REAL]=H1[i][REAL]*X[i][REAL]-H1[i][IMAG]*X[i][IMAG];
        Y[i][IMAG]=H1[i][REAL]*X[i][IMAG]+H1[i][IMAG]*X[i][REAL];
    }

    idft(Y,y);        //calcula y[n]

    for(int n=0;n<1024;++n){
        out_1[n]= y[n + 300][REAL];
        out_1[n]=(0.02)*(volumeGain)*(out_1[n]);//filtro de ganancia unitaria en banda pasante, se escala por 0.02 para ajustar la ganancia del slider
        out[n]=static_cast<float>(out_1[n]);// se hace conversion de double a float
    }

//        energia4000=FFT(blockSize,out_4);//se determina la energia de la banda

    delete out_1;
    delete in_1;
    fftw_free(x);
    fftw_free(X);
    fftw_free(Y);
    fftw_free(h1);
    fftw_free(H1);
}

//-------------------------------------------------FILTRO DE 2KHz------------------------------------------------------------------------------//
void controlVolume::filter_2k(int blockSize, int volumeGain, bool inicial, float *in, float *out){//filtro de 2kHz
    double s_2=0.095049796405053844;
    double a_0_2=-1.725609356082455;
    double a_1_2=0.86668619099358524;
    double b_0_2=1.9621261892748729;
    double a_2_2=-1.8843420257509356;
    double a_3_2=0.92568054774238073;
    double b_1_2=-1.9999837077911962;

    double *tmpout_3=new double[blockSize]; //temporal out de filtro 1KHz
    double *tmpout_4=new double[blockSize];
    double *x=new double[blockSize];
    double *out_1=new double[blockSize];



    if(inicial){
        _debug("filtro de 2kHz---------------------------------------------------------------------------------------------------" << std::endl);
        for(int n=0;n<blockSize;++n){

            if(n==0){
                tmpout_3[n]=(s_2)*in[n];
                x[n]=tmpout_3[n];
                tmpout_4[n]=(s_2)*x[n];
                out_1[n]=tmpout_4[n];
            }
            else if(n==1){
                tmpout_3[n]=(s_2)*in[n]-(a_0_2)*tmpout_3[n-1];
                x[n]=tmpout_3[n]+( b_0_2)*tmpout_3[n-1];
                tmpout_4[n]=(s_2)*x[n]-(a_2_2)*tmpout_4[n-1];
                out_1[n]=tmpout_4[n]+(b_1_2)*tmpout_4[n-1];
            }
            else{
                tmpout_3[n]=(s_2)*in[n]-(a_0_2)*tmpout_3[n-1]-(a_1_2)*tmpout_3[n-2];
                x[n]=tmpout_3[n]+( b_0_2)*tmpout_3[n-1]+tmpout_3[n-2];
                tmpout_4[n]=(s_2)*x[n]-(a_2_2)*tmpout_4[n-1]-(a_3_2)*tmpout_4[n-2];
                out_1[n]=tmpout_4[n]+(b_1_2)*tmpout_4[n-1]+tmpout_4[n-2];
            }
            out_1[n]=(0.02)*(volumeGain)*(out_1[n]);
            out[n]=static_cast<float>(out_1[n]);
        }

     }
    else{

            for(int n=0;n<blockSize;++n){

                if(n==0){
                    tmpout_3[n]=(s_2)*in[n]-(a_0_2)*o-(a_1_2)*m;
                    x[n]=tmpout_3[n]+( b_0_2)*o+m;
                    tmpout_4[n]=(s_2)*x[n]-(a_2_2)*q-(a_3_2)*p;
                    out_1[n]=tmpout_4[n]+(b_1_2)*q+p;
                }
                else if(n==1){
                    tmpout_3[n]=(s_2)*in[n]-(a_0_2)*tmpout_3[n-1]-(a_1_2)*o;
                    x[n]=tmpout_3[n]+( b_0_2)*tmpout_3[n-1]+o;
                    tmpout_4[n]=(s_2)*x[n]-(a_2_2)*tmpout_4[n-1]-(a_3_2)*q;
                    out_1[n]=tmpout_4[n]+(b_1_2)*tmpout_4[n-1]+q;
                }
                else{
                    tmpout_3[n]=(s_2)*in[n]-(a_0_2)*tmpout_3[n-1]-(a_1_2)*tmpout_3[n-2];
                    x[n]=tmpout_3[n]+( b_0_2)*tmpout_3[n-1]+tmpout_3[n-2];
                    tmpout_4[n]=(s_2)*x[n]-(a_2_2)*tmpout_4[n-1]-(a_3_2)*tmpout_4[n-2];
                    out_1[n]=tmpout_4[n]+(b_1_2)*tmpout_4[n-1]+tmpout_4[n-2];
                }
                out_1[n]=(0.02)*(volumeGain)*(out_1[n]);
                out[n]=static_cast<float>(out_1[n]);
            }
       }

//      energia2000=FFT(blockSize,out_1);
      m=tmpout_3[1022];//w1(-2)//e
      o=tmpout_3[1023];//w1(-1)//f
      p=tmpout_4[1022];//w2(-2)//g
      q=tmpout_4[1023];//w2(-1)//h

      delete tmpout_3;
      delete tmpout_4;
      delete x;
      delete out_1;

}
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

//      energia250=FFT(blockSize,out_4);
      i=tmpout_9[1022];//w1(-2)
      j=tmpout_9[1023];//w1(-1)
      k=tmpout_0[1022];//w2(-2)
      l=tmpout_0[1023];//w2(-1)

      delete tmpout_9;
      delete tmpout_0;
      delete z;
      delete out_4;

}
//-------------------------------------------------FILTRO DE 500Hz------------------------------------------------------------------------------//
void controlVolume::filter_500(int blockSize, int volumeGain, bool inicial, float *in, float *out){
    double s_500=0.026555215679496418;
    double a_0_500=-1.9551029756830791;
    double a_1_500=0.96448409180774719;
    double b_0_500=1.4498228212665976;
    double a_2_500=-1.9783358378384872;
    double a_3_500=0.98099508660482093;
    double b_1_500=-1.9999989745440212;
    double *tmpout_a=new double[blockSize];
    double *tmpout_b=new double[blockSize];
    double *y=new double[blockSize];
    double *out_6=new double[blockSize];



    if(inicial){
        _debug("filtro de 500Hz---------------------------------------------------------------------------------------------------" << std::endl);
        for(int n=0;n<blockSize;++n){

            if(n==0){
                tmpout_a[n]=(s_500)*in[n];
                y[n]=tmpout_a[n];
                tmpout_b[n]=(s_500)*y[n];
                out_6[n]=tmpout_b[n];
            }
            else if(n==1){
                tmpout_a[n]=(s_500)*in[n]-(a_0_500)*tmpout_a[n-1];
                y[n]=tmpout_a[n]+(b_0_500)*tmpout_a[n-1];
                tmpout_b[n]=(s_500)*y[n]-(a_2_500)*tmpout_b[n-1];
                out_6[n]=tmpout_b[n]+(b_1_500)*tmpout_b[n-1];
            }
            else{
                tmpout_a[n]=(s_500)*in[n]-(a_0_500)*tmpout_a[n-1]-(a_1_500)*tmpout_a[n-2];
                y[n]=tmpout_a[n]+(b_0_500)*tmpout_a[n-1]+tmpout_a[n-2];
                tmpout_b[n]=(s_500)*y[n]-(a_2_500)*tmpout_b[n-1]-(a_3_500)*tmpout_b[n-2];
                out_6[n]=tmpout_b[n]+(b_1_500)*tmpout_b[n-1]+tmpout_b[n-2];
            }
           out_6[n]=(0.02)*(volumeGain)*out_6[n];
           out[n]=static_cast<float>(out_6[n]);
        }

     }
    else{

            for(int n=0;n<blockSize;++n){

                if(n==0){
                    tmpout_a[n]=(s_500)*in[n]-(a_0_500)*ab-(a_1_500)*aa;
                    y[n]=tmpout_a[n]+(b_0_500)*ab+aa;
                    tmpout_b[n]=(s_500)*y[n]-(a_2_500)*ad-(a_3_500)*ac;
                    out_6[n]=tmpout_b[n]+(b_1_500)*ad+ac;
                }
                else if(n==1){
                    tmpout_a[n]=(s_500)*in[n]-(a_0_500)*tmpout_a[n-1]-(a_1_500)*ab;
                    y[n]=tmpout_a[n]+(b_0_500)*tmpout_a[n-1]+ab;
                    tmpout_b[n]=(s_500)*y[n]-(a_2_500)*tmpout_b[n-1]-(a_3_500)*ad;
                    out_6[n]=tmpout_b[n]+(b_1_500)*tmpout_b[n-1]+ad;
                }
                else{
                    tmpout_a[n]=(s_500)*in[n]-(a_0_500)*tmpout_a[n-1]-(a_1_500)*tmpout_a[n-2];
                    y[n]=tmpout_a[n]+(b_0_500)*tmpout_a[n-1]+tmpout_a[n-2];
                    tmpout_b[n]=(s_500)*y[n]-(a_2_500)*tmpout_b[n-1]-(a_3_500)*tmpout_b[n-2];
                    out_6[n]=tmpout_b[n]+(b_1_500)*tmpout_b[n-1]+tmpout_b[n-2];
                }
                out_6[n]=(0.02)*(volumeGain)*out_6[n];
                out[n]=static_cast<float>(out_6[n]);
            }
       }

//      energia500=FFT(blockSize,out_6);

      aa=tmpout_a[1022];//w1(-2)//i
      ab=tmpout_a[1023];//w1(-1)//j
      ac=tmpout_b[1022];//w2(-2)//k
      ad=tmpout_b[1023];//w2(-1)//l

      delete tmpout_a;
      delete tmpout_b;
      delete y;
      delete out_6;

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

//      energia125=FFT(blockSize,out_2);
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

 //     energia64=FFT(blockSize,out_3);
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

 //     energia31=FFT(blockSize,out_0);
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
    filter_16k(blockSize,g10,inicial,in,tmp9);

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

    int tamano;     // se define el tamano de la DFT
    tamano=blockSize_;
    double outtotal=0;  // variable para calcular la energia



    fftw_complex *out;
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * tamano);
    fftw_plan my_plan; // se crea variable para almacenar el metodo a ejecutar de la biblbioteca
    my_plan =fftw_plan_dft_r2c_1d(tamano, input, out,FFTW_ESTIMATE); // se dan a tributos al metodo a utilizar
    fftw_execute(my_plan); // se ejecuta la FFT, se ejecuta el metodo


    for(int i=0;i<(tamano/2)+1;i++) // Se recorre el arreglo y se calcula la energia
    {
        outtotal=outtotal+((pow(out[i][0],2)+pow(out[i][1],2))/tamano);    // Relacion  de Parseval para obtener la energia de cada banda.
    }

    outtotal=10*(2*outtotal-((pow(out[tamano/2][0],2)+pow(out[tamano/2][1],2))/tamano)); // se da una ganancia a la energia para que se visualice correctamente en las barras

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




