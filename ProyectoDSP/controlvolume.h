
/*
 * DSP Example is part of the DSP Lecture at TEC-Costa Rica
 * Copyright (C) 2017  Jose Miguel Barboza
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
 * \file   controlVolume.h
 *         Implements control volumen in the time domain
 * \author Jose Miguel Barboza
 * \date   2017.31.05
 *
 * $Id: controlVolume.cpp $
 */


#ifndef CONTROLVOLUME_H
#define CONTROLVOLUME_H
#include<fftw3.h>

/**
 * Control Volume class
 *
 * This filter has control in the volume gain of the signal.
 *
 * The filter follows the difference equation
 * \f[
 * y(n)=\cvGain x(n)
 * \f]
 */

class controlVolume {
public:
    /**
     * Constructor
     */
    controlVolume();

    /**
     * Destructor
     */
    ~controlVolume();

   /**
    * Filter the in buffer and leave the result in out
    */

    void eq(int blockSize,
                int volumeGain,int g1,int g2,int g3, int g4,int g5, int g6,int g7,int g8,int g9,int g10,bool inicial,
                float* in,
                float* out);
    void filter_31_5(int blockSize,
                int volumeGain,bool inicial,
                float* in,
                float* out);
    void filter_63(int blockSize,
                int volumeGain,bool inicial,
                float* in,
                float* out);
    void filter_125(int blockSize,
                int volumeGain,bool inicial,
                float* in,
                float* out);
    void filter_250(int blockSize,
                int volumeGain,bool inicial,
                float* in,
                float* out);

    void filter_500(int blockSize,
                int volumeGain,bool inicial,
                float* in,
                float* out);

   void filter_1k(int blockSize,
               int volumeGain,bool inicial,
               float* in,
               float* out);

   void filter_2k(int blockSize,
               int volumeGain,bool inicial,
               float* in,
               float* out);
   void filter_4k(int blockSize,
               int volumeGain,bool inicial,
               float* in,
               float* out);
   void filter_8k(int blockSize,
               int volumeGain,bool inicial,
               float* in,
               float* out);
   void filter_16k(int blockSize,
               int volumeGain,bool inicial,
               float* in,
               float* out);
   void filter(int blockSize,
               int volumeGain,bool inicial,
               float* in,
               float* out);


   //----------------------------------------------parte para DFT
   int ENERGIA(int N, fftw_complex *yk);

   int valorbarra1_();
   int valorbarra2_();
   int valorbarra3_();
   int valorbarra4_();
   int valorbarra5_();
   int valorbarra6_();
   int valorbarra7_();
   int valorbarra8_();
   int valorbarra9_();
   int valorbarra10_();

    int energia31;  // para cambiar los valores de la  barra dft
    int energia64;
    int energia125;
    int energia250;
    int energia500;
    int energia1000;
    int energia2000;
    int energia4000;
    int energia8000;
    int energia16000;
                            //Variables para etapa RC
    int y1;int y2;int y3;int y4;int y5;int y6;int y7;int y8;int y9;int y10;
//---------------------------------------------------------------------------variables
   float a=0;
   float b=0;
   float c=0;
   float d=0;
   double e=0;double f=0;double g=0;double h=0;double i=0;double j=0;double k=0; double l=0;
   double m=0;double q=0; double o=0;double p=0;

   double r=0;double s=0; double t=0; double u=0;
   double aa=0;double ab=0; double ac=0; double ad=0;
   double ba=0;double bb=0; double bc=0; double bd=0;
   double ca;double cb; double cc; double cd;
   double da;double db; double dc; double dd;
   double ea;double eb; double ec; double ed;

   bool haux16 = true;
   double sol16[300];

   bool haux8 = true;
   double sol8[15];

   bool haux4 = true;
   double sol4[30];

   bool haux2 = true;
   double sol2[35];

   bool haux1 = true;
   double sol1[50];

   bool haux500 = true;
   double sol500[100];

   bool haux250 = true;
   double sol250[220];

   bool haux125 = true;
   double sol125[300];

   bool haux64 = true;
   double sol64[350];
};


#endif // CONTROLVOLUME_H
