/*
 * DSP Example is part of the DSP Lecture at TEC-Costa Rica
 * Copyright (C) 2010  Pablo Alvarado
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



#include "dspsystem.h"
#include <cstring>

#undef _DSP_DEBUG
#define _DSP_DEBUG

#ifdef _DSP_DEBUG
#define _debug(x) std::cerr << x
#include <iostream>
#else
#define _debug(x)
#endif

/**
 * constructor dspSystem::dspSystem
 */
dspSystem::dspSystem()
  :sampleRate_(0),bufferSize_(0),cv_(0){       // se inician los valores de la libreria dspsystem
}

/**
 * destructor dspSystem::~dspSystem
 */
dspSystem::~dspSystem() {
    delete cv_;     // -------------------------cv_ significa control de volumen

}





//=======================================================================================================================actualiza valores de volumen=============================
void dspSystem::updateVolume(int value){        // ------------------------se modifica la variable de volumen general
   /*
    * Updating volume value
    */
   volumeGain_=value;

}
void dspSystem::updateg1(int value){        // ------------------------se modifica la variable de volumen banda 1
   /*
    * Updating volume value
    */
   g1_=value;

}

void dspSystem::updateg2(int value){        // ------------------------se modifica la variable de volumen banda 2
   /*
    * Updating volume value
    */
   g2_=value;

}


void dspSystem::updateg3(int value){        // ------------------------se modifica la variable de volumen banda 3
   /*
    * Updating volume value
    */
   g3_=value;

}

void dspSystem::updateg4(int value){        // ------------------------se modifica la variable de volumen banda 4
   /*
    * Updating volume value
    */
   g4_=value;

}

void dspSystem::updateg5(int value){        // ------------------------se modifica la variable de volumen banda 5
   /*
    * Updating volume value
    */
   g5_=value;

}

void dspSystem::updateg6(int value){        // ------------------------se modifica la variable de volumen banda 6
   /*
    * Updating volume value
    */
   g6_=value;

}

void dspSystem::updateg7(int value){        // ------------------------se modifica la variable de volumen banda 7
   /*
    * Updating volume value
    */
   g7_=value;

}
void dspSystem::updateg8(int value){        // ------------------------se modifica la variable de volumen banda 8
   /*
    * Updating volume value
    */
   g8_=value;

}
void dspSystem::updateg9(int value){        // ------------------------se modifica la variable de volumen banda 9
   /*
    * Updating volume value
    */
   g9_=value;

}
void dspSystem::updateg10(int value){        // ------------------------se modifica la variable de volumen banda 10
   /*
    * Updating volume value
    */
   g10_=value;

}
//=======================================================================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * bandera para saber cuando inicio una cancion :dspSystem::updateinicio
 * saber: valor de tipo booleano que indica en true el inicio de la cancion
 */

void dspSystem::updateinicio(bool saber){

    iniciorep=saber;

    _debug("cambio el inicio de la cancion" << std::endl);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------



/**
 * Initialization function for the current filter plan
 */
bool dspSystem::init(const int sampleRate,const int bufferSize) {
  _debug("dspSystem::init()" << std::endl);

  sampleRate_ = sampleRate;
  bufferSize_ = bufferSize;
  volumeGain_ = 0;// variables declaras en dspsystem.h, cuando se inicia dspsystem se da valor de 0 al volumen
  g1_=25;g2_=25;g3_=25;g4_=25;g5_=25;g6_=25;g7_=25;g8_=25;g9_=25;g10_=25;
  iniciorep=false;
  delete cv_;
  cv_=new controlVolume();             // se borra puntero y se pide a dios sistema operativo un nuevo puntro de tipo controlVolumen
  return true;
}


/**
 * Processing function
 */
bool dspSystem::process(float* in,float* out) {

  float* tmpIn = in;    // se crean dos elementos temporales punteros tipo flotantes para procesar los datos
  float* tmpOut = out;

  if (iniciorep){
      cv_->eq(bufferSize_,volumeGain_,g1_,g2_,g3_,g4_,g5_,g6_,g7_,g8_,g9_,g10_,iniciorep,tmpIn,tmpOut); // se llama al puntero de tipo controlVolumen y se llama a la funcion filtro.
      //cv_->filter(bufferSize_,g6_,iniciorep,tmpIn,tmpOut);
      iniciorep=false;
  }
  else{
      cv_->eq(bufferSize_,volumeGain_,g1_,g2_,g3_,g4_,g5_,g6_,g7_,g8_,g9_,g10_,iniciorep,tmpIn,tmpOut);// una vez identificado el primer frame se procede a esta linea

  }
  return true;
}


// ------------------------------------------------------------------------------para actualizar las barras----------------------------------
/**
 *dspSystem::valorbarraN, retorna a la ventana principal el valor de la energia
 *
 */

int dspSystem::valorbarra1(){

    return cv_->valorbarra1_();
}

int dspSystem::valorbarra2(){
    return cv_->valorbarra2_();
}

int dspSystem::valorbarra3(){
    return cv_->valorbarra3_();
}

int dspSystem::valorbarra4(){
    return cv_->valorbarra4_();
}

int dspSystem::valorbarra5(){
    return cv_->valorbarra5_();
}

int dspSystem::valorbarra6(){
    return cv_->valorbarra6_();
}

int dspSystem::valorbarra7(){
    return cv_->valorbarra7_();
}

int dspSystem::valorbarra8(){
    return cv_->valorbarra8_();
}

int dspSystem::valorbarra9(){
    return cv_->valorbarra9_();
}

int dspSystem::valorbarra10(){
    return cv_->valorbarra10_();
}
//--------------------------------------------------------------------------------------------------------------------------------------------





/**
 * Shutdown the processor
 */
bool dspSystem::shutdown() {
  return true;
}

/**
 * Set buffer size (call-back)
 */
int dspSystem::setBufferSize(const int bufferSize) {
  bufferSize_=bufferSize;
  return 1;
}

/**
 * Set sample rate (call-back)
 */
int dspSystem::setSampleRate(const int sampleRate) {
  sampleRate_=sampleRate;
  return 1;
}
