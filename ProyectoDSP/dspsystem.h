

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

/**
 * \file   freqFilter.h
 *         Implements filtering in the frequency domain
 * \author Pablo Alvarado/Jose Miguel Barboza
 * \date   2010.12.12/2017.05.26
 *
 * $Id: equalizer.cpp $
 */

#ifndef DSPSYSTEM_H
#define DSPSYSTEM_H

#include "processor.h"
#include "jack.h"
#include "controlvolume.h"

class dspSystem : public processor {
public:
  /**
   * Constructor
   */
  dspSystem();

  /**
   * Destructor
   */
  ~dspSystem();

  /**
   * Initialization function for the current filter plan
   */
  virtual bool init(const int frameRate,const int bufferSize);

  /**
   * Processing function
   */
  virtual bool process(float* in,float* out);

  /**
   * Shutdown the processor
   */
  virtual bool shutdown();

  /**
   * Set buffer size
   */
  virtual int setBufferSize(const int bufferSize);

  /**
   * Set frame rate
   */
  virtual int setSampleRate(const int sampleRate);

  void updateVolume(int value);       //---------------------------------se agrego para volumen

  void updateg1(int value);

  void updateg2(int value);

  void updateg3(int value);

  void updateg4(int value);

  void updateg5(int value);

  void updateg6(int value);//500 Hz

  void updateg7(int value);

  void updateg8(int value);

  void updateg9(int value);

  void updateg10(int value);
//---------------------------------------------------------------------------------------------------------------------------------------------
  int valorbarra1();  // para pedir el valores de las barras
  int valorbarra2();
  int valorbarra3();
  int valorbarra4();
  int valorbarra5();
  int valorbarra6();
  int valorbarra7();
  int valorbarra8();
  int valorbarra9();
  int valorbarra10();
  //--------------------------------------------------------------------------------------------------------------------------------------------
  virtual void updateinicio(bool saber);

protected:

  /**
   * Sample rate
   */
  int sampleRate_;

  /**
   * Buffer size
   */
  int bufferSize_;

  /**
   * VolumeGain
   */
//---------------------------------------------------------------------------------------------------------------------------------------------
  bool iniciorep;


  int volumeGain_;//------------------------------------------------------------ variable de control de volumen

  int g1_;  //variable de cada uno de los filtros
  int g2_;
  int g3_;
  int g4_;
  int g5_;
  int g6_;
  int g7_;
  int g8_;
  int g9_;
  int g10_;

  /**
   * control Volume
   */
  controlVolume* cv_;   //------------------------------------------ elemento de tipo control de volumen


};


#endif // DSPSYSTEM_H
