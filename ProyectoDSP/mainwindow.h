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
 * \file   mainwindow.h
 *         Implements the equalizer H(w) computation
 * \author Pablo Alvarado/Jose Miguel Barboza
 * \date   2010.12.12/2017.05.26
 *
 * $Id: mainwindow.h $
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <QFileDialog>
#include <QtGui>
#include <QtCore>
#include "dspsystem.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

    int volumeGain;


    /**
      *Tolerance value
      */
    static const float Epsilon;

     /**
      * Verbose flag
      */
     bool verbose_;

     /**
      * Timer used to recompute the filter once the user changes the
      * values
      */
     QTimer *timer_;

     /**
      * List of selected files so far
      */
     QStringList selectedFiles_;

     /**
      * Pointer to an inherited class of processor, which does
      * all the real work.
      */
     dspSystem* dsp_;

     /**
      *DSP change
      */
     bool dspChanged_;   //---------------------------------------------variable para saber si ocurrio un cambio de dsp

     bool iniciorep;    //---------------------------------------------variable para saber cuando inicia una cancion
     int gain0;
     int gain1;
     int gain2;
     int gain3;
     int gain4;
     int gain5;
     int gain6;
     int gain7;
     int gain8;
     int gain9;

protected:

     void paintEvent(QPaintEvent *e);

   private slots:
     void on_fileEdit_returnPressed();
     void on_fileButton_clicked();
     void on_volumeSlider_valueChanged(int value);


     void update();

     void on_g1_valueChanged(int value);
     void on_g2_valueChanged(int value);
     void on_g3_valueChanged(int value);
     void on_g4_valueChanged(int value);
     void on_g5_valueChanged(int value);
     void on_g6_valueChanged(int value);
     void on_g7_valueChanged(int value);
     void on_g8_valueChanged(int value);
     void on_g9_valueChanged(int value);
     void on_g10_valueChanged(int value);
     void on_preset_activated(int index);

     void barrachange();        // -------------------------------------funcion para pedir el valor de la barra
};




#endif // MAINWINDOW_H
