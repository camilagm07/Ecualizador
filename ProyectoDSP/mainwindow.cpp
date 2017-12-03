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
 * \file   mainwindow.cpp
 *         Implements the equalizer H(w) computation
 * \author Pablo Alvarado
 * \date   2010.12.12
 *
 * $Id: equalizer.cpp $
 */


#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "jack.h"
#include <string>


#undef _DSP_DEBUG
#define _DSP_DEBUG

#ifdef _DSP_DEBUG
#define _debug(x) std::cerr << x
#include <iostream>
#else
#define _debug(x)
#endif


/**
 * Precision used by trimming
 */
const float MainWindow::Epsilon = 0.001;


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    verbose_(false),
    dspChanged_(true),
    iniciorep(true),
    //se inicializan las ganancias de cada slider
    gain0(25),
    gain1(25),
    gain2(25),
    gain3(25),
    gain4(25),
    gain5(25),
    gain6(25),
    gain7(25),
    gain8(25),
    gain9(25)
{
    ui->setupUi(this);
    /*
     * Set up a timer 4 times in a second to check if the user
     * changed the equalizer values, and if so, then create a new
     * filter response
     */
    timer_ = new QTimer(this);
    connect(timer_, SIGNAL(timeout()), this, SLOT(update()));
    connect(timer_, SIGNAL(timeout()), this, SLOT(barrachange()));    // senal de tiempo para actualizar el valor de las barras
    timer_->start(25); // se cambio el valor, para que actualize la barra mas rapido


    dsp_ = new dspSystem;
    jack::init(dsp_);

    // parse some command line arguments
    QStringList argv(QCoreApplication::arguments());

    QStringList::const_iterator it(argv.begin());
    while(it!=argv.end()) {
      if ((*it)=="-v" || (*it)=="--verbose") {
        verbose_=true;
      } else if ((*it).indexOf(".wav",0,Qt::CaseInsensitive)>0) {
        ui->fileEdit->setText(*it);
        std::string tmp(qPrintable(*it));
        jack::playAlso(tmp.c_str());
      }
      ++it;
    }

}


MainWindow::~MainWindow()
{
    jack::close();
    delete timer_;
    delete ui;
    delete dsp_;
}

void MainWindow::update() {
    if(dspChanged_){
        _debug("Updating" << std::endl);

        dspChanged_=false;
    }

}

/**
 * MainWindow::on_volumeSlider_valueChanged
 *Controla el volumen general del sistema
 *  value: valor actual del slider volumen principal
 */
void MainWindow::on_volumeSlider_valueChanged(int value){ // cada vez que cambia se ejecuta la rutina
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateVolume(value);  // esto es un puntero de tipo dspsystem que modifica el valor de volumen, se ubica en dspsystem.h

}


void MainWindow::on_fileButton_clicked() {
  selectedFiles_ =
      QFileDialog::getOpenFileNames(this,
                                   "Select one or more audio files to open",
                                   ui->fileEdit->text(),
                                   "WAV Files (*.wav)");

  if (!selectedFiles_.empty()) {
    ui->fileEdit->setText(*selectedFiles_.begin());

    jack::stopFiles();
    QStringList::iterator it;
    for (it=selectedFiles_.begin();it!=selectedFiles_.end();++it) {


        dsp_->updateinicio(iniciorep);


      std::string tmp(qPrintable(*it));
      jack::playAlso(tmp.c_str());
    }
  }
}


void MainWindow::on_fileEdit_returnPressed() {
  jack::stopFiles();
  std::string tmp(qPrintable(ui->fileEdit->text()));
  if (!tmp.empty()) {
    jack::playAlso(tmp.c_str());
  }
}

/**
 * MainWindow::on_g1_valueChanged
 * actualiza el valor de ganancia de la primera banda de los filtros
 * El parametro value contiene el valor actual del slider
 */


void MainWindow::on_g1_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg1(value);  // esto es un puntero de tipo dspsystem que modifica el valor de volumen, se ubica en dspsystem.h
    MainWindow::gain0=value;
    this->repaint();
}

/**
 * MainWindow::on_g2_valueChanged
 * actualiza el valor de ganancia de la segunda banda de los filtros
 * El parametro value contiene el valor actual del slider
 */

void MainWindow::on_g2_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg2(value);
    MainWindow::gain1=value;
    this->repaint();

}

/**
 * MainWindow::on_g3_valueChanged
 * actualiza el valor de ganancia de la tercera banda de los filtros
 * El parametro value contiene el valor actual del slider
 */

void MainWindow::on_g3_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg3(value);
    MainWindow::gain2=value;
    this->repaint();
}

/**
 * MainWindow::on_g4_valueChanged
 * actualiza el valor de ganancia de la cuarta banda de los filtros
 * El parametro value contiene el valor actual del slider
 */

void MainWindow::on_g4_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg4(value);
    MainWindow::gain3=value;
    this->repaint();
}

/**
 * MainWindow::on_g5_valueChanged
 * actualiza el valor de ganancia de la quinta banda de los filtros
 * El parametro value contiene el valor actual del slider
 */

void MainWindow::on_g5_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg5(value);
    MainWindow::gain4=value;
    this->repaint();
}

/**
 * MainWindow::on_g6_valueChanged
 * actualiza el valor de ganancia de la sexta banda de los filtros
 * El parametro value contiene el valor actual del slider
 */
void MainWindow::on_g6_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg6(value);
    MainWindow::gain5=value;
    this->repaint();

}

/**
 * MainWindow::on_g7_valueChanged
 * actualiza el valor de ganancia de la setima banda de los filtros
 * El parametro value contiene el valor actual del slider
 */
void MainWindow::on_g7_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg7(value);
    MainWindow::gain6=value;
    this->repaint();
}

/**
 * MainWindow::on_8_valueChanged
 * actualiza el valor de ganancia de la octava banda de los filtros
 * El parametro value contiene el valor actual del slider
 */
void MainWindow::on_g8_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg8(value);
    MainWindow::gain7=value;
    this->repaint();
}

/**
 * MainWindow::on_g9_valueChanged
 * actualiza el valor de ganancia de la novena banda de los filtros
 * El parametro value contiene el valor actual del slider
 */
void MainWindow::on_g9_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg9(value);
    MainWindow::gain8=value;
    this->repaint();

}

/**
 * MainWindow::on_g10_valueChanged
 * actualiza el valor de ganancia de la decima banda de los filtros
 * El parametro value contiene el valor actual del slider
 */
void MainWindow::on_g10_valueChanged(int value)
{
    if (!dspChanged_){
        dspChanged_=true;   //le indica al programa que se va a realizar un dsp, se define en mainwindow.h
    }
    dsp_->updateg10(value);
    MainWindow::gain9=value;
    this->repaint();
}

/**
 * MainWindow::barrachange
 * Funcion que se ejecuta en un intervalo de tiempo para actualizar las barras que muestran la energia de cada banda
 */
void MainWindow::barrachange(){
    //se actualizan los valores de la barra de acuerdo a la energia de cada banda
    ui->barra1->setValue(1*(dsp_->valorbarra1())); //se utiliza una funcion en el elemenete tipo dspSystem que retorna el valor actual de cad aenergia
    ui->barra2->setValue(1*(dsp_->valorbarra2()));
    ui->barra3->setValue(1*(dsp_->valorbarra3()));
    ui->barra4->setValue(1*(dsp_->valorbarra4()));
    ui->barra5->setValue(1*(dsp_->valorbarra5()));
    ui->barra6->setValue(1*(dsp_->valorbarra6()));
    ui->barra7->setValue(1*(dsp_->valorbarra7()));
    ui->barra8->setValue(1*(dsp_->valorbarra8()));
    ui->barra9->setValue(1*(dsp_->valorbarra9()));
    ui->barra10->setValue(1*(dsp_->valorbarra10()));
}

/**
 * MainWindow::on_preset_activated
 * Funcion que modifica el valor de los sliders definidos por cada preset
 * devuelve el valor actual del indece del elemento activado en combobox en la variable index
 */
void MainWindow::on_preset_activated(int index)  // funcion para establecer los preset
{   // se definen los presets segun el indice de la interfaz grafica
    if(index==0){//flat
        ui->g1->setValue(25);
        ui->g2->setValue(25);
        ui->g3->setValue(25);
        ui->g4->setValue(25);
        ui->g5->setValue(25);
        ui->g6->setValue(25);
        ui->g7->setValue(25);
        ui->g8->setValue(25);
        ui->g9->setValue(25);
        ui->g10->setValue(25);
    }
    else if(index==1){//clasical
        ui->g1->setValue(25);
        ui->g2->setValue(25);
        ui->g3->setValue(25);
        ui->g4->setValue(25);
        ui->g5->setValue(25);
        ui->g6->setValue(25);
        ui->g7->setValue(17);
        ui->g8->setValue(17);
        ui->g9->setValue(17);
        ui->g10->setValue(17);
    }
    else if (index==2){//club
        ui->g1->setValue(25);
        ui->g2->setValue(25);
        ui->g3->setValue(33);
        ui->g4->setValue(31);
        ui->g5->setValue(31);
        ui->g6->setValue(31);
        ui->g7->setValue(27);
        ui->g8->setValue(25);
        ui->g9->setValue(25);
        ui->g10->setValue(25);
    }
    else if(index==3){//dance
        ui->g1->setValue(31);
        ui->g2->setValue(33);
        ui->g3->setValue(27);
        ui->g4->setValue(25);
        ui->g5->setValue(25);
        ui->g6->setValue(19);
        ui->g7->setValue(17);
        ui->g8->setValue(17);
        ui->g9->setValue(25);
        ui->g10->setValue(25);
    }
    else if(index==4){//full bass&trebble
        ui->g1->setValue(33);
        ui->g2->setValue(31);
        ui->g3->setValue(25);
        ui->g4->setValue(17);
        ui->g5->setValue(21);
        ui->g6->setValue(23);
        ui->g7->setValue(33);
        ui->g8->setValue(37);
        ui->g9->setValue(39);
        ui->g10->setValue(39);
    }
    else if(index==5){//full trebble
        ui->g1->setValue(13);
        ui->g2->setValue(13);
        ui->g3->setValue(13);
        ui->g4->setValue(21);
        ui->g5->setValue(29);
        ui->g6->setValue(39);
        ui->g7->setValue(43);
        ui->g8->setValue(43);
        ui->g9->setValue(43);
        ui->g10->setValue(45);
    }
    else if(index==6){//pop
        ui->g1->setValue(27);
        ui->g2->setValue(31);
        ui->g3->setValue(33);
        ui->g4->setValue(35);
        ui->g5->setValue(33);
        ui->g6->setValue(25);
        ui->g7->setValue(21);
        ui->g8->setValue(21);
        ui->g9->setValue(23);
        ui->g10->setValue(23);
    }
    else if(index==7){//reggae
        ui->g1->setValue(25);
        ui->g2->setValue(25);
        ui->g3->setValue(25);
        ui->g4->setValue(19);
        ui->g5->setValue(25);
        ui->g6->setValue(33);
        ui->g7->setValue(33);
        ui->g8->setValue(25);
        ui->g9->setValue(25);
        ui->g10->setValue(25);
    }
    else if(index==8){//rock
        ui->g1->setValue(35);
        ui->g2->setValue(31);
        ui->g3->setValue(17);
        ui->g4->setValue(15);
        ui->g5->setValue(21);
        ui->g6->setValue(27);
        ui->g7->setValue(35);
        ui->g8->setValue(39);
        ui->g9->setValue(39);
        ui->g10->setValue(39);
    }
    else if(index==9){//techno
        ui->g1->setValue(35);
        ui->g2->setValue(32);
        ui->g3->setValue(25);
        ui->g4->setValue(18);
        ui->g5->setValue(19);
        ui->g6->setValue(25);
        ui->g7->setValue(35);
        ui->g8->setValue(36);
        ui->g9->setValue(36);
        ui->g10->setValue(35);
    }
    else if(index==10){//zero
        ui->g1->setValue(0);
        ui->g2->setValue(0);
        ui->g3->setValue(0);
        ui->g4->setValue(0);
        ui->g5->setValue(0);
        ui->g6->setValue(0);
        ui->g7->setValue(0);
        ui->g8->setValue(0);
        ui->g9->setValue(0);
        ui->g10->setValue(0);
    }
    else{  //Valor por defecto
        ui->g1->setValue(25);
        ui->g2->setValue(25);
        ui->g3->setValue(25);
        ui->g4->setValue(25);
        ui->g5->setValue(25);
        ui->g6->setValue(25);
        ui->g7->setValue(25);
        ui->g8->setValue(25);
        ui->g9->setValue(25);
        ui->g10->setValue(25);
    }
}

void MainWindow::paintEvent(QPaintEvent *e)//funcion encargada de graficar el nivel de ganancia
{
    QPoint begin,Ctrl0,Ctrl1,Ctrl2,Ctrl3,Ctrl4,Ctrl5,Ctrl6,Ctrl7,Ctrl8,Ctrl9;//QPoint es una clase que define un par ordenado
    QPoint end0,end1,end2,end3,end4,end5,end6,end7,end8,end9;
    int Posy = 150; // posicion en y del total de la grafica, en otras palabras desplaza la grafia de arriba y abajo
    int ofsetYC = 275; // desplaza las maximos y minimos de la funcion sinuseidal
    int ofsetX = 150; // posicion en x del total de la grafica, en otras palabras desplaza la grafia de derecha a izquierda
    int scale = -5;  // escala las maximos y minimos de la funcion sinuseidal

    //coordenadas en la ventana para colocar la grafica
    begin.setX(-5+ofsetX);
    begin.setY(Posy);

    Ctrl0.setX(30+ofsetX);
    Ctrl0.setY(gain0*scale + ofsetYC);
    end0.setX(65+ofsetX);
    end0.setY(Posy);

    Ctrl1.setX(100+ofsetX);
    Ctrl1.setY(gain1*scale + ofsetYC);
    end1.setX(135+ofsetX);
    end1.setY(Posy);

    Ctrl2.setX(170+ofsetX);
    Ctrl2.setY(gain2*scale + ofsetYC);
    end2.setX(205+ofsetX);
    end2.setY(Posy);

    Ctrl3.setX(240+ofsetX);
    Ctrl3.setY(gain3*scale + ofsetYC);
    end3.setX(275+ofsetX);
    end3.setY(Posy);

    Ctrl4.setX(310+ofsetX);
    Ctrl4.setY(gain4*scale + ofsetYC);
    end4.setX(345+ofsetX);
    end4.setY(Posy);

    Ctrl5.setX(380+ofsetX);
    Ctrl5.setY(gain5*scale + ofsetYC);
    end5.setX(415+ofsetX);
    end5.setY(Posy);

    Ctrl6.setX(450+ofsetX);
    Ctrl6.setY(gain6*scale + ofsetYC);
    end6.setX(485+ofsetX);
    end6.setY(Posy);

    Ctrl7.setX(520+ofsetX);
    Ctrl7.setY(gain7*scale + ofsetYC);
    end7.setX(555+ofsetX);
    end7.setY(Posy);

    Ctrl8.setX(590+ofsetX);
    Ctrl8.setY(gain8*scale + ofsetYC);
    end8.setX(625+ofsetX);
    end8.setY(Posy);

    Ctrl9.setX(660+ofsetX);
    Ctrl9.setY(gain9*scale + ofsetYC);
    end9.setX(695+ofsetX);
    end9.setY(Posy);



    // se declara la clase QPainterPath encargada de generar curvas entre otras figuras
    QPainterPath myPath,myPath1,myPath2,myPath3,myPath4,myPath5,myPath6,myPath7,myPath8,myPath9;
    myPath.moveTo(begin);//mueve el punto de inicio de la grafica
    myPath.quadTo(Ctrl0,end0); // quadTo es el metodo encargado de generar una curva sinuseidal, Ctrl mueve la cresta y end determina el punto final

    myPath1.moveTo(end0);
    myPath1.quadTo(Ctrl1,end1);

    myPath2.moveTo(end1);
    myPath2.quadTo(Ctrl2,end2);

    myPath3.moveTo(end2);
    myPath3.quadTo(Ctrl3,end3);

    myPath4.moveTo(end3);
    myPath4.quadTo(Ctrl4,end4);

    myPath5.moveTo(end4);
    myPath5.quadTo(Ctrl5,end5);

    myPath6.moveTo(end5);
    myPath6.quadTo(Ctrl6,end6);

    myPath7.moveTo(end6);
    myPath7.quadTo(Ctrl7,end7);

    myPath8.moveTo(end7);
    myPath8.quadTo(Ctrl8,end8);

    myPath9.moveTo(end8);
    myPath9.quadTo(Ctrl9,end9);

    QLinearGradient myGradient; //genera un sombreado debajo de la grafica
    myGradient.setColorAt(0,QColor(0, 204, 102,200)); //le da el color deseado en formato RGB

    QPen myPen; // tipo de lines
    myPen.setColor(QColor(0, 153, 76));//color
    myPen.setWidth(2); // grueso


    QPainter painter(this); // Clase para pintar lo descrito anteriormente
    painter.setBrush(myGradient);
    painter.setPen(myPen);

    painter.drawPath(myPath);
    painter.drawPath(myPath1);
    painter.drawPath(myPath2);
    painter.drawPath(myPath3);
    painter.drawPath(myPath4);
    painter.drawPath(myPath5);
    painter.drawPath(myPath6);
    painter.drawPath(myPath7);
    painter.drawPath(myPath8);
    painter.drawPath(myPath9);

}
