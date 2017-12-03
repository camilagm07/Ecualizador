/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.9.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QSlider *volumeSlider;
    QSlider *g1;
    QSlider *g2;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QSlider *g3;
    QSlider *g4;
    QSlider *g5;
    QLabel *label_4;
    QSlider *g6;
    QSlider *g7;
    QSlider *g8;
    QSlider *g9;
    QSlider *g10;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;
    QLabel *label_9;
    QLabel *label_10;
    QLabel *label_11;
    QComboBox *preset;
    QProgressBar *barra1;
    QProgressBar *barra2;
    QProgressBar *barra4;
    QProgressBar *barra3;
    QProgressBar *barra6;
    QProgressBar *barra8;
    QProgressBar *barra10;
    QProgressBar *barra5;
    QProgressBar *barra7;
    QProgressBar *barra9;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLineEdit *fileEdit;
    QToolButton *fileButton;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QToolBar *toolBar;
    QToolBar *toolBar_2;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(981, 689);
        QPalette palette;
        QBrush brush(QColor(255, 255, 255, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::WindowText, brush);
        QBrush brush1(QColor(124, 121, 127, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Button, brush1);
        QBrush brush2(QColor(186, 182, 191, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Light, brush2);
        QBrush brush3(QColor(155, 151, 159, 255));
        brush3.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Midlight, brush3);
        QBrush brush4(QColor(62, 60, 63, 255));
        brush4.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Dark, brush4);
        QBrush brush5(QColor(82, 80, 84, 255));
        brush5.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Mid, brush5);
        palette.setBrush(QPalette::Active, QPalette::Text, brush);
        palette.setBrush(QPalette::Active, QPalette::BrightText, brush);
        palette.setBrush(QPalette::Active, QPalette::ButtonText, brush);
        QBrush brush6(QColor(5, 10, 103, 255));
        brush6.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Base, brush6);
        palette.setBrush(QPalette::Active, QPalette::Window, brush1);
        QBrush brush7(QColor(0, 0, 0, 255));
        brush7.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Shadow, brush7);
        QBrush brush8(QColor(0, 255, 0, 255));
        brush8.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Link, brush8);
        palette.setBrush(QPalette::Active, QPalette::AlternateBase, brush4);
        QBrush brush9(QColor(255, 255, 220, 255));
        brush9.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::ToolTipBase, brush9);
        palette.setBrush(QPalette::Active, QPalette::ToolTipText, brush7);
        palette.setBrush(QPalette::Inactive, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Button, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Light, brush2);
        palette.setBrush(QPalette::Inactive, QPalette::Midlight, brush3);
        palette.setBrush(QPalette::Inactive, QPalette::Dark, brush4);
        palette.setBrush(QPalette::Inactive, QPalette::Mid, brush5);
        palette.setBrush(QPalette::Inactive, QPalette::Text, brush);
        palette.setBrush(QPalette::Inactive, QPalette::BrightText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush6);
        palette.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Shadow, brush7);
        palette.setBrush(QPalette::Inactive, QPalette::Link, brush8);
        palette.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush4);
        palette.setBrush(QPalette::Inactive, QPalette::ToolTipBase, brush9);
        palette.setBrush(QPalette::Inactive, QPalette::ToolTipText, brush7);
        palette.setBrush(QPalette::Disabled, QPalette::WindowText, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::Button, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Light, brush2);
        palette.setBrush(QPalette::Disabled, QPalette::Midlight, brush3);
        palette.setBrush(QPalette::Disabled, QPalette::Dark, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::Mid, brush5);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::BrightText, brush);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Shadow, brush7);
        palette.setBrush(QPalette::Disabled, QPalette::Link, brush8);
        palette.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipBase, brush9);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipText, brush7);
        MainWindow->setPalette(palette);
        MainWindow->setCursor(QCursor(Qt::ArrowCursor));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        volumeSlider = new QSlider(centralWidget);
        volumeSlider->setObjectName(QStringLiteral("volumeSlider"));
        volumeSlider->setGeometry(QRect(80, 209, 20, 391));
        volumeSlider->setCursor(QCursor(Qt::ClosedHandCursor));
        volumeSlider->setMinimum(1);
        volumeSlider->setMaximum(50);
        volumeSlider->setValue(1);
        volumeSlider->setSliderPosition(1);
        volumeSlider->setOrientation(Qt::Vertical);
        volumeSlider->setTickPosition(QSlider::TicksBothSides);
        volumeSlider->setTickInterval(5);
        g1 = new QSlider(centralWidget);
        g1->setObjectName(QStringLiteral("g1"));
        g1->setGeometry(QRect(170, 209, 20, 191));
        g1->setCursor(QCursor(Qt::ClosedHandCursor));
        g1->setMinimum(0);
        g1->setMaximum(50);
        g1->setValue(25);
        g1->setSliderPosition(25);
        g1->setOrientation(Qt::Vertical);
        g1->setTickPosition(QSlider::TicksBothSides);
        g1->setTickInterval(5);
        g2 = new QSlider(centralWidget);
        g2->setObjectName(QStringLiteral("g2"));
        g2->setGeometry(QRect(240, 210, 20, 191));
        g2->setCursor(QCursor(Qt::ClosedHandCursor));
        g2->setMinimum(0);
        g2->setMaximum(50);
        g2->setSliderPosition(25);
        g2->setOrientation(Qt::Vertical);
        g2->setTickPosition(QSlider::TicksBothSides);
        g2->setTickInterval(5);
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(60, 600, 59, 14));
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(300, 410, 59, 14));
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(510, 410, 41, 16));
        g3 = new QSlider(centralWidget);
        g3->setObjectName(QStringLiteral("g3"));
        g3->setGeometry(QRect(310, 210, 20, 191));
        g3->setCursor(QCursor(Qt::ClosedHandCursor));
        g3->setMinimum(0);
        g3->setMaximum(50);
        g3->setValue(25);
        g3->setSliderPosition(25);
        g3->setOrientation(Qt::Vertical);
        g3->setInvertedControls(false);
        g3->setTickPosition(QSlider::TicksBothSides);
        g3->setTickInterval(5);
        g4 = new QSlider(centralWidget);
        g4->setObjectName(QStringLiteral("g4"));
        g4->setGeometry(QRect(380, 209, 20, 191));
        g4->setCursor(QCursor(Qt::ClosedHandCursor));
        g4->setMinimum(0);
        g4->setMaximum(50);
        g4->setValue(25);
        g4->setSliderPosition(25);
        g4->setOrientation(Qt::Vertical);
        g4->setTickPosition(QSlider::TicksBothSides);
        g4->setTickInterval(5);
        g5 = new QSlider(centralWidget);
        g5->setObjectName(QStringLiteral("g5"));
        g5->setGeometry(QRect(450, 209, 20, 191));
        g5->setCursor(QCursor(Qt::ClosedHandCursor));
        g5->setMinimum(0);
        g5->setMaximum(50);
        g5->setValue(25);
        g5->setSliderPosition(25);
        g5->setOrientation(Qt::Vertical);
        g5->setTickPosition(QSlider::TicksBothSides);
        g5->setTickInterval(5);
        label_4 = new QLabel(centralWidget);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(580, 410, 41, 16));
        g6 = new QSlider(centralWidget);
        g6->setObjectName(QStringLiteral("g6"));
        g6->setGeometry(QRect(520, 209, 20, 191));
        g6->setCursor(QCursor(Qt::ClosedHandCursor));
        g6->setMinimum(0);
        g6->setMaximum(50);
        g6->setValue(25);
        g6->setSliderPosition(25);
        g6->setOrientation(Qt::Vertical);
        g6->setTickPosition(QSlider::TicksBothSides);
        g6->setTickInterval(5);
        g7 = new QSlider(centralWidget);
        g7->setObjectName(QStringLiteral("g7"));
        g7->setGeometry(QRect(590, 210, 20, 191));
        g7->setCursor(QCursor(Qt::ClosedHandCursor));
        g7->setMinimum(0);
        g7->setMaximum(50);
        g7->setValue(25);
        g7->setSliderPosition(25);
        g7->setOrientation(Qt::Vertical);
        g7->setTickPosition(QSlider::TicksBothSides);
        g7->setTickInterval(5);
        g8 = new QSlider(centralWidget);
        g8->setObjectName(QStringLiteral("g8"));
        g8->setGeometry(QRect(660, 210, 20, 191));
        g8->setCursor(QCursor(Qt::ClosedHandCursor));
        g8->setMinimum(0);
        g8->setMaximum(50);
        g8->setValue(25);
        g8->setSliderPosition(25);
        g8->setOrientation(Qt::Vertical);
        g8->setTickPosition(QSlider::TicksBothSides);
        g8->setTickInterval(5);
        g9 = new QSlider(centralWidget);
        g9->setObjectName(QStringLiteral("g9"));
        g9->setGeometry(QRect(730, 209, 20, 191));
        g9->setCursor(QCursor(Qt::ClosedHandCursor));
        g9->setMinimum(0);
        g9->setMaximum(50);
        g9->setValue(25);
        g9->setSliderPosition(25);
        g9->setOrientation(Qt::Vertical);
        g9->setTickPosition(QSlider::TicksBothSides);
        g9->setTickInterval(5);
        g10 = new QSlider(centralWidget);
        g10->setObjectName(QStringLiteral("g10"));
        g10->setGeometry(QRect(800, 209, 20, 191));
        g10->setCursor(QCursor(Qt::ClosedHandCursor));
        g10->setMinimum(0);
        g10->setMaximum(50);
        g10->setValue(25);
        g10->setSliderPosition(25);
        g10->setOrientation(Qt::Vertical);
        g10->setTickPosition(QSlider::TicksBothSides);
        g10->setTickInterval(5);
        label_5 = new QLabel(centralWidget);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(230, 410, 41, 16));
        label_6 = new QLabel(centralWidget);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(160, 410, 51, 16));
        label_7 = new QLabel(centralWidget);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(370, 410, 59, 20));
        label_8 = new QLabel(centralWidget);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(440, 410, 51, 20));
        label_9 = new QLabel(centralWidget);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(650, 410, 41, 20));
        label_10 = new QLabel(centralWidget);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setGeometry(QRect(720, 410, 41, 20));
        label_11 = new QLabel(centralWidget);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setGeometry(QRect(790, 410, 51, 20));
        preset = new QComboBox(centralWidget);
        preset->setObjectName(QStringLiteral("preset"));
        preset->setGeometry(QRect(830, 210, 81, 22));
        preset->setEditable(false);
        barra1 = new QProgressBar(centralWidget);
        barra1->setObjectName(QStringLiteral("barra1"));
        barra1->setGeometry(QRect(170, 440, 16, 181));
        QPalette palette1;
        QBrush brush10(QColor(126, 126, 126, 255));
        brush10.setStyle(Qt::SolidPattern);
        palette1.setBrush(QPalette::Active, QPalette::Base, brush10);
        QBrush brush11(QColor(0, 170, 255, 255));
        brush11.setStyle(Qt::SolidPattern);
        palette1.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        QBrush brush12(QColor(108, 131, 128, 255));
        brush12.setStyle(Qt::SolidPattern);
        palette1.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette1.setBrush(QPalette::Active, QPalette::Link, brush8);
        palette1.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette1.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette1.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette1.setBrush(QPalette::Inactive, QPalette::Link, brush8);
        palette1.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        QBrush brush13(QColor(145, 141, 126, 255));
        brush13.setStyle(Qt::SolidPattern);
        palette1.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette1.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette1.setBrush(QPalette::Disabled, QPalette::Link, brush8);
        barra1->setPalette(palette1);
        barra1->setValue(24);
        barra1->setOrientation(Qt::Vertical);
        barra2 = new QProgressBar(centralWidget);
        barra2->setObjectName(QStringLiteral("barra2"));
        barra2->setGeometry(QRect(240, 440, 16, 181));
        QPalette palette2;
        palette2.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette2.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette2.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        QBrush brush14(QColor(0, 0, 255, 255));
        brush14.setStyle(Qt::SolidPattern);
        palette2.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette2.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette2.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette2.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette2.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette2.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette2.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette2.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        QBrush brush15(QColor(120, 255, 30, 255));
        brush15.setStyle(Qt::SolidPattern);
        palette2.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra2->setPalette(palette2);
        barra2->setValue(24);
        barra2->setOrientation(Qt::Vertical);
        barra4 = new QProgressBar(centralWidget);
        barra4->setObjectName(QStringLiteral("barra4"));
        barra4->setGeometry(QRect(380, 440, 16, 181));
        QPalette palette3;
        palette3.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette3.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette3.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette3.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette3.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette3.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette3.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette3.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette3.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette3.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette3.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette3.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra4->setPalette(palette3);
        barra4->setValue(24);
        barra4->setOrientation(Qt::Vertical);
        barra3 = new QProgressBar(centralWidget);
        barra3->setObjectName(QStringLiteral("barra3"));
        barra3->setGeometry(QRect(310, 440, 16, 181));
        QPalette palette4;
        palette4.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette4.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette4.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette4.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette4.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette4.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette4.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette4.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette4.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette4.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette4.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette4.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra3->setPalette(palette4);
        barra3->setValue(24);
        barra3->setOrientation(Qt::Vertical);
        barra6 = new QProgressBar(centralWidget);
        barra6->setObjectName(QStringLiteral("barra6"));
        barra6->setGeometry(QRect(520, 440, 16, 181));
        QPalette palette5;
        palette5.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette5.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette5.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette5.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette5.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette5.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette5.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette5.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette5.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette5.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette5.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette5.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra6->setPalette(palette5);
        barra6->setValue(24);
        barra6->setOrientation(Qt::Vertical);
        barra8 = new QProgressBar(centralWidget);
        barra8->setObjectName(QStringLiteral("barra8"));
        barra8->setGeometry(QRect(660, 440, 16, 181));
        QPalette palette6;
        palette6.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette6.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette6.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette6.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette6.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette6.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette6.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette6.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette6.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette6.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette6.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette6.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra8->setPalette(palette6);
        barra8->setValue(24);
        barra8->setOrientation(Qt::Vertical);
        barra10 = new QProgressBar(centralWidget);
        barra10->setObjectName(QStringLiteral("barra10"));
        barra10->setGeometry(QRect(800, 440, 16, 181));
        QPalette palette7;
        palette7.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette7.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette7.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette7.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette7.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette7.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette7.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette7.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette7.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette7.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette7.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette7.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra10->setPalette(palette7);
        barra10->setValue(24);
        barra10->setOrientation(Qt::Vertical);
        barra5 = new QProgressBar(centralWidget);
        barra5->setObjectName(QStringLiteral("barra5"));
        barra5->setGeometry(QRect(450, 440, 16, 181));
        QPalette palette8;
        palette8.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette8.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette8.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette8.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette8.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette8.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette8.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette8.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette8.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette8.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette8.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette8.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra5->setPalette(palette8);
        barra5->setValue(24);
        barra5->setOrientation(Qt::Vertical);
        barra7 = new QProgressBar(centralWidget);
        barra7->setObjectName(QStringLiteral("barra7"));
        barra7->setGeometry(QRect(590, 440, 16, 181));
        QPalette palette9;
        palette9.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette9.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette9.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette9.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette9.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette9.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette9.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette9.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette9.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette9.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette9.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette9.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra7->setPalette(palette9);
        barra7->setValue(24);
        barra7->setOrientation(Qt::Vertical);
        barra9 = new QProgressBar(centralWidget);
        barra9->setObjectName(QStringLiteral("barra9"));
        barra9->setGeometry(QRect(730, 440, 16, 181));
        QPalette palette10;
        palette10.setBrush(QPalette::Active, QPalette::Base, brush10);
        palette10.setBrush(QPalette::Active, QPalette::Highlight, brush11);
        palette10.setBrush(QPalette::Active, QPalette::HighlightedText, brush12);
        palette10.setBrush(QPalette::Active, QPalette::Link, brush14);
        palette10.setBrush(QPalette::Inactive, QPalette::Base, brush10);
        palette10.setBrush(QPalette::Inactive, QPalette::Highlight, brush11);
        palette10.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush12);
        palette10.setBrush(QPalette::Inactive, QPalette::Link, brush14);
        palette10.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette10.setBrush(QPalette::Disabled, QPalette::Highlight, brush13);
        palette10.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush12);
        palette10.setBrush(QPalette::Disabled, QPalette::Link, brush15);
        barra9->setPalette(palette10);
        barra9->setValue(24);
        barra9->setOrientation(Qt::Vertical);
        widget = new QWidget(centralWidget);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(60, 10, 899, 49));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setSpacing(3);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        fileEdit = new QLineEdit(widget);
        fileEdit->setObjectName(QStringLiteral("fileEdit"));
        QPalette palette11;
        QBrush brush16(QColor(95, 95, 95, 255));
        brush16.setStyle(Qt::SolidPattern);
        palette11.setBrush(QPalette::Active, QPalette::Base, brush16);
        palette11.setBrush(QPalette::Inactive, QPalette::Base, brush16);
        palette11.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        fileEdit->setPalette(palette11);

        horizontalLayout->addWidget(fileEdit);

        fileButton = new QToolButton(widget);
        fileButton->setObjectName(QStringLiteral("fileButton"));

        horizontalLayout->addWidget(fileButton);

        MainWindow->setCentralWidget(centralWidget);
        volumeSlider->raise();
        g1->raise();
        g2->raise();
        label->raise();
        label_2->raise();
        label_3->raise();
        g3->raise();
        g4->raise();
        g5->raise();
        label_4->raise();
        g6->raise();
        g7->raise();
        g8->raise();
        g9->raise();
        g10->raise();
        label_5->raise();
        label_6->raise();
        label_7->raise();
        label_8->raise();
        label_9->raise();
        label_10->raise();
        label_11->raise();
        preset->raise();
        barra1->raise();
        barra2->raise();
        barra4->raise();
        barra3->raise();
        barra6->raise();
        barra8->raise();
        barra10->raise();
        barra5->raise();
        barra7->raise();
        barra9->raise();
        fileEdit->raise();
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 981, 19));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);
        toolBar = new QToolBar(MainWindow);
        toolBar->setObjectName(QStringLiteral("toolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);
        toolBar_2 = new QToolBar(MainWindow);
        toolBar_2->setObjectName(QStringLiteral("toolBar_2"));
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar_2);

        retranslateUi(MainWindow);

        preset->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", Q_NULLPTR));
        label->setText(QApplication::translate("MainWindow", "Volumen", Q_NULLPTR));
        label_2->setText(QApplication::translate("MainWindow", " 125 Hz", Q_NULLPTR));
        label_3->setText(QApplication::translate("MainWindow", "1 kHz", Q_NULLPTR));
        label_4->setText(QApplication::translate("MainWindow", "2 kHz", Q_NULLPTR));
        label_5->setText(QApplication::translate("MainWindow", " 63 Hz", Q_NULLPTR));
        label_6->setText(QApplication::translate("MainWindow", " 31.5 Hz", Q_NULLPTR));
        label_7->setText(QApplication::translate("MainWindow", " 250 Hz", Q_NULLPTR));
        label_8->setText(QApplication::translate("MainWindow", " 500 Hz", Q_NULLPTR));
        label_9->setText(QApplication::translate("MainWindow", "4 kHz", Q_NULLPTR));
        label_10->setText(QApplication::translate("MainWindow", "8 kHz", Q_NULLPTR));
        label_11->setText(QApplication::translate("MainWindow", "16 kHz", Q_NULLPTR));
        preset->clear();
        preset->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "FLAT", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Clasica", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Club", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Dance", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Full Bass & Trebble", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Full Trebble", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Pop", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Reggae", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Rock", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Techno", Q_NULLPTR)
         << QApplication::translate("MainWindow", "Zero", Q_NULLPTR)
        );
        preset->setCurrentText(QApplication::translate("MainWindow", "FLAT", Q_NULLPTR));
        fileButton->setText(QApplication::translate("MainWindow", "...", Q_NULLPTR));
        toolBar->setWindowTitle(QApplication::translate("MainWindow", "toolBar", Q_NULLPTR));
        toolBar_2->setWindowTitle(QApplication::translate("MainWindow", "toolBar_2", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
