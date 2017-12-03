#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
 //   signal(SIGPIPE, signal_callback_handler);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
