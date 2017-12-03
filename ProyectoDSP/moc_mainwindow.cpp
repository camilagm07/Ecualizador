/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[20];
    char stringdata0[331];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 25), // "on_fileEdit_returnPressed"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 21), // "on_fileButton_clicked"
QT_MOC_LITERAL(4, 60, 28), // "on_volumeSlider_valueChanged"
QT_MOC_LITERAL(5, 89, 5), // "value"
QT_MOC_LITERAL(6, 95, 6), // "update"
QT_MOC_LITERAL(7, 102, 18), // "on_g1_valueChanged"
QT_MOC_LITERAL(8, 121, 18), // "on_g2_valueChanged"
QT_MOC_LITERAL(9, 140, 18), // "on_g3_valueChanged"
QT_MOC_LITERAL(10, 159, 18), // "on_g4_valueChanged"
QT_MOC_LITERAL(11, 178, 18), // "on_g5_valueChanged"
QT_MOC_LITERAL(12, 197, 18), // "on_g6_valueChanged"
QT_MOC_LITERAL(13, 216, 18), // "on_g7_valueChanged"
QT_MOC_LITERAL(14, 235, 18), // "on_g8_valueChanged"
QT_MOC_LITERAL(15, 254, 18), // "on_g9_valueChanged"
QT_MOC_LITERAL(16, 273, 19), // "on_g10_valueChanged"
QT_MOC_LITERAL(17, 293, 19), // "on_preset_activated"
QT_MOC_LITERAL(18, 313, 5), // "index"
QT_MOC_LITERAL(19, 319, 11) // "barrachange"

    },
    "MainWindow\0on_fileEdit_returnPressed\0"
    "\0on_fileButton_clicked\0"
    "on_volumeSlider_valueChanged\0value\0"
    "update\0on_g1_valueChanged\0on_g2_valueChanged\0"
    "on_g3_valueChanged\0on_g4_valueChanged\0"
    "on_g5_valueChanged\0on_g6_valueChanged\0"
    "on_g7_valueChanged\0on_g8_valueChanged\0"
    "on_g9_valueChanged\0on_g10_valueChanged\0"
    "on_preset_activated\0index\0barrachange"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   94,    2, 0x08 /* Private */,
       3,    0,   95,    2, 0x08 /* Private */,
       4,    1,   96,    2, 0x08 /* Private */,
       6,    0,   99,    2, 0x08 /* Private */,
       7,    1,  100,    2, 0x08 /* Private */,
       8,    1,  103,    2, 0x08 /* Private */,
       9,    1,  106,    2, 0x08 /* Private */,
      10,    1,  109,    2, 0x08 /* Private */,
      11,    1,  112,    2, 0x08 /* Private */,
      12,    1,  115,    2, 0x08 /* Private */,
      13,    1,  118,    2, 0x08 /* Private */,
      14,    1,  121,    2, 0x08 /* Private */,
      15,    1,  124,    2, 0x08 /* Private */,
      16,    1,  127,    2, 0x08 /* Private */,
      17,    1,  130,    2, 0x08 /* Private */,
      19,    0,  133,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void, QMetaType::Int,   18,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MainWindow *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->on_fileEdit_returnPressed(); break;
        case 1: _t->on_fileButton_clicked(); break;
        case 2: _t->on_volumeSlider_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->update(); break;
        case 4: _t->on_g1_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->on_g2_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->on_g3_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->on_g4_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: _t->on_g5_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->on_g6_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->on_g7_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->on_g8_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: _t->on_g9_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: _t->on_g10_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: _t->on_preset_activated((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: _t->barrachange(); break;
        default: ;
        }
    }
}

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow.data,
      qt_meta_data_MainWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 16)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 16;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
