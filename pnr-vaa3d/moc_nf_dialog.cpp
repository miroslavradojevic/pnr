/****************************************************************************
** Meta object code from reading C++ file 'nf_dialog.h'
**
** Created: Wed May 17 00:19:18 2017
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "nf_dialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'nf_dialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_CommonDialog[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x0a,
      23,   13,   13,   13, 0x0a,
      32,   13,   13,   13, 0x0a,
      41,   13,   13,   13, 0x0a,
      55,   13,   13,   13, 0x0a,
      66,   13,   13,   13, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_CommonDialog[] = {
    "CommonDialog\0\0accept()\0reject()\0"
    "update()\0setFilePath()\0showHelp()\0"
    "showHistory()\0"
};

const QMetaObject CommonDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_CommonDialog,
      qt_meta_data_CommonDialog, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &CommonDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *CommonDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *CommonDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_CommonDialog))
        return static_cast<void*>(const_cast< CommonDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int CommonDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: accept(); break;
        case 1: reject(); break;
        case 2: update(); break;
        case 3: setFilePath(); break;
        case 4: showHelp(); break;
        case 5: showHistory(); break;
        default: ;
        }
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
