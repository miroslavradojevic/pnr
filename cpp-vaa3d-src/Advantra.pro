TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external
INCLUDEPATH	+= ./
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/common_lib/include

HEADERS	+= Advantra_plugin.h \
    connected.h \
    frangi.h \
    nf_dialog.h \
    node.h \
    seed.h \
    toolbox.h \
    tracker.h
SOURCES	+= Advantra_plugin.cpp \
    frangi.cpp \
    node.cpp \
    seed.cpp \
    toolbox.cpp \
    tracker.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp

TARGET	= $$qtLibraryTarget(Advantra)
DESTDIR	= $$VAA3DPATH/bin/plugins/neuron_tracing/Advantra/
