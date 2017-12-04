TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    ../backprojection/backprojection.c

LIBS += -lfftw3

HEADERS += \
    ../backprojection/backprojection.h

INCLUDEPATH += ../backprojection

QMAKE_CFLAGS += -O3
QMAKE_CFLAGS += -fopenmp

QMAKE_LFLAGS += -fopenmp
