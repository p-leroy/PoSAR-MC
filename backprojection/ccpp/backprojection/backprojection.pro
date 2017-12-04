TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    backprojection.c

HEADERS += \
    backprojection.h

LIBS += -lfftw3

QMAKE_CFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
