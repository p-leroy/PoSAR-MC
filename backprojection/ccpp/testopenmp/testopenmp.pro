TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c

QMAKE_CFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
