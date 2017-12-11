TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    backprojection.c

HEADERS += \
    backprojection.h

LIBS += -lfftw3

#QMAKE_CFLAGS += -fopenmp -fprofile-arcs -ftest-coverage
#QMAKE_LFLAGS += -fopenmp -fprofile-arcs -ftest-coverage

QMAKE_CFLAGS += -fopenmp --coverage
QMAKE_LFLAGS += -fopenmp --coverage
