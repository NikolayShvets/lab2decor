TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        custommodel.cpp \
        integrator.cpp \
        linearalgebra.cpp \
        main.cpp \
        mathmodel.cpp \
        simplealgorithms.cpp

HEADERS += \
    custommodel.h \
    integrator.h \
    linearalgebra.h \
    mathmodel.h \
    simplealgorithms.h

DISTFILES += \
    umlmodel.qmodel
