#-------------------------------------------------
#
# Project created by QtCreator 2013-09-21T21:39:19
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = curveeditor
TEMPLATE = app


SOURCES +=\
        mainwindow.cpp \
    aboutdialog.cpp \
    svzalgebra.cpp \
    spline.cpp \
    qgraphicsitem.cpp \
    node.cpp \
    graphwidget.cpp \
    graphicsscene.cpp \
    edge.cpp \
    main.cpp

HEADERS  += mainwindow.h \
    aboutdialog.h \
    ui_mainwindow.h \
    ui_aboutdialog.h \
    SvzMatrix.h \
    SvzAlgebra.h \
    spline.h \
    qgraphicsitem.h \
    node.h \
    graphwidget.h \
    graphicsscene.h \
    enums.h \
    edge.h

FORMS    += mainwindow.ui \
    aboutdialog.ui

RESOURCES += \
    Resource.qrc
