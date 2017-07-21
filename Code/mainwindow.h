#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "image.h"
#include "imageprocessing.h"
#include <QLabel>
#include <QSlider>
#include <QLayout>
#include <QPainter>
#include <QFileDialog>


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget * parent = 0);
    ~MainWindow();

protected:
    void resizeEvent ( QResizeEvent * event );

private slots:
    void display();

    void ouvrirEXR();
    void enregistrerEXR();
    void selectMasqueMethode(int index);
    void setMasqueRayonDistRay();
    void setMasqueAngle();
    void setXYZTot();
    void setXYZMasque();
    void avantApresPreview();
    void corrigerImage();

private:
    Ui::MainWindow *ui;

    QString exrFileName;
    QString exrSaveName;

    Image exrImage;
    Image exrCorrigee;

    QImage imagePreview;

    QVector<float> XYZTot;
    QVector<float> XYZMasquee;

    QVector<unsigned int> xySoleil;
    QVector<unsigned int> xySoleilPreview;

    QPainter afficheSol;

    bool afficherSoleil;
    bool afficherMasque;
    bool afficherAvant;

    float gain;
    float angle;

    double pasPhi;
    double pasTeta;

    QPen colorPen;

};

#endif // MAINWINDOW_H
