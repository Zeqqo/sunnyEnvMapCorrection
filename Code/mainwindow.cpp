#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    XYZMasquee.resize(3);
    XYZTot.resize(3);
    xySoleil.resize(2);
    xySoleilPreview.resize(2);
    afficherAvant = true;

    colorPen.setColor(QColor(255,0,0));



    QObject::connect(ui->gainSlider, SIGNAL(valueChanged(int)), this, SLOT(display()));
    QObject::connect(ui->actionOuvrir, SIGNAL(triggered(bool)), this, SLOT(ouvrirEXR()));
    QObject::connect(ui->actionEnregistrer_sous, SIGNAL(triggered(bool)), this, SLOT(enregistrerEXR()));
    QObject::connect(ui->rayonAngleComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(selectMasqueMethode(int)));
    QObject::connect(ui->distBox, SIGNAL(valueChanged(double)) , this , SLOT(setMasqueRayonDistRay()));
    QObject::connect(ui->rayBox, SIGNAL(valueChanged(double)) , this , SLOT(setMasqueRayonDistRay()));
    QObject::connect(ui->angleBox, SIGNAL(valueChanged(double)) , this , SLOT(setMasqueAngle()));
    QObject::connect(ui->XTotBox, SIGNAL(valueChanged(double)) , this , SLOT(setXYZTot()));
    QObject::connect(ui->YTotBox, SIGNAL(valueChanged(double)) , this , SLOT(setXYZTot()));
    QObject::connect(ui->ZTotBox, SIGNAL(valueChanged(double)) , this , SLOT(setXYZTot()));
    QObject::connect(ui->XMasqueBox, SIGNAL(valueChanged(double)) , this , SLOT(setXYZMasque()));
    QObject::connect(ui->YMasqueBox, SIGNAL(valueChanged(double)) , this , SLOT(setXYZMasque()));
    QObject::connect(ui->ZMasqueBox, SIGNAL(valueChanged(double)) , this , SLOT(setXYZMasque()));
    QObject::connect(ui->affPosCheckBox, SIGNAL(stateChanged(int)) , this ,SLOT(display()));
    QObject::connect(ui->affMasqueBox, SIGNAL(toggled(bool)) , this , SLOT(display()));
    QObject::connect(ui->rechargerButton, SIGNAL(clicked(bool)) , this , SLOT(corrigerImage()));
    QObject::connect(ui->AvantButton, SIGNAL(toggled(bool)) , this , SLOT(avantApresPreview()));

}

MainWindow::~MainWindow()
{

}

void MainWindow::resizeEvent(QResizeEvent *event)
{
   if(!exrImage.isEmpty())
        display();
}

void MainWindow::display()
{
    QSize sizeLabel = ui->imageLabel->size();
    int g = ui->gainSlider->value();

    switch (afficherAvant) {
    case true:
        imagePreview =  exrImage.preview(sizeLabel,g);
        break;
    case false:
        imagePreview =  exrCorrigee.preview(sizeLabel,g);
        break;
    }

    float pasx = (float) sizeLabel.width()/(float) exrImage.width();
    float pasy = (float) sizeLabel.height()/(float) exrImage.height();

    xySoleilPreview[0] = xySoleil[0]*pasx;
    xySoleilPreview[1] = xySoleil[1]*pasy;


    if(ui->affMasqueBox->isChecked())
    {
        double hMoinsUn = sizeLabel.height()-1;
        double pasPhi   = (2.0*M_PI) / (double) sizeLabel.width();
        double pasTeta  = M_PI/hMoinsUn;

        float tetaMask = xySoleilPreview[1]*pasTeta;
        float phiMask  = xySoleilPreview[0]*pasPhi;

       for( int yPrev=0; yPrev<sizeLabel.height(); yPrev++)
       {
           for( int xPrev=0; xPrev<sizeLabel.width(); xPrev++)
           {
               float tetaVec = (float) yPrev*pasTeta;
               float phiVec  = (float) xPrev*pasPhi;

               float cosAlpha = sin(tetaVec)*sin(tetaMask)*cos(phiVec-phiMask)+cos(tetaVec)*cos(tetaMask);

               if(acos(cosAlpha)<angle)
               {
                   imagePreview.setPixel(xPrev, yPrev, qRgb(0, 0, 0));
               }
           }
       }
    }

    if(ui->affPosCheckBox->isChecked())
    {
        afficheSol.begin(&imagePreview);
        afficheSol.setPen(colorPen);
        afficheSol.drawRect((int) xySoleilPreview[0]-5, (int) xySoleilPreview[1]-5 , 11 , 11);
        afficheSol.end();
    }

    ui->imageLabel->setPixmap(QPixmap::fromImage(imagePreview));
}

void MainWindow::ouvrirEXR()
{
   exrFileName = QFileDialog::getOpenFileName(this, "Ouvrir une image exr", QString(),"Images (*.exr *.EXR);;Tout les formats (*.*)");

   this->setCursor( QCursor( Qt::WaitCursor ) ) ;

   if(exrImage.load(exrFileName.toStdString().c_str()))
   {
       ui->gainSlider->setEnabled(true);
       ui->rechargerButton->setEnabled(true);
       ui->affMasqueBox->setEnabled(true);
       ui->affPosCheckBox->setEnabled(true);

       xySoleil = ImageProcessing::barycentrePow(exrImage, 3 );



       unsigned int w = exrImage.width();
       unsigned int h = exrImage.height();

       double hMoinsUn = h-1;
       pasPhi   = (2.0*M_PI) / (double) w;
       pasTeta  = M_PI/hMoinsUn;

       display();
   }
   else
   {
       qWarning("Image non chargee");
   }

    this->setCursor( QCursor( Qt::ArrowCursor ) ) ;

}

void MainWindow::enregistrerEXR()
{
    exrSaveName = QFileDialog::getSaveFileName(this, "Enregistrer l'image corrigée sous:", QString(),"Images (*.exr *.EXR);;Tout les formats (*.*)");
    this->setCursor( QCursor( Qt::WaitCursor ) ) ;
    if(!exrCorrigee.writeXYZ(exrSaveName.toStdString().c_str()))
        qWarning("Erreur enregistrement!");
    this->setCursor( QCursor( Qt::ArrowCursor ) ) ;

}

void MainWindow::selectMasqueMethode( int index )
{
   if(index)
   {
       ui->angleBox->setDisabled(true);
       ui->labelAngle->setDisabled(true);

       ui->labelDist->setDisabled(false);
       ui->labelRay->setDisabled(false);
       ui->distBox->setDisabled(false);
       ui->rayBox->setDisabled(false);
   } else
   {
       ui->angleBox->setDisabled(false);
       ui->labelAngle->setDisabled(false);

       ui->labelDist->setDisabled(true);
       ui->labelRay->setDisabled(true);
       ui->distBox->setDisabled(true);
       ui->rayBox->setDisabled(true);
   }
}

void MainWindow::setMasqueRayonDistRay()
{
    double dist =  ui->distBox->value();
    double ray  =  ui->rayBox->value();

    angle = atan2(ray,dist);

    display();
}

void MainWindow::setMasqueAngle()
{
    double angleDeg =  ui->angleBox->value();

    angle = M_PI*angleDeg/180.0;

    display();
}

void MainWindow::setXYZTot()
{
    XYZTot[0] = ui->XTotBox->value();
    XYZTot[1] = ui->YTotBox->value();
    XYZTot[2] = ui->ZTotBox->value();
}

void MainWindow::setXYZMasque()
{
    XYZMasquee[0] = ui->XMasqueBox->value();
    XYZMasquee[1] = ui->YMasqueBox->value();
    XYZMasquee[2] = ui->ZMasqueBox->value();
}

void MainWindow::corrigerImage()
{

   this->setCursor( QCursor( Qt::WaitCursor ) ) ;
   //qDebug("teta phi: %f %f" , xySoleil[1]*pasTeta , xySoleil[0]*pasPhi);
   exrCorrigee = ImageProcessing::correctionDisqueSolaire( exrImage ,xySoleil[0] ,xySoleil[1] ,xySoleil[1]*pasTeta ,xySoleil[0]*pasPhi ,angle ,XYZTot ,XYZMasquee );
   this->setCursor( QCursor( Qt::ArrowCursor ) ) ;
   if (!exrCorrigee.isEmpty())
   {
       ui->ApresButton->setEnabled(true);
       ui->actionEnregistrer_sous->setEnabled(true);
   }
}

void MainWindow::avantApresPreview()
{
    afficherAvant = 1 - afficherAvant;
    ui->affMasqueBox->setChecked(afficherAvant);
    ui->affPosCheckBox->setChecked(afficherAvant);
    display();
    //qDebug("Affiche avant: %d" , afficherAvant);
}
