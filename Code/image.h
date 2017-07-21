#ifndef IMAGE_H
#define IMAGE_H

#include <istream>
#include <ImathBox.h>
#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <QtGlobal>
#include <ImfChannelList.h>
#include <QVector>
#include <cstdint>
#include <QImage>
#include <QColor>
#include <math.h>
#include <QProgressDialog>

class Image
{
public:
    Image();
    QImage preview(int w , int h, int gain);
    QImage preview(QSize size, int gain) { return preview(size.width(), size.height(), gain); }
    QVector<float> getPixelValue( int width , int height);
    bool load(const char fileName[]);
    bool writeXYZ(const char fileName[]);
    bool isEmpty();
    void setPixelValue(unsigned int x, unsigned int y, float R, float G , float B);
    void addPixelValue(unsigned int x, unsigned int y, float R, float G , float B);
    float* getPixelAddress(unsigned int width , unsigned int height);
    unsigned int width();
    unsigned int height();


private:
    unsigned int sizeU;
    unsigned int sizeV;
    int previewU;
    int previewV;
    QVector<float> values;
    QVector<float> previewImageValues;
    bool empty;
    bool previewEmpty;

};

#endif // IMAGE_H
