#include "image.h"

#include <QtGlobal>

Image::Image()
{
 empty = true;
 previewEmpty = true;

 previewU = 0;
 previewV = 0;

 sizeU = 0;
 sizeV = 0;

 values.resize(0);
 previewImageValues.resize(0);
}



bool Image::load(const char fileName[])
{
    //Load exr data line by line and copy it into our buffers
    try
    {

        Imf::InputFile file(fileName);
        Imath::Box2i dw = file.header().dataWindow();

        int64_t w = dw.max.x - dw.min.x + 1;
        int64_t h = dw.max.y - dw.min.y + 1;

        //Check image size
        if(w<0 || h<0 || w>(1<<20) || h>(1<<20))
        {
            qWarning("OpenEXR image size too large");
            return false;
        }

        sizeU = (unsigned int) w;
        sizeV = (unsigned int) h;


        //Check channels
        const Imf::ChannelList & chl = file.header().channels();
        if(!chl.findChannel("R") || !chl.findChannel("G") || !chl.findChannel("B"))
        {
            qWarning("Channels RGB not found in OpenEXR image");
            return false;
        }



        //Allocate values buffer
        values.resize((uint64_t)sizeU*(uint64_t)sizeV*3);


        //Allocate a temporary line buffer for reading one line of data
        QVector<float> line (sizeU*3);

        QProgressDialog progress("Ouverture de l'EXR...", "Annuler", dw.min.y+0, dw.max.y+1 );
        progress.setWindowModality(Qt::WindowModal);

        //Load image lines inverting y axis
        for(int y=dw.min.y ; y<=dw.max.y ; y++)
        {

            progress.setValue(y);
            if (progress.wasCanceled())
                return false;

            //int imgY = dw.max.y - y;

            Imf::FrameBuffer frameBuffer;
            frameBuffer.insert (	"R", // name
                        Imf::Slice (Imf::FLOAT, // type
                                (char *) (&line[0] - y*w*3 - dw.min.x*3 ),
                            sizeof (float) * 3,		// xStride
                            sizeof (float) * w*3));	// yStride

            frameBuffer.insert (	"G", // name
                        Imf::Slice (Imf::FLOAT, // type
                                (char *) (&line[1] - y*w*3 - dw.min.x*3 ),
                            sizeof (float) * 3,		// xStride
                            sizeof (float) * w*3));	// yStride

            frameBuffer.insert (	"B", // name
                        Imf::Slice (Imf::FLOAT, // type
                                (char *) (&line[2] - y*w*3 - dw.min.x*3 ),
                            sizeof (float) * 3,		// xStride
                            sizeof (float) * w*3));	// yStride

            file.setFrameBuffer (frameBuffer);
            file.readPixels (y);

            for(unsigned int x=0;x<sizeU;x++)
            {
                values[(y*sizeU+x)*3]=line[x*3];
                values[(y*sizeU+x)*3+1]=line[x*3+1];
                values[(y*sizeU+x)*3+2]=line[x*3+2];
            }
        }

         progress.setValue(dw.max.y+1);
    }
    catch(const Iex::BaseExc & exc)
    {
        qCritical("OpenEXR - %s", exc.what());
        return false;
    }
    empty = false;
    return true;
}

bool Image::writeXYZ (const char fileName[])
{
    try
    {
        float *xPixels = values.data() ;
        float *yPixels = values.data() + 1 ;
        float *zPixels = values.data() + 2 ;

        unsigned int width  = sizeU;
        unsigned int height = sizeV;

        Imf::Header header (width, height);
        header.channels().insert ("X", Imf::Channel (Imf::FLOAT));
        header.channels().insert ("Y", Imf::Channel (Imf::FLOAT));
        header.channels().insert ("Z", Imf::Channel (Imf::FLOAT));

        Imf::OutputFile file(fileName, header);

        Imf::FrameBuffer frameBuffer;
        frameBuffer.insert ("X", // name
                    Imf::Slice (Imf::FLOAT, // type
                        (char *) xPixels, // basse
                        sizeof (float) * 3, // xStride
                        sizeof (float) * width * 3)); // yStride
        frameBuffer.insert ("Y", // name
                    Imf::Slice (Imf::FLOAT, // type
                        (char *) yPixels, // basse
                        sizeof (float) * 3, // xStride
                        sizeof (float) * width * 3)); // yStride
        frameBuffer.insert ("Z", // name
                    Imf::Slice (Imf::FLOAT, // type
                        (char *) zPixels, // basse
                        sizeof (float) * 3, // xStride
                        sizeof (float) * width * 3)); // yStride

        file.setFrameBuffer (frameBuffer);
        file.writePixels (height);
    }
    catch(const Iex::BaseExc & exc)
    {
        qCritical("OpenEXR - %s", exc.what());
        return false;
    }

    return true;
}

QVector<float> Image::getPixelValue( int width , int height)
{

   QVector<float> pixel(3); //Allocate vector for R,G,B values
   uint64_t pixelstart = (width+sizeU*height)*3; //Calculate offset of the pixel stored as lines
   pixel[0]= values[pixelstart];
   pixel[1]= values[pixelstart+1];
   pixel[2]= values[pixelstart+2];

   return pixel ;
}

bool Image::isEmpty()
{
    return empty;
}

QImage Image::preview(int w ,int h, int gain)
{
    QImage previewImage(sizeU,sizeV,QImage::Format_RGB32);
    QRgb rgb;
    int R ,G ,B;
    float unsurgama = 1.0/2.2;
    float powGain = pow(2,gain*unsurgama/100);




 /*  double hMoinsUn = sizeV-1;
   double pasPhi   = (2.0*M_PI) / (double) sizeU;
   double pasTeta  = M_PI/hMoinsUn;*/

    previewImageValues.resize(sizeV*sizeU*3);

    for(unsigned int y =0; y<sizeV; y++)
    {
        for(unsigned int x =0; x<sizeU; x++)
        {

            int valueOffset = (int) (x+sizeU*y)*3;

            previewImageValues[valueOffset]     = sqrt(values[valueOffset])   * 255;
            previewImageValues[valueOffset + 1] = sqrt(values[valueOffset+1]) * 255;
            previewImageValues[valueOffset + 2] = sqrt(values[valueOffset+2]) * 255;

            R = qBound(0,(int)(previewImageValues[valueOffset]   * powGain), 255);
            G = qBound(0,(int)(previewImageValues[valueOffset+1] * powGain), 255);
            B = qBound(0,(int)(previewImageValues[valueOffset+2] * powGain), 255);

         /*  //Test pour afficher le masque
          * float tetaVec = (float) y*pasTeta;
            float phiVec  = (float) x*pasPhi;
            float tetaMask = M_PI/4.0;
            float phiMask  = M_PI;

            float cosAlpha = sin(tetaVec)*sin(tetaMask)*cos(phiVec-phiMask)+cos(tetaVec)*cos(tetaMask);
            //qDebug("offset %d value %f R=%d", valueOffset, values[valueOffset], R);
            if(acos(cosAlpha)<M_PI/18.0)
            {
            rgb= qRgb(0, 0, 0);
            previewImage.setPixel(x, y, rgb);
            }else
            {*/
            rgb= qRgb(R, G, B);
            previewImage.setPixel(x, y, rgb);

        }
    }
    previewImage = previewImage.scaled(w,h,Qt::IgnoreAspectRatio,Qt::FastTransformation);
    previewEmpty = false;
    return previewImage;
}

unsigned int Image::width()
{
    return sizeU;
}

unsigned int Image::height()
{
    return sizeV;
}

float* Image::getPixelAddress(unsigned int width ,unsigned int height)
{
   return (width < sizeU && height < sizeV) ? values.data() +(width+height*sizeU)*3:NULL;
}

void Image::setPixelValue(unsigned int x, unsigned int y,float R, float G, float B)
{

    float* pAdr = getPixelAddress(x,y);

    pAdr[0] = R;
    pAdr[1] = G;
    pAdr[2] = B;

    return;
}

void Image::addPixelValue(unsigned int x, unsigned int y, float R, float G, float B)
{
    float* pAdr = getPixelAddress(x,y);

    pAdr[0] += R;
    pAdr[1] += G;
    pAdr[2] += B;

    return;
}


