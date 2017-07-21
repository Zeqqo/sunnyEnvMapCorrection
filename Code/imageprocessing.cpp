#include "imageprocessing.h"

bool  ImageProcessing::isVecInMask(float tetaVec , float phiVec , float tetaMask , float phiMask ,float alphaMask )
{
    float cosAlpha = sin(tetaVec)*sin(tetaMask)*cos(phiVec-phiMask)+cos(tetaVec)*cos(tetaMask);
    return acos(cosAlpha)<alphaMask;
}

float ImageProcessing::isVecInMask(float tetaVec , float phiVec , float tetaMask , float phiMask ,float alphaMask1 , float alphaMask2  )
{
    float cosAlpha = sin(tetaVec)*sin(tetaMask)*cos(phiVec-phiMask)+cos(tetaVec)*cos(tetaMask);
    float alpha = acos(cosAlpha);
    if (alpha < alphaMask1)
        return 1.0;
    else if(alpha > alphaMask2)
        return 0.0;
    else
        return (alpha-alphaMask2)/(alphaMask1-alphaMask2);
}

float ImageProcessing::RGB2X(float * pixelAddress)
{
    return 0.4124564*(pixelAddress[0])+0.3575761*(pixelAddress[1])+0.1804375*(pixelAddress[2]);
}

float ImageProcessing::RGB2Y(float * pixelAddress)
{
    return 0.21267296*(pixelAddress[0])+0.7151522*(pixelAddress[1])+0.0721750*(pixelAddress[2]);
}

float ImageProcessing::RGB2Z(float * pixelAddress)
{
    return 0.0193339*(pixelAddress[0])+0.1191920*(pixelAddress[1])+0.9503041*(pixelAddress[2]);
}

float ImageProcessing::eclairementHorizontal( Image envMap ,int XYZ )
{
    /* Compute Horizontal srbnsqeir:
     *  envMap is an exr image */

    unsigned int w = envMap.width();
    unsigned int h = envMap.height();
    double hSurDeuxMoisUn = h/2-1;
    float unSurW= 1.0/((float)w);
    float eHori = 0;
    float radiance = -1;
    float contribHaute;
    float contribBasse;
    float unsursix= 1.0/6.0;
    float unsurtrois= 1.0/3.0;
    QVector<float> z(h/2+2);

    z[0]=1; // Setting z[0]= to 1 so for the first pixel the "top contribution" is null
    for(unsigned int i=0; i<h/2; i++ )
    {
        z[i+1]=cos((double) i*M_PI/(2.0*hSurDeuxMoisUn));
    }
    z[h/2+1]=0; // Setting z[h/2+1]= to 0 so for the last pixel the "bottom contribution" is null

    for(unsigned int y=0; y<h/2; y++)
    {
        radiance     =   0 ;
        contribHaute =   (unsursix*z[y]*(z[y]+z[y+1])-unsurtrois*z[y+1]*z[y+1]);
        contribBasse =  -(unsursix*z[y+2]*(z[y+2]+z[y+1])-unsurtrois*z[y+1]*z[y+1]);

        for(unsigned int x=0; x<w; x++)
        {
            radiance     += RGB2XYZ(envMap.getPixelAddress(x,y), XYZ );
        }
        eHori += radiance*(contribBasse+contribHaute)*2.0*M_PI*unSurW;

    }
    return eHori;
}

float ImageProcessing::RGB2XYZ( float * pixelAddress , int XYZ )
{
    switch (XYZ) {
    case 1:
        return RGB2X( pixelAddress );
        break;
    case 2:
        return RGB2Y( pixelAddress );
        break;
    case 3:
        return RGB2Z( pixelAddress );
        break;
    case -1:
        return pixelAddress[0];
        break;
    case -2:
        return pixelAddress[1];
        break;
    case -3:
        return pixelAddress[2];
        break;
    default:
        return RGB2Y( pixelAddress );
        break;
    }

}

QVector<float> ImageProcessing::RGB2XYZ( float * pixelAddress )
{
    QVector<float> XYZ(3,-1);

    XYZ[0] = RGB2X(pixelAddress);
    XYZ[1] = RGB2Y(pixelAddress);
    XYZ[2] = RGB2Z(pixelAddress);

    return XYZ;
}

QVector<float> ImageProcessing::XYZ2RGB( float X , float Y , float Z )
{
    QVector<float> RGB(3,-1);
    RGB[0] =  3.2404542*X -1.5371385*Y -0.4985314*Z;
    RGB[1] = -0.9692660*X +1.8760108*Y +0.0415560*Z;
    RGB[2] =  0.0556434*X -0.2040259*Y +1.0572252*Z;
    return RGB;
}

QVector<float> ImageProcessing::MatBradfordXYZs2XYZd(float Xs , float Ys , float Zs,float Xd , float Yd, float Zd)
{
    float rhoS  =  0.8951000*Xs +0.2664000*Ys -0.1614000*Zs;
    float gamaS = -0.7502000*Xs +1.7135000*Ys +0.0367000*Zs;
    float betaS =  0.0389000*Xs -0.0685000*Ys +1.0296000*Zs;

    float rhoD  =  0.8951000*Xd +0.2664000*Yd -0.1614000*Zd;
    float gamaD = -0.7502000*Xd +1.7135000*Yd +0.0367000*Zd;
    float betaD =  0.0389000*Xd -0.0685000*Yd +1.0296000*Zd;

    QVector<float> Mtemp(9,-1);

    Mtemp[0] =  0.8951000*(rhoD/rhoS);
    Mtemp[1] =  0.2664000*(rhoD/rhoS);
    Mtemp[2] = -0.1614000*(rhoD/rhoS);

    Mtemp[3] = -0.7502000*(gamaD/gamaS);
    Mtemp[4] =  1.7135000*(gamaD/gamaS);
    Mtemp[5] =  0.0367000*(gamaD/gamaS);

    Mtemp[6] =  0.0389000*(betaD/betaS);
    Mtemp[7] = -0.0685000*(betaD/betaS);
    Mtemp[8] =  1.0296000*(betaD/betaS);

    QVector<float> M(9,-1);

    M[0] =  0.9869929*Mtemp[0] -0.1470543*Mtemp[3] +0.1599627*Mtemp[6];
    M[1] =  0.9869929*Mtemp[1] -0.1470543*Mtemp[4] +0.1599627*Mtemp[7];
    M[2] =  0.9869929*Mtemp[2] -0.1470543*Mtemp[5] +0.1599627*Mtemp[8];

    M[3] =  0.4323053*Mtemp[0] +0.5183603*Mtemp[3] +0.0492912*Mtemp[6];
    M[4] =  0.4323053*Mtemp[1] +0.5183603*Mtemp[4] +0.0492912*Mtemp[7];
    M[5] =  0.4323053*Mtemp[2] +0.5183603*Mtemp[5] +0.0492912*Mtemp[8];

    M[6] = -0.0085287*Mtemp[0] +0.0400428*Mtemp[3] +0.9684867*Mtemp[6];
    M[7] = -0.0085287*Mtemp[1] +0.0400428*Mtemp[4] +0.9684867*Mtemp[7];
    M[8] = -0.0085287*Mtemp[2] +0.0400428*Mtemp[5] +0.9684867*Mtemp[8];

    return M;
}

QVector<float> ImageProcessing::MFoisXYZ( QVector<float> M, float X , float Y , float Z)
{
    QVector<float> MXYZ(3,-1);

    MXYZ[0] = M[0]*X + M[1]*Y + M[2]*Z;
    MXYZ[1] = M[3]*X + M[4]*Y + M[5]*Z;
    MXYZ[2] = M[6]*X + M[7]*Y + M[8]*Z;

    return MXYZ;
}

QVector<float> ImageProcessing::eclairementHorizontal( Image envMap , float tetaMask , float phiMask ,float alphaMask1, float alphaMask2 , int XYZ )
{
    /* Compute Horizontal ligthning with one part of the sky dome masked by a disk
     *
     * Need: -envMap an image
     *  -(tetaMask, phiMask) indicate the center of the mask in spherical coordinates
     *  -alphaMask the angle between (tetaMask, phiMask,1) and the border of the disk
     *
     * Return: QVector of 2 values
     *  -Horizontal ligthning of the image whithout pixels inside the mask
     *  -Horizontal ligthning produced by pixels inside the mask
     *
     * NOTE: the spherical vector (pi/2,0) is the left down corner */

    unsigned int w = envMap.width();
    unsigned int h = envMap.height();

    double hSurDeuxMoisUn = h/2-1;

    float unSurW         =  1.0/(float)w;
    float radiance       = -1;
    float radianceMask   = -1;
    float contribHaute   = -1;
    float contribBasse   = -1;
    float unsursix       =  1.0/6.0;
    float unsurtrois     =  1.0/3.0;
    float tetaVec,phiVec = -1;
    float vecInMask      = -1;

    QVector<float> z(h/2+2);
    QVector<float> eHori(2);


    z[0]=1; // Setting z[0]= to 1 so for the first pixel the "top contribution" is null
    for(unsigned int i=0; i<h/2; i++ )
    {
        z[i+1]=cos((double) i*M_PI/(2.0*hSurDeuxMoisUn));
    }
    z[h/2+1]=0; // Setting z[h/2+1]= to 0 so for the last pixel the "bottom contribution" is null

    for(unsigned int y=0; y<h/2; y++)
    {
        radiance     =   0 ;
        radianceMask =   0 ;
        contribHaute =   (unsursix*z[y]*(z[y]+z[y+1])-unsurtrois*z[y+1]*z[y+1]);
        contribBasse =  -(unsursix*z[y+2]*(z[y+2]+z[y+1])-unsurtrois*z[y+1]*z[y+1]);

        for(unsigned int x=0; x<w; x++)
        {
            tetaVec = (float) y*M_PI/(2.0*hSurDeuxMoisUn);
            phiVec  = (float) x*M_PI*2.0*unSurW;

            vecInMask     = isVecInMask(tetaVec ,phiVec ,tetaMask ,phiMask ,alphaMask1 ,alphaMask2);

            radianceMask += vecInMask*RGB2XYZ(envMap.getPixelAddress(x,y), XYZ);
            radiance     += (1-vecInMask)*RGB2XYZ(envMap.getPixelAddress(x,y) , XYZ);


        }
        eHori[0] += radiance*(contribBasse+contribHaute)*2.0*M_PI*unSurW ;
        eHori[1] += radianceMask*(contribBasse+contribHaute)*2.0*M_PI*unSurW ;

    }
    return eHori;
}

QVector<unsigned int> ImageProcessing::barycentrePow( Image envMap, int puissance )
{
    unsigned int w = envMap.width();
    unsigned int h = envMap.height();

    double hMoinsUn = h-1;
    double pasPhi   = (2.0*M_PI) / (double) w;
    double pasTeta  = M_PI/hMoinsUn;

    double radiance = 0;
    double sommeRadiance = 0;
    double sommeX  = 0;
    double sommeY  = 0;
    double sommeZ  = 0;

    QVector <unsigned int> barycentre(2,0);

    for(unsigned int y=0; y<h/2; y++)
    {
        for(unsigned int x=0; x<w; x++)
        {
            radiance       = pow(RGB2Y(envMap.getPixelAddress(x,y)), puissance);

            sommeRadiance += radiance;

            sommeX += radiance*sin(pasTeta*y)*cos(pasPhi*x)*sin(pasTeta*y);
            sommeY += radiance*sin(pasTeta*y)*sin(pasPhi*x)*sin(pasTeta*y);
            sommeZ += radiance*cos(pasTeta*y)*sin(pasTeta*y);

        }
    }

    sommeX /= sommeRadiance;
    sommeY /= sommeRadiance;
    sommeZ /= sommeRadiance;

    double rho  = sqrt(sommeX*sommeX+sommeY*sommeY+sommeZ*sommeZ);
    double teta = acos(sommeZ/rho);
    double phi  = atan2(sommeY,sommeX);

    if (phi < 0)
        phi += 2.0*M_PI;

    barycentre[0] = round(phi/pasPhi);
    barycentre[1] = round(teta/pasTeta);

    return barycentre;
}

QVector<unsigned int> ImageProcessing::barycentrePow(Image envMap, int puissance, float tetaMask , float phiMask ,float alphaMask )
{
    unsigned int w = envMap.width();
    unsigned int h = envMap.height();

    double hMoinsUn = h-1;
    double pasPhi   = (2.0*M_PI) / (double) w;
    double pasTeta  = M_PI/hMoinsUn;

    double radiance = 0;
    double sommeRadiance = 0;
    double sommeX  = 0;
    double sommeY  = 0;
    double sommeZ  = 0;

    QVector <unsigned int> barycentre(2,0);

    for(unsigned int y=0; y<h/2; y++)
    {
        for(unsigned int x=0; x<w; x++)
        {
            if(isVecInMask((float) y*pasTeta ,(float) x*pasPhi ,tetaMask ,phiMask ,alphaMask))
            {
                radiance       = pow(RGB2Y(envMap.getPixelAddress(x,y)), puissance);

                sommeRadiance += radiance;

                sommeX += radiance*sin(pasTeta*y)*cos(pasPhi*x)*sin(pasTeta*y);
                sommeY += radiance*sin(pasTeta*y)*sin(pasPhi*x)*sin(pasTeta*y);
                sommeZ += radiance*cos(pasTeta*y)*sin(pasTeta*y);
            }

        }
    }

    sommeX /= sommeRadiance;
    sommeY /= sommeRadiance;
    sommeZ /= sommeRadiance;

    double rho  = sqrt(sommeX*sommeX+sommeY*sommeY+sommeZ*sommeZ);
    double teta = acos(sommeZ/rho);
    double phi  = atan2(sommeY,sommeX);

    if (phi < 0)
        phi += 2.0*M_PI;

    barycentre[0] = round(phi/pasPhi);
    barycentre[1] = round(teta/pasTeta);

    return barycentre;
}

Image ImageProcessing::correctionDisqueSolaire(Image envMap, int xSoleil, int ySoleil, float tetaMask, float phiMask, float alphaMask, QVector<float> XYZTotal, QVector<float> XYZMasquee )
{
    unsigned int w = envMap.width();
    unsigned int h = envMap.height();

    QProgressDialog progress("Calcul de l'image corigee", "Annuler", 0, h+ ceil(h/2));
    progress.setWindowModality(Qt::WindowModal);
    progress.setValue(0);

    Image imgCorigee = envMap;

    float XHOutMask  = -1;
    float XHInMask   = -1;
    float YHOutMask  = -1;
    float YHInMask   = -1;
    float ZHOutMask  = -1;
    float ZHInMask   = -1;


    float rayAngSoleil = 0.53*M_PI/180.0;
    float deltaAlpha = 2.0*2.0*M_PI/(float) envMap.width();

    float diffX = -1;
    float diffY = -1;
    float diffZ = -1;

    double pixValX = -1;
    double pixValY = -1;
    double pixValZ = -1;


    QVector<float>  XH(2,-1);
    QVector<float>  YH(2,-1);
    QVector<float>  ZH(2,-1);

    QVector<float>  pixXYZ(3,-1);

    QVector<float>  MBradford(9,-1);



    XH = eclairementHorizontal( envMap, tetaMask, phiMask, alphaMask-deltaAlpha, alphaMask+deltaAlpha,1);
    XHOutMask = XH[0];
    XHInMask  = XH[1];

    YH = eclairementHorizontal( envMap, tetaMask, phiMask, alphaMask-deltaAlpha, alphaMask+deltaAlpha,2);
    YHOutMask = YH[0];
    YHInMask  = YH[1];

    ZH = eclairementHorizontal( envMap, tetaMask, phiMask, alphaMask-deltaAlpha, alphaMask+deltaAlpha,3);
    ZHOutMask = ZH[0];
    ZHInMask  = ZH[1];

    MBradford=MatBradfordXYZs2XYZd( XHOutMask ,YHOutMask ,ZHOutMask ,XYZMasquee[0] ,XYZMasquee[1] ,XYZMasquee[2]);

    QVector<float> XYZSoleilAdapte = MFoisXYZ( MBradford , XHInMask , YHInMask , ZHInMask );


    //début integrale dS*cos(teta) sur la calotte du disque solaire

    double hSurDeuxMoisUn = h/2-1;

    double unSurW         =  1.0/(float)w;
    double contribHaute   = -1;
    double contribBasse   = -1;
    double unsursix       =  1.0/6.0;
    double unsurtrois     =  1.0/3.0;
    double tetaVec,phiVec = -1;
    double radiance       = -1;
    double eclCalotteSphe = 0;

    double tetaSoleil =  (double) ySoleil*M_PI/(2.0*hSurDeuxMoisUn);
    double phiSoleil  =  (double) xSoleil*M_PI*2.0*unSurW;;

    QVector<double> z(h/2+2);

    z[0]=1; // Setting z[0]= to 1 so for the first pixel the "top contribution" is null
    for(unsigned int i=0; i<h/2; i++ )
    {
        z[i+1]=cos((double) i*M_PI/(2.0*hSurDeuxMoisUn));
    }
    z[h/2+1]=0; // Setting z[h/2+1]= to 0 so for the last pixel the "bottom contribution" is null

    for(unsigned int y=0; y<h/2; y++)
    {
        progress.setValue(y);//set value barre de chargement
        if (progress.wasCanceled())
        {
          Image cancel;
          return cancel;
        }
        radiance     =   0 ;
        contribHaute =   (unsursix*z[y]*(z[y]+z[y+1])-unsurtrois*z[y+1]*z[y+1]);
        contribBasse =  -(unsursix*z[y+2]*(z[y+2]+z[y+1])-unsurtrois*z[y+1]*z[y+1]);

        for(unsigned int x=0; x<w; x++)
        {


            tetaVec = (double) y*M_PI/(2.0*hSurDeuxMoisUn);
            phiVec  = (double) x*M_PI*2.0*unSurW;

            if(isVecInMask(tetaVec ,phiVec ,tetaSoleil ,phiSoleil ,rayAngSoleil))
                radiance     += 1;

        }
        eclCalotteSphe += radiance*(contribBasse+contribHaute)*2.0*M_PI*unSurW ;
    }

    //fin integrale

    //difference  Eclairement

    diffX= (XYZTotal[0]-XYZMasquee[0]) - XYZSoleilAdapte[0];
    diffY= (XYZTotal[1]-XYZMasquee[1]) - XYZSoleilAdapte[1];
    diffZ= (XYZTotal[2]-XYZMasquee[2]) - XYZSoleilAdapte[2];

    //valeur pixels X,Y,Z a rajouter sur le disque solaire

    pixValX = (double) diffX / eclCalotteSphe;
    pixValY = (double) diffY / eclCalotteSphe;
    pixValZ = (double) diffZ / eclCalotteSphe;

    //Correction de l'image

    int progressValue = ceil(h/2);

    for(unsigned int y=0; y<h; y++)
    {
        progress.setValue(y+progressValue);//set value barre de chargement
        if (progress.wasCanceled())
        {
            Image cancel;
            return cancel;
        }
        for(unsigned int x=0; x<w; x++)
        {


            pixXYZ = RGB2XYZ( envMap.getPixelAddress(x,y) );
            pixXYZ = MFoisXYZ( MBradford , pixXYZ);

            imgCorigee.setPixelValue( x ,y ,pixXYZ[0] ,pixXYZ[1] ,pixXYZ[2]);

            tetaVec = (double) y*M_PI/(2.0*hSurDeuxMoisUn);
            phiVec  = (double) x*M_PI*2.0*unSurW;

            if(isVecInMask(tetaVec ,phiVec ,tetaSoleil ,phiSoleil ,rayAngSoleil))
            {
               imgCorigee.addPixelValue(x ,y ,(float) pixValX ,(float) pixValY ,(float) pixValZ);
            }
        }
    }



return imgCorigee;
}


