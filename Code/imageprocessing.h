#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

#include "image.h"
#include <math.h>

namespace ImageProcessing
{

bool  isVecInMask( float tetaVec , float phiVec , float tetaMask , float phiMask ,float alphamax );

float isVecInMask( float tetaVec , float phiVec , float tetaMask , float phiMask ,float alphamax1 ,float alphamax2 );
float RGB2X( float * pixelAddress );
float RGB2Y( float * pixelAddress );
float RGB2Z( float * pixelAddress );
float RGB2XYZ( float * pixelAddress , int channel );
float eclairementHorizontal( Image envMap, int XYZ );
inline float eclairementHorizontal( Image envMap )
               { return eclairementHorizontal( envMap, 2 ); }

QVector<float> RGB2XYZ( float * pixelAddress );
QVector<float> XYZ2RGB(float X , float Y , float Z );
inline QVector<float> XYZ2RGB( QVector<float> XYZ )
               { return XYZ2RGB( XYZ[0] ,XYZ[1] ,XYZ[2]);}
QVector<float> MatBradfordXYZs2XYZd(float Xs , float Ys , float Zs , float Xd , float Yd , float Zd);
inline QVector<float> MatBradfordXYZs2XYZd( QVector<float> XYZs ,QVector<float> XYZd)
               {return MatBradfordXYZs2XYZd(XYZs[0] ,XYZs[1] ,XYZs[2] ,XYZd[0] ,XYZd[1] ,XYZd[2] );}
QVector<float> MFoisXYZ( QVector<float> M, float X , float Y , float Z);
inline QVector<float> MFoisXYZ( QVector<float> M, QVector<float> XYZ)
               { return MFoisXYZ( M , XYZ[0], XYZ[1], XYZ[2] );}
QVector<float> eclairementHorizontal(Image envMap , float tetaMask , float phiMask , float alphaMask1 , float alphaMask2, int XYZ );
inline QVector<float> eclairementHorizontal(Image envMap , float tetaMask , float phiMask , float alphaMask1 , float alphaMask2)
               { return eclairementHorizontal(envMap , tetaMask , phiMask , alphaMask1 , alphaMask2 , 2);}
inline QVector<float> eclairementHorizontal(Image envMap , float tetaMask , float phiMask , float alphaMask )
               { return eclairementHorizontal(envMap , tetaMask , phiMask , alphaMask , alphaMask ,2);}
inline QVector<float> eclairementHorizontal(Image envMap , float tetaMask , float phiMask , float alphaMask ,int XYZ)
               { return eclairementHorizontal(envMap , tetaMask , phiMask , alphaMask , alphaMask , XYZ);}



QVector<unsigned int> barycentrePow(Image envMap, int puissance );
QVector<unsigned int> barycentrePow(Image envMap, int puissance , float tetaMask , float phiMask ,float alphaMask);

Image correctionDisqueSolaire(Image envMap, int xSoleil, int ySoleil, float tetaMask, float phiMask, float alphaMask, QVector<float> XYZTotal, QVector<float> XYZMasquee );
}

#endif // IMAGEPROCESSING_H
