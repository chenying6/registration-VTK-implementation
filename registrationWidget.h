#ifndef MAININTERFACE_H
#define MAININTERFACE_H
/*
*Author: Chen Ying
*Version:1.3
*Date: 2019-4-17
*/

#include "ui_registration.h"
#include <qwidget.h>
#include <QVTKWidget.h>
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkActor.h"
#include "vtkConeSource.h"
#include "vtkAxesActor.h"
#include "vtkRenderWindow.h"
#include "vtkMarchingCubes.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "vtkImageData.h"
#include "vtkDataSetMapper.h"

using namespace std;
namespace Ui {
	class marchingCubeAPP : public Ui_Form {};
}
class registrationWidget : public QWidget {
	
	Q_OBJECT
public:
	registrationWidget(QWidget *parent = 0);
	~registrationWidget();
public slots:
	void setPresentStates();
	void mapCT2Toumo();
	void mapCT2Marker();
private:
	vtkMatrix4x4 * setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	vtkMatrix4x4 * setTransformation_left(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void readCase(std::string fileName);
	void writeOBJCase();
	void writeSTLCase(vtkMarchingCubes *data);
	void getXYZRotationAngles(vtkMatrix4x4 *m);
	void outputMatrix(vtkMatrix4x4 *m);
	void constructCompositeModel();
	vtkMatrix4x4*  objTrans(vtkMatrix4x4 *m);
	vtkMatrix4x4* setCurrentMatrix(char matrix[4][4]);
	vtkMatrix4x4* getmarker2CToriginMatrix();

private:
	Ui_Form *ui;
	double *m_CTcenter;
	vtkTransform *m_CTorigin2center = vtkTransform::New();
	vtkMatrix4x4 *m_CTcenter2origin = vtkMatrix4x4::New();
	vtkPolyDataMapper *m_CTMapper = vtkPolyDataMapper::New();
	vtkActor *m_CTActor = vtkActor::New();
	vtkPolyDataMapper *m_ToumoMapper = vtkPolyDataMapper::New();
	vtkActor *m_ToumoOriginActor = vtkActor::New();
	vtkAxesActor *m_MarkerActor = vtkAxesActor::New();
	vtkAxesActor *m_worldActor = vtkAxesActor::New();
	vtkMatrix4x4 *m_Markermatrix, *m_CToriginmatrix, *m_Toumomatrix;
	vtkRenderWindow *m_renderWindow = vtkRenderWindow::New();
	vtkRenderWindow *m_exportWindow = vtkRenderWindow::New();
	vtkRenderWindowInteractor *m_renderWindowInteractor = vtkRenderWindowInteractor::New();
	vtkPolyData *compositepolydata = vtkPolyData::New();
	int firstClick2Marker=0;
	int polyOrImage = 0;
	double outMatrix[20][4];
	int index = 0;

	vtkImageData *niftiImage = vtkImageData::New();
	vtkDataSetMapper *imageMapper = vtkDataSetMapper::New();

};
#endif