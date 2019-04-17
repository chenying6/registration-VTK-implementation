#ifndef MAININTERFACE_H
#define MAININTERFACE_H
/*
*Author: Chen Ying
*Version:1.3
*Date: 2019-4-17
*/

#include "ui_transparency.h"
#include <qwidget.h>
#include <QVTKWidget.h>
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkConeSource.h"
#include "vtkAxesActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
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
	void setCurrentTransformation();
	void mapCT2Toumo();
	void mapCT2Marker();
private:
	void setTransformation_right(const int flag, const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void setTransformation_left(const int flag, const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void readCase(std::string fileName);
	void writeOBJCase();
	void getXYZRotationAngles(vtkMatrix4x4 *m);
	void outputMatrix(vtkMatrix4x4 *m);

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
	vtkConeSource *m_conedata = vtkConeSource::New();
	vtkPolyDataMapper *m_coneMapper = vtkPolyDataMapper::New();
	vtkActor *m_coneActor = vtkActor::New();
	vtkRenderWindow *m_renderWindow = vtkRenderWindow::New();
	vtkRenderWindow *m_exportWindow = vtkRenderWindow::New();
	vtkRenderWindowInteractor *m_renderWindowInteractor = vtkRenderWindowInteractor::New();
	int firstClick2Marker=0;
	int polyOrImage = 0;
	double outMatrix[20][4];
	int index = 0;
};
#endif