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
#include "vtkRenderWindowInteractor.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "vtkImageData.h"
#include "vtkDataSetMapper.h"
#include "qstring.h"
using namespace std;
class InputOutput;
class Transformation;
namespace Ui {
	class marchingCubeAPP : public Ui_Form {};
}
class registrationWidget : public QWidget {
	Q_OBJECT
public:
	registrationWidget(QWidget *parent = 0);
	~registrationWidget();
public slots:
	void initiateWindow(vtkActor *CT, vtkActor *Toumo,vtkAxesActor *marker, vtkAxesActor *world, vtkRenderWindow *renWindow);
	//设置当前头模、marker板的位姿
	void setPresentStates();
	//将CT影像配准到头颅模型上
	void mapCT2Toumo();
	void mapCT2Marker();
	//将CT影像在Marker板坐标系中的表示表现出来
	void getCoorsInMarker(vtkMatrix4x4 *CT, vtkMatrix4x4 *matrix, vtkMatrix4x4 *trans);
	void setPresentStatesByMatrix();
	//改变Toumo的位姿
	void transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transToumo(vtkMatrix4x4* matrix);
	//改变CT的位姿
	void transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transCT(vtkMatrix4x4* matrix);
	//改变marker板的位姿
	void transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transMarker(vtkMatrix4x4* matrix);
	//从Marker2变换到Marker1的变换矩阵
	void TransMatrix1to2(vtkMatrix4x4* trans1To2);
private:	
	void getFloatFromQString(QString s, float*& array);
private:
	Ui_Form *ui;
	double *m_CTcenter;
	vtkTransform *m_CTorigin2center = vtkTransform::New();
	vtkMatrix4x4 *m_CTcenter2origin = vtkMatrix4x4::New();
	vtkActor *m_CTActor = vtkActor::New();
	vtkActor *m_ToumoOriginActor = vtkActor::New();
	vtkAxesActor *m_MarkerActor = vtkAxesActor::New();
	vtkAxesActor *m_worldActor = vtkAxesActor::New();
	vtkMatrix4x4 *m_Markermatrix, *m_CToriginmatrix, *m_Toumomatrix;
	vtkRenderWindow *m_renderWindow = vtkRenderWindow::New();
	vtkRenderWindow *m_exportWindow = vtkRenderWindow::New();
	vtkRenderWindowInteractor *m_renderWindowInteractor = vtkRenderWindowInteractor::New();
	vtkPolyData *compositepolydata = vtkPolyData::New();
	void printToUI(std::string t);

	vtkImageData *niftiImage = vtkImageData::New();
	vtkDataSetMapper *imageMapper = vtkDataSetMapper::New();
	InputOutput *inputOutput;
	Transformation *m_transformation;
};
#endif