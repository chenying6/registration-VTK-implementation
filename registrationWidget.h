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
	//���õ�ǰͷģ��marker���λ��
	void setPresentStates();
	//��CTӰ����׼��ͷ­ģ����
	void mapCT2Toumo();
	void mapCT2Marker();
	//��CTӰ����Marker������ϵ�еı�ʾ���ֳ���
	void getCoorsInMarker(vtkMatrix4x4 *matrix, vtkMatrix4x4 *CT, vtkMatrix4x4 *trans);
	//���A�任��B�þ���Ҳ��B��A�ı任��B��A����ϵ�µı�ʾ�����У�AB��ʾͬһ������ϵ�£�����AB�ı�ʾ
	void getMatrixTransAtoB(vtkMatrix4x4 *A, vtkMatrix4x4 *B, vtkMatrix4x4 *trans);
	void setPresentStatesByMatrix();
	//�ı�Toumo��λ��
	void transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transToumo(vtkMatrix4x4* matrix);
	//�ı�CT��λ��
	void transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transCT(vtkMatrix4x4* matrix);
	//�ı�marker���λ��
	void transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transMarker(vtkMatrix4x4* matrix);
	//��Marker2�任��Marker1�ı任����
	void TransMatrix1to2(vtkMatrix4x4* trans1To2);
	void getFloatFromQString(QString s, float*& array);
	void printToUI(std::string t);
private:
	Ui_Form *ui;
	vtkActor *m_CTActor = vtkActor::New();
	vtkActor *m_ToumoOriginActor = vtkActor::New();
	vtkAxesActor *m_MarkerActor = vtkAxesActor::New();
	vtkAxesActor *m_worldActor = vtkAxesActor::New();
	vtkRenderWindow *m_renderWindow = vtkRenderWindow::New();
	vtkRenderWindow *m_exportWindow = vtkRenderWindow::New();
	vtkRenderWindowInteractor *m_renderWindowInteractor = vtkRenderWindowInteractor::New();

	InputOutput *inputOutput;
	Transformation *m_transformation;
};
#endif