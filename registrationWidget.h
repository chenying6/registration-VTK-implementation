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
#include "qstring.h"
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
	//���õ�ǰͷģ��marker���λ��
	void setPresentStates();
	//��CTӰ����׼��ͷ­ģ����
	void mapCT2Toumo();
	//��CTӰ����Marker������ϵ�еı�ʾ���ֳ���
	void mapCT2Marker();
	void setPresentStatesByMatrix();
private:
	//����������ϵ�ķ�ʽ����λ��
	vtkMatrix4x4 * setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	//����������ϵ�ķ�ʽ����λ�ˣ�ʹ����Unity�еı�ʾ��������һ��
	vtkMatrix4x4 * setTransformation_left(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	//�ı�Toumo��λ��
	void transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transToumo(vtkMatrix4x4* matrix);
	//�ı�CT��λ��
	void transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transCT(vtkMatrix4x4* matrix);
	//�ı�marker���λ��
	void transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transMarker(vtkMatrix4x4* matrix);
	//����������ļ������ͣ������ļ��Ķ�ȡ
	void readCase(std::string fileName);
	//��obj��ʽ����ģ��
	void writeOBJCase();
	//��STL��ʽ����ģ��
	void writeSTLCase(vtkMarchingCubes *data);
	//VTK�ж��ھ�ϵ��Ĭ�ϵı任˳�������ƽ�ƣ���Z,X,Y����ת˳�򣬼����Ƹ������ת�Ƕȣ�ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	void getZXYRotationAngles(vtkMatrix4x4 *m);
	//Unity�ж��ھ�ϵ��Ĭ�ϵı任˳�򣬸�����ƽ�ƣ���Y,X,Z����ת˳�򣬼����Ƹ������ת�Ƕȣ�ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	void getYXZRotationAngles(vtkMatrix4x4 *m);
	//�Զ���������ģ��
	void constructCompositeModel();
	//�����ģ�����
	void exportCompositeModel();
	//Ϊ���CTģ����׼�������ʼ��ʼ����������ÿһ�η�������ʱ��׼�����ᱨ��
	void prepareExportCTModel();
	void prepareExportMarkerModel();
	vtkMatrix4x4* setCurrentMatrix(double matrix[4][4]);
	vtkMatrix4x4* setCurrentMatrix(float* r1, float* r2, float *r3, float* r4);
	//�õ���CT�任��marker�ı任���󣬼�CT��marker����ϵ�еı�ʾ
	vtkMatrix4x4* getmarker2CToriginMatrix();
	void getFloatFromQString(QString s, float*& array);
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
	int index = 0;

	vtkImageData *niftiImage = vtkImageData::New();
	vtkDataSetMapper *imageMapper = vtkDataSetMapper::New();

};
#endif