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
	//设置当前头模、marker板的位姿
	void setPresentStates();
	//将CT影像配准到头颅模型上
	void mapCT2Toumo();
	//将CT影像在Marker板坐标系中的表示表现出来
	void mapCT2Marker();
	void setPresentStatesByMatrix();
private:
	//以右手坐标系的方式设置位姿
	vtkMatrix4x4 * setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	//以左手坐标系的方式设置位姿，使其在Unity中的表示与期望的一样
	vtkMatrix4x4 * setTransformation_left(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	//改变Toumo的位姿
	void transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transToumo(vtkMatrix4x4* matrix);
	//改变CT的位姿
	void transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transCT(vtkMatrix4x4* matrix);
	//改变marker板的位姿
	void transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	void transMarker(vtkMatrix4x4* matrix);
	//根据输入的文件的类型，进行文件的读取
	void readCase(std::string fileName);
	//以obj格式导出模型
	void writeOBJCase();
	//以STL格式导出模型
	void writeSTLCase(vtkMarchingCubes *data);
	//VTK中对于静系，默认的变换顺序根据先平移，再Z,X,Y的旋转顺序，计算绕各轴的旋转角度；注意，使用哪一个函数，应与setTransformation_right对应
	void getZXYRotationAngles(vtkMatrix4x4 *m);
	//Unity中对于静系，默认的变换顺序，根据先平移，再Y,X,Z的旋转顺序，计算绕各轴的旋转角度；注意，使用哪一个函数，应与setTransformation_right对应
	void getYXZRotationAngles(vtkMatrix4x4 *m);
	//自定义绘制组合模型
	void constructCompositeModel();
	//将组合模型输出
	void exportCompositeModel();
	//为输出CT模型做准备，在最开始初始化，而不是每一次发生调用时做准备（会报错）
	void prepareExportCTModel();
	void prepareExportMarkerModel();
	vtkMatrix4x4* setCurrentMatrix(double matrix[4][4]);
	vtkMatrix4x4* setCurrentMatrix(float* r1, float* r2, float *r3, float* r4);
	//得到由CT变换到marker的变换矩阵，即CT在marker坐标系中的表示
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