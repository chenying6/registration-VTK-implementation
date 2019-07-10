#ifndef MAININTERFACE_H
#define MAININTERFACE_H
/*
*Author: Chen Ying
*Content: The main interface.It has the following functions:
*Date: 2019-4-17
*/
#include <vector>
#include "ui_registration.h"
#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
class InputOutput;
class Transformation;
namespace Ui {
	class registrationWidget : public Ui_Form {};
}
class registrationWidget : public QWidget {
	Q_OBJECT

public:
	registrationWidget(QWidget *parent = 0);
	~registrationWidget();

public slots:
	//设置当前头模、marker板的位姿
	void on_presentButton_clicked();
	void on_mpresentButton_clicked();
	//将CT影像配准到头颅模型上
	void on_transformButton_clicked();
	void on_ct2axisButton_clicked();
	//以obj格式导出模型
	void on_writeOBJButton_clicked();
	void on_txtPresentButton_clicked();
private:
	void initiateWindow(vtkActor *CT, vtkActor *Toumo,vtkAxesActor *marker, vtkAxesActor *world, vtkRenderWindow *renWindow);
	//将CT影像在Marker板坐标系中的表示表现出来
	void getCoorsInMarker(vtkMatrix4x4 *matrix, vtkMatrix4x4 *CT, vtkMatrix4x4 *trans);
	//求把A变换到B得矩阵，也即B到A的变换，B在A坐标系下的表示。其中，AB表示同一个坐标系下，对象AB的表示
	void getMatrixTransAtoB(vtkMatrix4x4 *A, vtkMatrix4x4 *B, vtkMatrix4x4 *trans);
	void transModel(vtkActor* model, std::vector<int> array, std::string modelName);
	void transModel(vtkActor* model, vtkMatrix4x4* matrix, std::string modelName);
	//从Marker2变换到Marker1的变换矩阵
	void TransMatrix1to2(vtkMatrix4x4* trans1To2);
	void GetMatrix1Angles();
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

	InputOutput *m_pInputOutput;
	Transformation *m_transformation;
};
#endif MAININTERFACE_H