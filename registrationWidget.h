#ifndef MAININTERFACE_H
#define MAININTERFACE_H
/*******************************************************************************
*Author: Chen Ying
*Content: The main interface.It has the following functions:
		1. transform CT and Marker panel according to the input transformation matrix;
		2. get the Euler rotation angles according to the rotation routine of Unity
		3. transform the regular Marker panel coordinate to the Marker panel coordinate in Unity
		4. export the obj model
*Date: 2019-4-17
********************************************************************************/
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
	//�ԽǶȷ�ʽ�����õ�ǰͷģ��marker���λ�ˡ�ע�⣬��ʱ�����ӦΪ����ϵ�ж�Ӧ����ת�ǶȺ�ƽ����
	//�÷���������֤����õ���ŷ�����Ƿ���Unity�ж�ӦCT��ŷ����ƥ�䡣
	void on_presentButton_clicked();
	//ʹ�õ�ǰ����ı任����ֱ������marker���CTģ�͵�λ��
	void on_mpresentButton_clicked();
	//ʹ��txt�ļ��еľ���ֱ������marker���CTģ�͵�λ��
	void on_txtPresentButton_clicked();
	//��CTӰ����׼��ͷ­ģ����
	void on_transformButton_clicked();
	//��ʾCTģ����Marker������ϵ�е�λ��
	void on_ct2axisButton_clicked();
	//��obj��ʽ����ģ��
	void on_writeOBJButton_clicked();
private:
	void initiateWindow(vtkActor *CT, vtkActor *Toumo,vtkAxesActor *marker, vtkAxesActor *world, vtkRenderWindow *renWindow);
	//��CTӰ����Marker������ϵ�еı�ʾ��ʾ����
	void getCoorsInMarker(vtkMatrix4x4 *matrix, vtkMatrix4x4 *CT, vtkMatrix4x4 *trans);
	//���A�任��B�þ���Ҳ��B��A�ı任��B��A����ϵ�µı�ʾ�����У�AB��ʾͬһ������ϵ�£�����AB�ı�ʾ
	void getMatrixTransAtoB(vtkMatrix4x4 *A, vtkMatrix4x4 *B, vtkMatrix4x4 *trans);
	void transModel(vtkActor* model, std::vector<int> array, std::string modelName);
	void transModel(vtkActor* model, vtkMatrix4x4* matrix, std::string modelName);
	//��Marker2�任��Marker1�ı任����
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