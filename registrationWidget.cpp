/*
*Author: Chen Ying
*Version:1.3
*Date: 2019-4-17
*History change description:
change to a universal model whose origin and center are not in the same place;
*/
#include "registrationWidget.h"
#include "InputOutput.h"
#include "Transformation.h"
#include "QApplication.h"
#include <vtkProperty.h>
#include <vtkColor.h>
#include "vtksys/SystemTools.hxx"
#include <vtkRenderer.h>
#include <iostream>
#include <fstream>
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPointPicker.h"
#include "PointPickerInteractorStyle.h"
#define LINE_LEN 1000
#define PI 3.14159265
#define DELIMITER ','
std::string s_dicomName = "E:\\VTK_Project\\registration_test\\data\\DICOM\\20190408\\10470000";
std::string s_modelName = "E:\\VTK_Project\\registration_test\\build\\10470000origin.obj";
//std::string s_modelName = "E:\\AR\\Vessel\\registration\\Assets\\zidingyi.obj";
std::string s_stlName = "C:\\Users\\29477\\Desktop\\registrationTest.stl";
registrationWidget::registrationWidget(QWidget *parent)
	: QWidget(parent)
{
	ui = new Ui_Form;
	ui->setupUi(this);
	connect(ui->present, SIGNAL(clicked()), this, SLOT(setPresentStates()));
	connect(ui->transform, SIGNAL(clicked()), this, SLOT(mapCT2Toumo()));
	connect(ui->ct2axis, SIGNAL(clicked()), this, SLOT(mapCT2Marker()));
	connect(ui->mpresent, SIGNAL(clicked()), this, SLOT(setPresentStatesByMatrix()));
	inputOutput = new InputOutput();
	m_transformation = new Transformation();
	vtkMapper* mapper= inputOutput->readCase(s_modelName);
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(s_modelName);
	if(extension==" "){
		m_CTActor->SetMapper((vtkDataSetMapper *)mapper);
	    m_ToumoOriginActor->SetMapper((vtkDataSetMapper *)mapper);
	}
	else {
		m_CTActor->SetMapper((vtkPolyDataMapper *)mapper);
		m_ToumoOriginActor->SetMapper((vtkPolyDataMapper *)mapper);
	}	
	initiateWindow(m_CTActor, m_ToumoOriginActor, m_MarkerActor, m_worldActor, m_renderWindow);
	inputOutput->prepareExportModel(m_CTActor, m_exportWindow);
	vtkMatrix4x4 *trans1To2 = vtkMatrix4x4::New();
	TransMatrix1to2(trans1To2);
	double Toumo2CameraArray[4][4] = {
		0.998,  -0.021, 0.048, -157.206,
		-0.051,  0.015, 0.998, -187.029,
		-0.022, -1.012,     0, 1057.975,
		0,			 0,		0,		  1
	};
	double Toumo2MarkerArray[4][4] = {
		0.693,  -0.003, -0.723, 105.555,
		0.721,  -0.013,  0.689, -22.503,
		-0.011, -1.013, -0.020,  11.072,
		0,			 0,		 0,		  1
	};
	double Camera2MarkerArray[4][4] = {
		0.654,  -0.756, -0.022, 90.609,
		0.756,  0.654,  0.006, 211.993,
		0.010, -0.021, 1.000,  -1049.021,
		-0,		 0,	     -0,		  1
	};
	vtkMatrix4x4 *Toumo2CameraMatrix = m_transformation->setCurrentMatrix(Toumo2CameraArray);
	cout <<"头模到相机的变换矩阵"<< endl;
	//m_ToumoOriginActor->SetUserMatrix(Toumo2CameraMatrix);
	//writeOBJCase();
	ui->display->insertPlainText(QStringLiteral("头模到相机的变换矩阵:\n"));
	Toumo2CameraMatrix->Print(cout);
	std::string t= m_transformation->getYXZRotationAngles(Toumo2CameraMatrix);
	ui->display->insertPlainText("the rotation angles: ");
	printToUI(t);
	vtkMatrix4x4 *Camera2markerMatrix = m_transformation->setCurrentMatrix(Camera2MarkerArray);
	cout << "相机到Marker板的变换矩阵:" << endl;
	Camera2markerMatrix->Print(cout);
	vtkMatrix4x4* Toumo2MarkerMatrix = m_transformation->setCurrentMatrix(Toumo2MarkerArray);
	vtkMatrix4x4* Marker2ToumoMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4::Invert(Toumo2MarkerMatrix, Marker2ToumoMatrix);
	vtkMatrix4x4 *markerMatrix = vtkMatrix4x4::New();
	markerMatrix->Multiply4x4(Toumo2CameraMatrix, Marker2ToumoMatrix, markerMatrix);
	cout << "marker2到相机的变换矩阵:" << endl;
	markerMatrix->Print(cout);
	auto Marker2CameraMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
	vtkMatrix4x4::Invert(Camera2markerMatrix, Marker2CameraMatrix);
	cout << "marker2到相机的变换矩阵:" << endl;
	Marker2CameraMatrix->Print(cout);
	//m_MarkerActor->SetUserMatrix(Marker2CameraMatrix);
	//m_renderWindow->Render();
	ui->display->insertPlainText(QStringLiteral("Marker2到相机的变换的旋转角度:\n"));
	t=m_transformation->getYXZRotationAngles(Marker2CameraMatrix);
	printToUI(t);
	Marker2CameraMatrix->Multiply4x4(Marker2CameraMatrix,trans1To2,  Marker2CameraMatrix);
	cout << "Marker1到相机的变换矩阵:" << endl;
	Marker2CameraMatrix->Print(cout);
	ui->display->insertPlainText(QStringLiteral("Marker1到相机的变换的旋转角度:\n"));
	t=m_transformation->getYXZRotationAngles(Marker2CameraMatrix);
	printToUI(t);
}
registrationWidget::~registrationWidget()
{
	m_CTActor->Delete();
	m_ToumoOriginActor->Delete();
	m_MarkerActor->Delete();
	m_worldActor->Delete();
	m_renderWindow->Delete();
	m_renderWindowInteractor->Delete();
}
void registrationWidget::initiateWindow(vtkActor *CT, vtkActor *Toumo, vtkAxesActor *marker, vtkAxesActor *world, vtkRenderWindow *renWindow) {
	CT->GetProperty()->SetColor(1, 1, 0.9412);
	Toumo->GetProperty()->SetOpacity(0.9);
	Toumo->GetProperty()->SetColor(0.1, 1, 0.1);
	marker->SetTotalLength(100.6, 100.6, 100.6);
	marker->SetShaftType(0);
	marker->SetAxisLabels(0);
	marker->SetCylinderRadius(0.02);
	world->SetPosition(0, 0, 0);
	world->SetTotalLength(LINE_LEN, LINE_LEN, LINE_LEN);
	world->SetShaftType(0);
	world->SetAxisLabels(0);
	world->SetCylinderRadius(0.02);
	vtkRenderer *renderer = vtkRenderer::New();
	renWindow->SetSize(800, 800);
	renWindow->AddRenderer(renderer);
	renderer->AddActor(world);
	renderer->AddActor(CT);
	renderer->AddActor(Toumo);
	renderer->AddActor(marker);
	renderer->SetBackground(.3, .3, .5);
	ui->qvtkWidget->SetRenderWindow(renWindow);
	vtkPointPicker* pointPicker = vtkPointPicker::New();
	m_renderWindowInteractor->SetPicker(pointPicker);
	m_renderWindowInteractor->SetRenderWindow(ui->qvtkWidget->GetRenderWindow());
	PointPickerInteractorStyle* style = PointPickerInteractorStyle::New();
	m_renderWindowInteractor->SetInteractorStyle(style);
	m_renderWindowInteractor->Initialize();
}
void registrationWidget::transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = m_transformation->setTransformation_left(x, y, z, rx, ry, rz);
	transToumo(matrix);
}

void registrationWidget::transToumo(vtkMatrix4x4* matrix) {
	m_Toumomatrix = m_ToumoOriginActor->GetUserMatrix();
	matrix->Multiply4x4(matrix, m_Toumomatrix, matrix);
	m_ToumoOriginActor->SetUserMatrix(matrix);
	ui->display->insertPlainText(QStringLiteral("the present transformation of Toumo:\n"));
	std::string t=m_transformation->getYXZRotationAngles(matrix);
	printToUI(t);
}

void registrationWidget::transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	//matrix = setTransformation_left(-x, y, z, rx, -ry, -rz);
	matrix = m_transformation->setTransformation_right(x, y, z, rx, ry, rz);
	transCT(matrix);
}

void registrationWidget::transCT(vtkMatrix4x4* matrix) {
	cout << "变换后的头模的位姿：" << endl;
	matrix->Print(cout);
	m_CTActor->SetUserMatrix(matrix);
	ui->display->insertPlainText(QStringLiteral("the present transformation of CT:\n"));
	std::string t=m_transformation->getYXZRotationAngles(matrix);
	printToUI(t);
	//exportCompositeModel();
	//writeOBJCase();
}

void registrationWidget::transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz){
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = m_transformation->setTransformation_left(x, y, z, rx, ry, rz);
	transMarker(matrix);
}

void registrationWidget::transMarker(vtkMatrix4x4* matrix) {
	ui->display->insertPlainText(QStringLiteral("Marker板当前变换的旋转角度:\n"));
	std::string t=m_transformation->getYXZRotationAngles(matrix);
	printToUI(t);
	vtkTransform *transformation = vtkTransform::New();
	transformation->SetMatrix(matrix);
	m_MarkerActor->SetUserTransform(transformation);
	//writeOBJCase();
}

void registrationWidget::setPresentStates()
{
	QString qrx=ui->rx->toPlainText(); 
	float rx = atof(qrx.toStdString().c_str());
	cout << "rx的值为：" << rx<<endl;
	QString qry = ui->ry->toPlainText();
	float ry = atof(qry.toStdString().c_str());
	QString qrz = ui->rz->toPlainText();
	float rz = atof(qrz.toStdString().c_str());

	QString qx = ui->x->toPlainText();
	float x = atof(qx.toStdString().c_str());
	QString qy = ui->y->toPlainText();
	float y = atof(qy.toStdString().c_str());
	QString qz = ui->z->toPlainText();
	float z = atof(qz.toStdString().c_str());

	//transToumo(x, y, z, rx, ry, rz);
	transCT(x, y, z, rx, ry, rz);

	QString qrx_2 = ui->rx_2->toPlainText();
	float rx_2 = atof(qrx_2.toStdString().c_str());
	QString qry_2 = ui->ry_2->toPlainText();
	float ry_2 = atof(qry_2.toStdString().c_str());
	QString qrz_2 = ui->rz_2->toPlainText();
	float rz_2 = atof(qrz_2.toStdString().c_str());

	QString qx_2 = ui->x_2->toPlainText();
	float x_2 = atof(qx_2.toStdString().c_str());
	QString qy_2 = ui->y_2->toPlainText();
	float y_2 = atof(qy_2.toStdString().c_str());
	QString qz_2 = ui->z_2->toPlainText();
	float z_2 = atof(qz_2.toStdString().c_str());

	transMarker(x_2, y_2, z_2, rx_2, ry_2, rz_2);

	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::setPresentStatesByMatrix() {
	QString t2c1 = ui->t2c1->toPlainText();
	float* ft2c1 = new float[4];
	getFloatFromQString(t2c1,ft2c1);
	QString t2c2 = ui->t2c2->toPlainText();
	float* ft2c2 = new float[4];
	getFloatFromQString(t2c2, ft2c2);
	QString t2c3 = ui->t2c3->toPlainText();
	float* ft2c3 = new float[4];
	getFloatFromQString(t2c3, ft2c3);
	QString t2c4 = ui->t2c4->toPlainText();
	float* ft2c4 = new float[4];
	getFloatFromQString(t2c4, ft2c4);
	vtkMatrix4x4* matrixT2C= m_transformation->setCurrentMatrix(ft2c1, ft2c2, ft2c3, ft2c4);
	matrixT2C->Print(cout);
	//transToumo(matrixT2C);
	transCT(matrixT2C);
	QString c2m1 = ui->c2m1->toPlainText();
	cout << c2m1.toStdString() << endl;
	float* fc2m1 = new float[4];
	getFloatFromQString(c2m1, fc2m1);
	QString c2m2 = ui->c2m2->toPlainText();
	float* fc2m2 = new float[4];
	getFloatFromQString(c2m2, fc2m2);
	QString c2m3 = ui->c2m3->toPlainText();
	float* fc2m3 = new float[4];
	getFloatFromQString(c2m3, fc2m3);
	QString c2m4 = ui->t2c4->toPlainText();
	float* fc2m4 = new float[4];
	getFloatFromQString(c2m4,fc2m4);
	vtkMatrix4x4* matrixC2M = m_transformation->setCurrentMatrix(fc2m1, fc2m2, fc2m3, fc2m4);
	matrixC2M->Print(cout);
	transMarker(matrixC2M);

	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::mapCT2Toumo()
{
#pragma region test
	ui->display->insertPlainText(QStringLiteral("The marker pistion and rotation:\n"));
	std::string t=m_transformation->getYXZRotationAngles(m_MarkerActor->GetMatrix());
	printToUI(t);
	vtkMatrix4x4 *Toumo2MarkerMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4 *CT2MarkerMatrix = vtkMatrix4x4::New();
	getCoorsInMarker(m_CTActor->GetMatrix(),m_MarkerActor->GetMatrix(),CT2MarkerMatrix);
	ui->display->insertPlainText(QStringLiteral("the transformation of CT in marker coordinates:\n"));
	t=m_transformation->getYXZRotationAngles(CT2MarkerMatrix);
	printToUI(t);
	vtkMatrix4x4 *Marker2CTMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4::Invert(CT2MarkerMatrix, Marker2CTMatrix);
	//此时的Toumo2CTMatrix，是指在Marker板坐标系中，由CT变换到头模的变换矩阵
	vtkMatrix4x4 *Toumo2CTMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4::Multiply4x4(Marker2CTMatrix, Toumo2MarkerMatrix, Toumo2CTMatrix);
	ui->display->insertPlainText(QStringLiteral("the CT2Toumo matrix is:\n"));
	t=m_transformation->getYXZRotationAngles(Toumo2CTMatrix);
	printToUI(t);
	vtkMatrix4x4 *transformedCTMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4::Multiply4x4(CT2MarkerMatrix, Toumo2CTMatrix, transformedCTMatrix);
	vtkMatrix4x4::Multiply4x4(m_MarkerActor->GetMatrix(),transformedCTMatrix, transformedCTMatrix);
	vtkTransform *CTtrans = vtkTransform::New();
	CTtrans->SetMatrix(transformedCTMatrix);
	m_CTActor->SetUserTransform(CTtrans);
#pragma endregion test
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::mapCT2Marker() {
	vtkMatrix4x4* CT2MarkerMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4* backCTmatrix = vtkMatrix4x4::New();
	getCoorsInMarker(m_CTActor->GetMatrix(), m_MarkerActor->GetMatrix(), CT2MarkerMatrix);
	CT2MarkerMatrix->Invert();
	vtkMatrix4x4::Multiply4x4(CT2MarkerMatrix, m_CTActor->GetMatrix(), backCTmatrix);
	m_CTActor->SetUserMatrix(backCTmatrix);
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::getCoorsInMarker(vtkMatrix4x4 *model, vtkMatrix4x4 *marker, vtkMatrix4x4 *trans)
{
		vtkMatrix4x4* camera2MarkerMatrix = vtkMatrix4x4::New();
		vtkMatrix4x4::Invert(marker,camera2MarkerMatrix);
		ui->display->insertPlainText(QStringLiteral("the angles of model in marker coordiates:\n"));
		vtkMatrix4x4::Multiply4x4(camera2MarkerMatrix, model, trans);
		std::string t=m_transformation->getYXZRotationAngles(trans);
		printToUI(t);
}

void registrationWidget::TransMatrix1to2(vtkMatrix4x4* trans1To2) {
	trans1To2->Zero();
	trans1To2->SetElement(0, 0, -1);
	trans1To2->SetElement(1, 2, -1);
	trans1To2->SetElement(2, 1, -1);
	trans1To2->SetElement(3, 3, 1);
}
void registrationWidget::getFloatFromQString(QString s, float*& array) {
	string originS = s.toStdString().c_str();
	int index0 = originS.find(DELIMITER);
	string sub0 = originS.substr(0, index0);
	array[0] = atof(sub0.c_str());
	string sub00 = originS.substr(index0 + 1,originS.size());
	int index1 = sub00.find(DELIMITER);
	string sub1 = sub00.substr(0, index1);
	array[1] = atof(sub1.c_str());
	string sub11 = sub00.substr(index1+1, sub00.size());
	int index2 = sub11.find(DELIMITER);
	string sub2 = sub11.substr(0, index2);
	array[2] = atof(sub2.c_str());
	string sub22 = sub11.substr(index2 + 1, sub11.size());
	int index3 = sub22.find(DELIMITER);
	string sub3 =sub22.substr(0,index3);
	array[3] = atof(sub3.c_str());
}
void registrationWidget::printToUI(std::string t) {
	ui->display->insertPlainText(QString::fromStdString(t));
	ui->display->insertPlainText("\n");
}
int main(int argc, char *argv[]) {
	QApplication app(argc, argv);
	registrationWidget registrationWidget;
	registrationWidget.show();
	return app.exec();
}