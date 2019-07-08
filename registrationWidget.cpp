/*
*Author: Chen Ying
*Version:1.3
*Date: 2019-4-17
*History change description:
change to a universal model whose origin and center are not in the same place;
*/
#include <iostream>
#include <fstream>
#include "vtkProperty.h"
#include "vtkColor.h"
#include "vtksys/SystemTools.hxx"
#include "vtkPointPicker.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkDataSetMapper.h"
#include "vtkPolyDataMapper.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "registrationWidget.h"
#include "InputOutput.h"
#include "Transformation.h"
#include "PointPickerInteractorStyle.h"
#define LINE_LEN 1000
#define PI 3.14159265
#define DELIMITER ','
std::string s_dicomName = "E:\\VTK_Project\\registration_test\\data\\DICOM\\20190408\\10470000";
std::string s_modelName = "E:\\VTK_Project\\registration_test\\build\\10470000origin.obj";
registrationWidget::registrationWidget(QWidget *parent)
	: QWidget(parent)
{
	ui = new Ui_Form;
	ui->setupUi(this);
	connect(ui->present, SIGNAL(clicked()), this, SLOT(on_presentButton_clicked()));
	connect(ui->transform, SIGNAL(clicked()), this, SLOT(on_transformButton_clicked()));
	connect(ui->ct2axis, SIGNAL(clicked()), this, SLOT(on_ct2axisButton_clicked()));
	connect(ui->mpresent, SIGNAL(clicked()), this, SLOT(on_mpresentButton_clicked()));
	m_pInputOutput = new InputOutput();
	m_transformation = new Transformation();
	vtkMapper* mapper= m_pInputOutput->readCase(s_modelName);
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(s_modelName);
	if(extension==" "|| extension==".nii"){
		m_CTActor->SetMapper((vtkDataSetMapper *)mapper);
	    m_ToumoOriginActor->SetMapper((vtkDataSetMapper *)mapper);
	}
	else {
		m_CTActor->SetMapper((vtkPolyDataMapper *)mapper);
		m_ToumoOriginActor->SetMapper((vtkPolyDataMapper *)mapper);
	}	
	initiateWindow(m_CTActor, m_ToumoOriginActor, m_MarkerActor, m_worldActor, m_renderWindow);
	m_pInputOutput->prepareExportModel(m_CTActor, m_exportWindow);
	double Toumo2CameraArray[4][4] = {
		1.001,  -0.015, -0.002, 47.835,
		0.000,  -0.094, 0.998, -61.530,
		-0.013, -0.998,  -0.093,1361.201,
		0,			 0,		0,		  1
	};
	double Marker2CameraArray[4][4] = {
		-0.75227,  -0.65882, 0.00721, -170.5557,
		0.63956,  -0.72757,  0.24821, 239.1999,
		-0.15828, 0.19134, 0.96868,  1265.3409,
		-0,		 0,	     -0,		  1
	};	
	vtkMatrix4x4 *Toumo2CameraMatrix = m_transformation->setCurrentMatrix(Toumo2CameraArray);
	m_ToumoOriginActor->SetUserMatrix(Toumo2CameraMatrix);
	//inputOutput->writeOBJCase(m_exportWindow);
	ui->display->insertPlainText(QStringLiteral("头模到相机的变换矩阵:\n"));
	cout << "头模到相机的变换矩阵:\n" << endl;
	Toumo2CameraMatrix->Print(cout);
	std::string t= m_transformation->getYXZRotationAngles(Toumo2CameraMatrix);
	printToUI(t);
	vtkMatrix4x4 *Marker2CameraMatrix = m_transformation->setCurrentMatrix(Marker2CameraArray);
	cout << "Marker板到相机的变换矩阵:" << endl;
	Marker2CameraMatrix->Print(cout);
	m_MarkerActor->SetUserMatrix(Marker2CameraMatrix);
	m_renderWindow->Render();
	ui->display->insertPlainText(QStringLiteral("Marker2到相机的变换的旋转角度:\n"));
	t=m_transformation->getYXZRotationAngles(Marker2CameraMatrix);
	printToUI(t);
	vtkMatrix4x4 *trans1To2 = vtkMatrix4x4::New();
	TransMatrix1to2(trans1To2);
	Marker2CameraMatrix->Multiply4x4(Marker2CameraMatrix,trans1To2, Marker2CameraMatrix);
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
	delete m_pInputOutput;
	delete m_transformation;
}
void registrationWidget::initiateWindow(vtkActor *CT, vtkActor *Toumo, vtkAxesActor *marker, vtkAxesActor *world, vtkRenderWindow *renWindow) {
	CT->GetProperty()->SetColor(1, 1, 0.9412);
	Toumo->GetProperty()->SetColor(0.1, 1, 0.1);
	marker->SetTotalLength(LINE_LEN, LINE_LEN, LINE_LEN);
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

void registrationWidget::transModel(vtkActor* model, std::vector<int> array, std::string modelName) {
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = m_transformation->setTransformation_left(array[0],array[1], array[2], array[3], array[4], array[5]);
	transModel(model,matrix,modelName);
}

void registrationWidget::transModel(vtkActor* actor, vtkMatrix4x4* matrix, std::string modelName) {
	std::string str= modelName + " current rotation angles:\n";
	ui->display->insertPlainText(QString(QString::fromStdString(str)));
	std::string t = m_transformation->getYXZRotationAngles(matrix);
	printToUI(t);
	actor->SetUserMatrix(matrix);
}

void registrationWidget::on_presentButton_clicked()
{
	vector<int> array;
	QString qrx=ui->rx->toPlainText(); 
	array.push_back(atof(qrx.toStdString().c_str()));
	QString qry = ui->ry->toPlainText();
	array.push_back(atof(qry.toStdString().c_str()));
	QString qrz = ui->rz->toPlainText();
	array.push_back(atof(qrz.toStdString().c_str()));

	QString qx = ui->x->toPlainText();
	array.push_back(atof(qx.toStdString().c_str()));
	QString qy = ui->y->toPlainText();
	array.push_back(atof(qy.toStdString().c_str()));
	QString qz = ui->z->toPlainText();
	array.push_back(atof(qz.toStdString().c_str()));

	transModel(m_CTActor,array,"CT");

	vector<int> markerArray;
	QString qrx_2 = ui->rx_2->toPlainText();
	markerArray.push_back(atof(qrx_2.toStdString().c_str()));
	QString qry_2 = ui->ry_2->toPlainText();
	markerArray.push_back(atof(qry_2.toStdString().c_str()));
	QString qrz_2 = ui->rz_2->toPlainText();
	markerArray.push_back(atof(qrz_2.toStdString().c_str()));

	QString qx_2 = ui->x_2->toPlainText();
	markerArray.push_back(atof(qx_2.toStdString().c_str()));
	QString qy_2 = ui->y_2->toPlainText();
	markerArray.push_back(atof(qy_2.toStdString().c_str()));
	QString qz_2 = ui->z_2->toPlainText();
	markerArray.push_back(atof(qz_2.toStdString().c_str()));

	transModel((vtkActor*)m_MarkerActor, markerArray, std::string("Marker"));

	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::on_mpresentButton_clicked() {
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
	transModel(m_CTActor, matrixT2C,"CT");
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
	transModel((vtkActor*)m_MarkerActor, matrixC2M,"Marker");
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::on_transformButton_clicked()
{
	ui->display->insertPlainText(QStringLiteral("The marker pistion and rotation:\n"));
	std::string t=m_transformation->getYXZRotationAngles(m_MarkerActor->GetMatrix());
	printToUI(t);
	vtkMatrix4x4 *Toumo2MarkerMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4 *CT2MarkerMatrix = vtkMatrix4x4::New();
	getCoorsInMarker(m_MarkerActor->GetMatrix(), m_CTActor->GetMatrix(),CT2MarkerMatrix);
	getCoorsInMarker(m_MarkerActor->GetMatrix(), m_ToumoOriginActor->GetMatrix(), Toumo2MarkerMatrix);
	ui->display->insertPlainText(QStringLiteral("the transformation of CT in marker coordinates:\n"));
	t=m_transformation->getYXZRotationAngles(CT2MarkerMatrix);
	printToUI(t);
	//此时的Toumo2CTMatrix，是指在Marker板坐标系中，由CT变换到头模的变换矩阵
	vtkMatrix4x4 *Toumo2CTMatrix = vtkMatrix4x4::New();
	getMatrixTransAtoB(CT2MarkerMatrix, Toumo2MarkerMatrix, Toumo2CTMatrix);
	ui->display->insertPlainText(QStringLiteral("the CT2Toumo matrix is:\n"));
	t=m_transformation->getYXZRotationAngles(Toumo2CTMatrix);
	printToUI(t);
	vtkMatrix4x4 *transformedCTMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4::Multiply4x4(CT2MarkerMatrix, Toumo2CTMatrix, transformedCTMatrix);
	vtkMatrix4x4::Multiply4x4(m_MarkerActor->GetMatrix(),transformedCTMatrix, transformedCTMatrix);
	m_CTActor->SetUserMatrix(transformedCTMatrix);
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::on_ct2axisButton_clicked() {
	vtkMatrix4x4* CT2MarkerMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4* backCTmatrix = vtkMatrix4x4::New();
	getCoorsInMarker(m_MarkerActor->GetMatrix(), m_CTActor->GetMatrix(), CT2MarkerMatrix);
	CT2MarkerMatrix->Invert();
	vtkMatrix4x4::Multiply4x4(CT2MarkerMatrix, m_CTActor->GetMatrix(), backCTmatrix);
	m_CTActor->SetUserMatrix(backCTmatrix);
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}
void registrationWidget::getMatrixTransAtoB(vtkMatrix4x4 *A, vtkMatrix4x4 *B, vtkMatrix4x4 *trans) {
	vtkSmartPointer<vtkMatrix4x4> invertA = vtkSmartPointer<vtkMatrix4x4>::New();
	vtkMatrix4x4::Invert(A, invertA);
	vtkMatrix4x4::Multiply4x4(invertA, B, trans);
}
void registrationWidget::getCoorsInMarker(vtkMatrix4x4 *marker, vtkMatrix4x4 *model, vtkMatrix4x4 *trans)
{
		ui->display->insertPlainText(QStringLiteral("the angles of model in marker coordiates:\n"));
		getMatrixTransAtoB(marker, model, trans);
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
	size_t index0 = originS.find(DELIMITER);
	string sub0 = originS.substr(0, index0);
	array[0] = atof(sub0.c_str());
	string sub00 = originS.substr(index0 + 1,originS.size());
	size_t index1 = sub00.find(DELIMITER);
	string sub1 = sub00.substr(0, index1);
	array[1] = atof(sub1.c_str());
	string sub11 = sub00.substr(index1+1, sub00.size());
	size_t index2 = sub11.find(DELIMITER);
	string sub2 = sub11.substr(0, index2);
	array[2] = atof(sub2.c_str());
	string sub22 = sub11.substr(index2 + 1, sub11.size());
	size_t index3 = sub22.find(DELIMITER);
	string sub3 =sub22.substr(0,index3);
	array[3] = atof(sub3.c_str());
}
void registrationWidget::printToUI(std::string t) {
	ui->display->insertPlainText(QString::fromStdString(t));
	ui->display->insertPlainText("\n");
}