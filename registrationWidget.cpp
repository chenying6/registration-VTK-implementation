/*
*Copyright(C), 2019, Remebot
*Author: Chen Ying
*Version:1.1
*Date: 2019-4-3
*History change description:
将原先center与origin重合的特殊模型改为通用模型，即直接由VTK生成的obj格式文件
添加由旋转矩阵得各轴向旋转角度的函数
*/
#include "registrationWidget.h"
#include "ui_transparency.h"
#include <vtkProperty.h>
#include <vtkColor.h>
#include <vtkOBJExporter.h>
#include <vtkDecimatePro.h>
#include <math.h>
#include <fstream>
#define LINE_LEN 1.
#define PI 3.14159265
std::string s_fileName = "E:\\Remebot\\dcm";
std::string s_modelName = "C:\\Users\\29477\\Desktop\\Decimate.obj";

registrationWidget::registrationWidget(QWidget *parent)
	: QWidget(parent)
{
	ui = new Ui_Form;
	ui->setupUi(this);
	connect(ui->present, SIGNAL(clicked()), this, SLOT(setCurrentTransformation()));
	connect(ui->transform, SIGNAL(clicked()), this, SLOT(mapCT2Toumo()));
	connect(ui->ct2axis, SIGNAL(clicked()), this, SLOT(mapCT2Marker()));
	m_reader->SetFileName(s_modelName.data());
	m_reader->Update();
	m_CTdata = m_reader->GetOutput();
	m_Toumodata = m_reader->GetOutput();
	//this->readCase();
	m_CTMapper->SetInputData(m_CTdata);
	m_ToumoMapper->SetInputData(m_Toumodata);
	m_CTActor->SetMapper(m_CTMapper);
	m_CTActor->GetProperty()->SetColor(1, 1, 0.9412);
	m_CTcenter = m_CTActor->GetCenter();
	m_CTorigin2center->Translate(m_CTcenter[0], m_CTcenter[1], m_CTcenter[2]);
	m_CTcenter2origin->Invert(m_CTorigin2center->GetMatrix(), m_CTcenter2origin);
	m_ToumoOriginActor->SetMapper(m_ToumoMapper);
	m_ToumoOriginActor->GetProperty()->SetOpacity(0.5);
	m_ToumoOriginActor->GetProperty()->SetColor(0.1, 1, 0.1);
	m_MarkerActor->SetTotalLength(LINE_LEN, LINE_LEN, LINE_LEN);
	m_MarkerActor->SetShaftType(0);
	m_MarkerActor->SetAxisLabels(0);
	m_MarkerActor->SetCylinderRadius(0.02);
	m_worldActor->SetPosition(0, 0, 0);
	m_worldActor->SetTotalLength(LINE_LEN, LINE_LEN, LINE_LEN);
	m_worldActor->SetShaftType(0);
	m_worldActor->SetAxisLabels(0);
	m_worldActor->SetCylinderRadius(0.02);
	m_renderWindow->SetSize(800, 800);
	m_renderWindow->AddRenderer(m_renderer);
	m_renderWindowInteractor->SetRenderWindow(m_renderWindow);
	m_conedata->SetAngle(30);
	m_conedata->SetHeight(.2);
	m_conedata->SetRadius(.1);
	m_conedata->SetResolution(10);
	m_coneMapper->SetInputConnection(m_conedata->GetOutputPort());
	m_coneActor->SetMapper(m_coneMapper);
	m_exportrenderer->AddActor(m_CTActor);
	m_renderer->AddActor(m_coneActor);
	m_renderer->AddActor(m_worldActor);
	m_renderer->AddActor(m_CTActor);
	m_renderer->AddActor(m_ToumoOriginActor);
	m_renderer->AddActor(m_MarkerActor);
	m_renderer->SetBackground(.3, .3, .5);
	m_exportWindow->AddRenderer(m_exportrenderer);
	m_renderWindowInteractor->Initialize();
	//writeOBJCase();
}

registrationWidget::~registrationWidget()
{
	m_reader->Delete();
	m_CTdata->Delete();
	m_CTMapper->Delete();
	m_CTActor->Delete();
	m_Toumodata->Delete();
	m_ToumoMapper->Delete();
	m_ToumoOriginActor->Delete();
	m_MarkerActor->Delete();
	m_worldActor->Delete();
	m_conedata->Delete();
	m_coneMapper->Delete();
	m_coneActor->Delete();
	m_renderer->Delete();
	m_renderWindow->Delete();
	m_renderWindowInteractor->Delete();
}

void registrationWidget::setTransformation_right(const int flag, const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkTransform *transformation = vtkTransform::New();
	vtkTransform *conetrans = vtkTransform::New();
	vtkMatrix4x4 *coneMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4 *transMatrix = vtkMatrix4x4::New();
	transformation->Identity();
	transformation->PostMultiply();
	transformation->RotateWXYZ(ry, 0, 1, 0);
	transformation->RotateWXYZ(rx, 1, 0, 0);
	transformation->RotateWXYZ(rz, 0, 0, 1);
	transformation->Translate(x, y, z);
	transMatrix = transformation->GetMatrix();
	switch (flag) {
	case 1:
		coneMatrix = transMatrix;
		conetrans->SetMatrix(coneMatrix);
		m_coneActor->SetUserTransform(conetrans);
		m_Toumomatrix = m_ToumoOriginActor->GetMatrix();
		transMatrix->Multiply4x4(m_CTcenter2origin, transMatrix, transMatrix);
		transMatrix->Multiply4x4(transMatrix, m_Toumomatrix, transMatrix);
		transformation->SetMatrix(transMatrix);
		m_ToumoOriginActor->SetUserTransform(transformation);
		break;
	case 2:
		transformation->SetMatrix(transMatrix);
		m_MarkerActor->SetUserTransform(transformation);
		break;
	}
	transformation = NULL;
	coneMatrix = NULL;
	conetrans = NULL;
	transMatrix = NULL;
}

void registrationWidget::setTransformation_left(const int flag, const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	setTransformation_right(flag, x, y, -z, -rx, -ry, rz);
}

void registrationWidget::setCurrentTransformation()
{
	this->setTransformation_right(1, .1, .2, .3, 60, 0, 60);
	this->setTransformation_right(2, .3, .1, .2, 55, 35, 15);
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

void registrationWidget::mapCT2Toumo()
{	
	m_Markermatrix = m_MarkerActor->GetMatrix();
	m_CToriginmatrix = m_CTActor->GetMatrix();
	m_Toumomatrix = m_ToumoOriginActor->GetMatrix();
#pragma region test
	vtkMatrix4x4 *world2MarkerMatrix = m_Markermatrix;
	printf("The marker pistion and rotation:\n");
	world2MarkerMatrix->Print(cout);
	getXYZRotationAngles(world2MarkerMatrix);
	vtkMatrix4x4 *Marker2ToumoMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4 *Marker2worldMatrix = vtkMatrix4x4::New();
	Marker2worldMatrix->Invert(world2MarkerMatrix, Marker2worldMatrix);
	vtkMatrix4x4 *Marker2CTMatrix = vtkMatrix4x4::New();
	Marker2CTMatrix->Multiply4x4(m_CToriginmatrix, Marker2worldMatrix, Marker2CTMatrix);
	printf("The CT in marker coordinates:\n");
	Marker2CTMatrix->Print(cout);
	getXYZRotationAngles(Marker2CTMatrix);
	Marker2ToumoMatrix->Multiply4x4(m_Toumomatrix, Marker2worldMatrix, Marker2ToumoMatrix);
	printf("The Toumo in marker coordinates:\n");
	Marker2ToumoMatrix->Print(cout);
	getXYZRotationAngles(Marker2ToumoMatrix);
	vtkMatrix4x4 *CT2MarkerMatrix = vtkMatrix4x4::New();
	CT2MarkerMatrix->Invert(Marker2CTMatrix, CT2MarkerMatrix);
	vtkMatrix4x4 *CT2ToumoMatrix = vtkMatrix4x4::New();
	CT2ToumoMatrix->Multiply4x4(Marker2ToumoMatrix, CT2MarkerMatrix, CT2ToumoMatrix);
	printf("the CT2Toumo matrix is:\n");
	CT2ToumoMatrix->Print(cout);
	getXYZRotationAngles(CT2ToumoMatrix);
	vtkMatrix4x4 *transformedCTMatrix = vtkMatrix4x4::New();
	transformedCTMatrix->Multiply4x4(CT2ToumoMatrix, Marker2CTMatrix, transformedCTMatrix);
	transformedCTMatrix->Multiply4x4(transformedCTMatrix, world2MarkerMatrix, transformedCTMatrix);
	vtkTransform *CTtrans = vtkTransform::New();
	CTtrans->SetMatrix(transformedCTMatrix);
	m_CTActor->SetUserTransform(CTtrans);
#pragma endregion test
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();

	ofstream outFile("matrix.txt", fstream::trunc);
	for (int i = 0; i < index; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			outFile<<outMatrix[i][j]<<" ";
		}
		outFile << endl;
	}
	outFile.close();
}

void registrationWidget::mapCT2Marker()
{
	printf("\n%d\n", firstClick2Marker);
	if (firstClick2Marker == 0)
	{
		printf("%lf,%lf,%lf", m_CTcenter[0], m_CTcenter[1], m_CTcenter[2]);
		m_Markermatrix = m_MarkerActor->GetMatrix();
		m_CToriginmatrix = m_CTActor->GetMatrix();
		vtkMatrix4x4 *inverseMarkerMatrix = vtkMatrix4x4::New();//标记板坐标的逆矩阵
		inverseMarkerMatrix->Invert(m_Markermatrix, inverseMarkerMatrix);
		vtkMatrix4x4 *CTmarker2originMatrix = vtkMatrix4x4::New();
		CTmarker2originMatrix->Multiply4x4(m_CToriginmatrix, inverseMarkerMatrix, CTmarker2originMatrix);
		vtkMatrix4x4 *CTorigin2markerMatrix = vtkMatrix4x4::New();
		CTorigin2markerMatrix->Invert(CTmarker2originMatrix, CTorigin2markerMatrix);
		vtkMatrix4x4 *CTtrans = vtkMatrix4x4::New();
		CTtrans->Multiply4x4(CTorigin2markerMatrix, m_CToriginmatrix, CTtrans);
		CTtrans->Multiply4x4(CTtrans, m_CTcenter2origin, CTtrans);
		vtkTransform *trans = vtkTransform::New();
		trans->SetMatrix(CTtrans);
		//trans->Scale(0.001, 0.001, 0.001);
		m_CTActor->SetUserTransform(trans);
		m_renderWindow->Render();
		m_renderWindowInteractor->Start();
	}
	firstClick2Marker = 1;
}

//void registrationWidget::readCase()
//{
//	this->m_reader = vtkDICOMImageReader::New();
//	this->m_reader->SetDataByteOrderToLittleEndian();
//	this->m_reader->SetDirectoryName(s_fileName.data());
//	this->m_reader->ReleaseDataFlagOn();
//	this->m_reader->SetDataScalarTypeToUnsignedShort();
//	this->m_reader->Update();
//	this->m_reader->GlobalWarningDisplayOff();
//	vtkImageData *imageData = m_reader->GetOutput();
//	this->m_shrink = vtkImageShrink3D::New();
//	this->m_shrink->SetInputData((vtkDataObject*)imageData);
//	this->m_shrink->SetShrinkFactors(1, 1, 1);
//	this->m_shrink->AveragingOn();
//	this->m_shrink->Update();
//	this->m_CTdata = vtkMarchingCubes::New();
//	this->m_CTdata->SetInputData((vtkDataSet*)this->m_shrink->GetOutput());
//	this->m_CTdata->ComputeNormalsOn();
//	this->m_CTdata->SetValue(0, -150);
//	this->m_CTdata->Update();
//	this->m_Toumodata = vtkMarchingCubes::New();
//	this->m_Toumodata->SetInputData((vtkDataSet*)this->m_shrink->GetOutput());
//	this->m_Toumodata->ComputeNormalsOn();
//	this->m_Toumodata->SetValue(0, -150);
//	this->m_Toumodata->Update();
//	vtkSmartPointer<vtkTransform> rescaleCT = vtkSmartPointer<vtkTransform>::New();
//	rescaleCT->Scale(0.001, 0.001, 0.001);
//	m_CTActor->SetUserTransform(rescaleCT);
//	vtkSmartPointer<vtkTransform> rescaleToumo = vtkSmartPointer<vtkTransform>::New();
//	rescaleToumo->Scale(0.001, 0.001, 0.001);
//	m_ToumoActor->SetUserTransform(rescaleToumo);
//	decimate->SetInputData(m_CTdata);
//	decimate->SetTargetReduction(0.8);
//	decimate->Update();
//	m_CTMapper->SetInputConnection(decimate->GetOutputPort());
//	m_ToumoMapper->SetInputConnection(decimate->GetOutputPort());
//}

void registrationWidget::writeOBJCase()
{
	vtkOBJExporter *objExporter = vtkOBJExporter::New();
	objExporter->SetFilePrefix("registrationTest");
	objExporter->SetInput(m_exportWindow);
	objExporter->Write();
}

void registrationWidget::getXYZRotationAngles(vtkMatrix4x4 *m)
{
	double anglealpha, anglebeta, anglegamma;
	double radianalpha, radianbeta, radiangamma;
	anglealpha = asin(m->GetElement(2, 1))*180.0 / PI;
	radianalpha = asin(m->GetElement(2, 1));
	anglebeta = -asin((m->GetElement(2, 0)*1.0 / cos(radianalpha)))*180.0 / PI;
	anglegamma = -asin((m->GetElement(0, 1)*1.0 / cos(radianalpha)))*180.0 / PI;
	radianbeta = anglebeta * 1.0 / 180.0 * PI;
	radiangamma = anglegamma * 1.0 / 180.0 * PI;
	cout << "the angles are:" << endl;
	printf("%lf,%lf,%lf\n", anglealpha, anglebeta, anglegamma);
}

void registrationWidget::outputMatrix(vtkMatrix4x4 *m)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			outMatrix[index][j] = m->GetElement(i, j);
		}
		index++;
	}
}

int main(int argc, char *argv[]) {
	QApplication app(argc, argv);
	registrationWidget registrationWidget;
	registrationWidget.show();
	return app.exec();
}
