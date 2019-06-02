/*
*Author: Chen Ying
*Version:1.3
*Date: 2019-4-17
*History change description:
change to a universal model whose origin and center are not in the same place;
*/
#include "registrationWidget.h"
#include <vtkProperty.h>
#include <vtkColor.h>
#include <vtkOBJExporter.h>
#include <vtkDecimatePro.h>
#include <vtksys/SystemTools.hxx>
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkNIFTIImageReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageShrink3D.h>
#include <vtkMarchingCubes.h>
#include <vtkRenderer.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "vtkSTLWriter.h"
#define LINE_LEN 100.
#define PI 3.14159265
std::string s_dicomName = "E:\\VTK_Project\\registration_test\\data\\DICOM\\20190408\\10470000";
std::string s_modelName = "E:\\AR\\Vessel\\registration\\Assets\\10470000origin.obj";
std::string s_stlName = "E:\\VTK\\VTK-8.1.1_build\\ExternalData\\Testing\\Data\\42400-IDGH.stl";
registrationWidget::registrationWidget(QWidget *parent)
	: QWidget(parent)
{
	ui = new Ui_Form;
	ui->setupUi(this);
	connect(ui->present, SIGNAL(clicked()), this, SLOT(setPresentStates()));
	connect(ui->transform, SIGNAL(clicked()), this, SLOT(mapCT2Toumo()));
	connect(ui->ct2axis, SIGNAL(clicked()), this, SLOT(mapCT2Marker()));
	this->readCase(s_stlName);
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(s_stlName);
	if(extension==" "){
		m_CTActor->SetMapper(imageMapper);
	    m_ToumoOriginActor->SetMapper(imageMapper);
	}
	else {
		m_CTActor->SetMapper(m_CTMapper);
		m_ToumoOriginActor->SetMapper(m_ToumoMapper);
	}	
	m_CTActor->GetProperty()->SetColor(1, 1, 0.9412);
	m_CTcenter = m_CTActor->GetCenter();
	m_CTorigin2center->Translate(m_CTcenter[0], m_CTcenter[1], m_CTcenter[2]);
	m_CTcenter2origin->Invert(m_CTorigin2center->GetMatrix(), m_CTcenter2origin);
	m_ToumoOriginActor->GetProperty()->SetOpacity(0.9);
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
	vtkRenderer *renderer = vtkRenderer::New();
	m_renderWindow->SetSize(800, 800);
	m_renderWindow->AddRenderer(renderer);
	m_renderWindowInteractor->SetRenderWindow(m_renderWindow);
	m_conedata->SetAngle(30);
	m_conedata->SetHeight(.2);
	m_conedata->SetRadius(.1);
	m_conedata->SetResolution(10);
	m_coneMapper->SetInputConnection(m_conedata->GetOutputPort());
	m_coneActor->SetMapper(m_coneMapper);
	renderer->AddActor(m_coneActor);
	renderer->AddActor(m_worldActor);
	renderer->AddActor(m_CTActor);
	renderer->AddActor(m_ToumoOriginActor);
	renderer->AddActor(m_MarkerActor);
	renderer->SetBackground(.3, .3, .5);
	m_renderWindowInteractor->Initialize();
	char Toumo[4][4] = {0.6380,  -0.770,  0.026,  -0.12188,
						-0.012,   0.028,  0.999,  -1.23499,
						-0.780,  -0.646, -0.006,  0.841138,
							 0,	      0,      0,         1
	};
	vtkMatrix4x4 *matrix = setCurrentMatrix(Toumo);
	vtkMatrix4x4 *rotationM = vtkMatrix4x4::New();
	rotationM->Zero();
	rotationM->SetElement(0, 0, -1);
	rotationM->SetElement(1, 2, -1);
	rotationM->SetElement(2, 1, -1);
	rotationM->SetElement(3, 3, 1);
	vtkMatrix4x4 *newMatrix = vtkMatrix4x4::New();
	newMatrix->Multiply4x4(rotationM, matrix, newMatrix);
	//getXYZRotationAngles(newMatrix);	


}

registrationWidget::~registrationWidget()
{
	m_CTMapper->Delete();
	m_CTActor->Delete();
	m_ToumoMapper->Delete();
	m_ToumoOriginActor->Delete();
	m_MarkerActor->Delete();
	m_worldActor->Delete();
	m_conedata->Delete();
	m_coneMapper->Delete();
	m_coneActor->Delete();
	m_renderWindow->Delete();
	m_renderWindowInteractor->Delete();
}

vtkMatrix4x4 * registrationWidget::setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	std::ofstream outTXT("C:\\Users\\cy\\Desktop\\matrix.txt", ios::out);
	vtkTransform *transformation = vtkTransform::New();
	vtkMatrix4x4 *transMatrix = vtkMatrix4x4::New();
	transformation->Identity();
	transformation->Translate(x, y, z);
	transformation->RotateZ(rz);
	transformation->RotateX(rx);
	transformation->RotateY(ry);
	transMatrix = transformation->GetMatrix();
	transMatrix->Print(outTXT);
	outTXT.close();
	return transMatrix;
}

vtkMatrix4x4 * registrationWidget::setTransformation_left(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = setTransformation_right(-x, y, z, rx, -ry, -rz);
	return matrix;
}

void registrationWidget::transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = setTransformation_left(x, y, z, rx, ry, rz);
	vtkTransform *transformation = vtkTransform::New();
	transformation->SetMatrix(matrix);
	m_ToumoOriginActor->SetUserTransform(transformation);
	transformation->Delete();
	getXYZRotationAngles(matrix);
	//writeOBJCase();
}

vtkMatrix4x4 *  registrationWidget::objTrans(vtkMatrix4x4 *m){
	//In VTK, the obj model should firstly rotate around X with 90 degree, then back-direct the x axis
	
	//to do the back-direct
	m->SetElement(0, 1, -m->GetElement(0, 1));
	m->SetElement(0, 2, -m->GetElement(0, 2));
	m->SetElement(1, 2, -m->GetElement(1, 2));
	m->SetElement(2, 1, -m->GetElement(2, 1));
	//now, set the translation along x axis with negative translation values
	m->SetElement(0, 3, -m->GetElement(0, 3));
	return m;
}

void registrationWidget::transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz){
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = setTransformation_right(x, y, z, rx, ry, rz);
	vtkTransform *transformation = vtkTransform::New();
	transformation->SetMatrix(matrix);
	m_MarkerActor->SetUserTransform(transformation);
}

void registrationWidget::setPresentStates()
{
	QString qrx=ui->rx->toPlainText(); 
	float rx = atof(qrx.toStdString().c_str());
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

	vtkMatrix4x4 *ToumoInitial = vtkMatrix4x4::New();
	transToumo(100*x, 100*y, 100*z, rx, ry, rz);
	//transMarker(.1, .2, .3, 10, 20, 30);
	m_renderWindow->Render();
	m_renderWindowInteractor->Start();
}

vtkMatrix4x4 * registrationWidget::setCurrentMatrix(char m[4][4]) {
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix->SetElement(0, 0, m[0][0]);
	matrix->SetElement(0, 1, m[0][1]);
	matrix->SetElement(0, 2, m[0][2]);
	matrix->SetElement(0, 3, m[0][3]);
	matrix->SetElement(1, 0, m[1][0]);
	matrix->SetElement(1, 1, m[1][1]);
	matrix->SetElement(1, 2, m[1][2]);
	matrix->SetElement(1, 3, m[1][3]);
	matrix->SetElement(2, 0, m[2][0]);
	matrix->SetElement(2, 1, m[2][1]);
	matrix->SetElement(2, 2, m[2][2]);
	matrix->SetElement(2, 3, m[2][3]);
	matrix->SetElement(3, 0, m[3][0]);
	matrix->SetElement(3, 1, m[3][1]);
	matrix->SetElement(3, 2, m[3][2]);
	matrix->SetElement(3, 3, m[3][3]);
	return matrix;
}

void registrationWidget::mapCT2Toumo()
{
	std::ofstream outTXT("C:\\Users\\cy\\Desktop\\matrix.txt", ios::out);
	m_Markermatrix = m_MarkerActor->GetMatrix();
	m_CToriginmatrix = m_CTActor->GetMatrix();
	m_Toumomatrix = m_ToumoOriginActor->GetMatrix();
#pragma region test
	vtkMatrix4x4 *world2MarkerMatrix = m_Markermatrix;
	printf("The marker pistion and rotation:\n");
	world2MarkerMatrix->Print(outTXT);
	getXYZRotationAngles(world2MarkerMatrix);
	vtkMatrix4x4 *Marker2ToumoMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4 *Marker2worldMatrix = vtkMatrix4x4::New();
	Marker2worldMatrix->Invert(world2MarkerMatrix, Marker2worldMatrix);
	vtkMatrix4x4 *Marker2CTMatrix = vtkMatrix4x4::New();
	Marker2CTMatrix->Multiply4x4(m_CToriginmatrix, Marker2worldMatrix, Marker2CTMatrix);
	printf("The CT in marker coordinates:\n");
	Marker2CTMatrix->Print(outTXT);
	getXYZRotationAngles(Marker2CTMatrix);
	Marker2ToumoMatrix->Multiply4x4(m_Toumomatrix, Marker2worldMatrix, Marker2ToumoMatrix);
	printf("The Toumo in marker coordinates:\n");
	Marker2ToumoMatrix->Print(outTXT);
	getXYZRotationAngles(Marker2ToumoMatrix);
	vtkMatrix4x4 *CT2MarkerMatrix = vtkMatrix4x4::New();
	CT2MarkerMatrix->Invert(Marker2CTMatrix, CT2MarkerMatrix);
	vtkMatrix4x4 *CT2ToumoMatrix = vtkMatrix4x4::New();
	CT2ToumoMatrix->Multiply4x4(Marker2ToumoMatrix, CT2MarkerMatrix, CT2ToumoMatrix);
	printf("the CT2Toumo matrix is:\n");
	CT2ToumoMatrix->Print(outTXT);
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

	outTXT.close();
}

void registrationWidget::mapCT2Marker()
{
	printf("\n%d\n", firstClick2Marker);
	std::ofstream outTXT("C:\\Users\\cy\\Desktop\\matrix.txt", ios::out);
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
		CTtrans->Print(outTXT);
		getXYZRotationAngles(CTtrans);
		vtkTransform *trans = vtkTransform::New();
		trans->SetMatrix(CTtrans);
		if(polyOrImage==1)
			trans->Scale(0.001, 0.001, 0.001);
		m_CTActor->SetUserTransform(trans);
		m_renderWindow->Render();
		m_renderWindowInteractor->Start();
	}
	firstClick2Marker = 1;
	outTXT.close();
}

void registrationWidget::readCase(std::string fileName)
{
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(fileName);
	ui->display->insertPlainText("the extension of the filename: ");
	ui->display->insertPlainText(QString::fromStdString(extension));
	ui->display->insertPlainText(QString::fromStdString("\n"));
	vtkAlgorithm *reader = vtkAlgorithm::New();
	vtkAbstractPolyDataReader *polyreader;
	vtkDICOMImageReader *dicomReader = vtkDICOMImageReader::New();
	vtkPolyData *CTdata = vtkPolyData::New();
	vtkPolyData *Toumodata = vtkPolyData::New();
	vtkNIFTIImageReader *nifftiReader = vtkNIFTIImageReader::New();
	if (extension == ".obj")
	{
		polyreader = vtkOBJReader::New();
		polyreader->SetFileName(fileName.data());
		polyreader->Update();
		CTdata = polyreader->GetOutput();
		Toumodata = polyreader->GetOutput();
		m_CTMapper->SetInputData(CTdata);
		m_ToumoMapper->SetInputData(Toumodata);
	}
	else if (extension == ".stl")
	{
		polyreader = vtkSTLReader::New();
		polyreader->SetFileName(fileName.data());
		polyreader->Update();
		CTdata = polyreader->GetOutput();
		Toumodata = polyreader->GetOutput();
		m_CTMapper->SetInputData(CTdata);
		m_ToumoMapper->SetInputData(Toumodata);
	}
	else if (extension == ".nii")
	{
		nifftiReader = vtkNIFTIImageReader::New();
		nifftiReader->SetFileName(fileName.data());
		nifftiReader->Update();
		niftiImage = nifftiReader->GetOutput();
		imageMapper->SetInputData(niftiImage);
	}
	else
	{
		dicomReader->SetDataByteOrderToLittleEndian();
		dicomReader->SetDirectoryName(fileName.data());
		dicomReader->ReleaseDataFlagOn();
		dicomReader->SetDataScalarTypeToUnsignedShort();
		dicomReader->Update();
		dicomReader->GlobalWarningDisplayOff();
		vtkImageShrink3D * shrink= vtkImageShrink3D::New();
		vtkMarchingCubes * CTdata= vtkMarchingCubes::New();
		vtkDecimatePro *decimate = vtkDecimatePro::New();
		vtkImageData *imageData = dicomReader->GetOutput();
		shrink->SetInputData((vtkDataObject*)imageData);
		shrink->SetShrinkFactors(1, 1, 1);
		shrink->AveragingOn();
		shrink->Update();
		CTdata->SetInputData((vtkDataSet*)shrink->GetOutput());
		CTdata->ComputeNormalsOn();
		CTdata->SetValue(0, -150);
		CTdata->Update();
		decimate->SetInputData(CTdata->GetOutputDataObject(0));
		decimate->SetTargetReduction(0.8);
		decimate->Update();
		vtkSmartPointer<vtkTransform> rescaleCT = vtkSmartPointer<vtkTransform>::New();
		rescaleCT->Scale(0.001, 0.001, 0.001);
		m_CTActor->SetUserTransform(rescaleCT);
		m_ToumoOriginActor->SetUserTransform(rescaleCT);
		m_CTMapper->SetInputConnection(decimate->GetOutputPort());
		m_ToumoMapper->SetInputConnection(decimate->GetOutputPort());
		polyOrImage = 1;
	}	
}

void registrationWidget::writeOBJCase()
{
	vtkRenderer *exportrenderer = vtkRenderer::New();
	exportrenderer->AddActor(m_ToumoOriginActor);
	m_exportWindow->AddRenderer(exportrenderer);
	vtkOBJExporter *objExporter = vtkOBJExporter::New();
	objExporter->SetFilePrefix("registrationTest");
	objExporter->SetInput(m_exportWindow);
	objExporter->Write();
}

void registrationWidget::writeSTLCase(vtkMarchingCubes *data)
{
	vtkSTLWriter *stlExporter = vtkSTLWriter::New();
	//stlExporter->SetInputData(data);
	stlExporter->SetInputConnection(data->GetOutputPort());
	stlExporter->SetFileTypeToBinary();
	stlExporter->SetFileName("C:\\Users\\29477\\Desktop\\registrationTest.stl");
	stlExporter->Update();
	stlExporter->Write();
}

void registrationWidget::getXYZRotationAngles(vtkMatrix4x4 *matrix)
{
	double y1 = 0, x1 = 0, z1 = 0;
	y1 = atan2(-matrix->GetElement(2, 0), matrix->GetElement(2, 2))*180.0 / PI;
	x1 = atan2(matrix->GetElement(2, 1), sqrt(matrix->GetElement(2, 0)*matrix->GetElement(2, 0) + matrix->GetElement(2, 2)*matrix->GetElement(2, 2)))*180.0 / PI;
	z1 = atan2(-matrix->GetElement(0, 1), matrix->GetElement(1, 1))*180.0 / PI;
	printf("%lf, %lf, %lf", x1, y1, z1);
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