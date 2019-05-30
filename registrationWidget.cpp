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
#include <fstream>
#include "vtkSTLWriter.h"
#define LINE_LEN 100.
#define PI 3.14159265
std::string s_dicomName = "E:\\VTK_Project\\registration_test\\data\\DICOM\\20190408\\10470000";
std::string s_modelName = "E:\\AR\\Vessel\\registration\\Assets\\10470000y180.obj";
std::string s_stlName = "C:\\Users\\29477\\Desktop\\registrationTest.stl";
registrationWidget::registrationWidget(QWidget *parent)
	: QWidget(parent)
{
	ui = new Ui_Form;
	ui->setupUi(this);
	connect(ui->present, SIGNAL(clicked()), this, SLOT(setCurrentTransformation()));
	connect(ui->transform, SIGNAL(clicked()), this, SLOT(mapCT2Toumo()));
	connect(ui->ct2axis, SIGNAL(clicked()), this, SLOT(mapCT2Marker()));
	this->readCase(s_stlName);
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(s_modelName);
	m_CTActor->SetMapper(m_CTMapper);
	m_ToumoOriginActor->SetMapper(m_ToumoMapper);
	//m_CTActor->SetMapper(imageMapper);
	//m_ToumoOriginActor->SetMapper(imageMapper);
	m_CTActor->GetProperty()->SetColor(1, 1, 0.9412);
	m_CTcenter = m_CTActor->GetCenter();
	m_CTorigin2center->Translate(m_CTcenter[0], m_CTcenter[1], m_CTcenter[2]);
	m_CTcenter2origin->Invert(m_CTorigin2center->GetMatrix(), m_CTcenter2origin);
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
	vtkTransform *transformation = vtkTransform::New();
	renderer->AddActor(m_coneActor);
	renderer->AddActor(m_worldActor);
	renderer->AddActor(m_CTActor);
	renderer->AddActor(m_ToumoOriginActor);
	renderer->AddActor(m_MarkerActor);
	renderer->SetBackground(.3, .3, .5);
	m_renderWindowInteractor->Initialize();
	transformation->Delete();

	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix->SetElement(0, 0, 0.638);
	matrix->SetElement(0, 1, -0.77);
	matrix->SetElement(0, 2, 0.026);
	matrix->SetElement(0, 3, -0.12188);
	matrix->SetElement(2, 0, -0.78);
	matrix->SetElement(2, 1, -0.646);
	matrix->SetElement(2, 2, 0.006);
	matrix->SetElement(2, 3, 0.841138);
	matrix->SetElement(1, 0, -0.012);
	matrix->SetElement(1, 1, 0.028);
	matrix->SetElement(1, 2, 0.999);
	matrix->SetElement(1, 3, -1.23499);
	matrix->SetElement(3, 0, 0);
	matrix->SetElement(3, 1, 0);
	matrix->SetElement(3, 2, 0);
	matrix->SetElement(3, 3, 1);

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
	vtkTransform *transformation = vtkTransform::New();
	vtkMatrix4x4 *transMatrix = vtkMatrix4x4::New();
	transformation->Identity();
	transformation->Translate(x, y, z);
	transformation->RotateZ(rz);
	transformation->RotateX(rx);
	transformation->RotateY(ry);
	transMatrix = transformation->GetMatrix();
	transMatrix->Print(cout);
	return transMatrix;
}
void registrationWidget::transToumo(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = setTransformation_right(x, y, z, rx, ry, rz);
	m_Toumomatrix = m_ToumoOriginActor->GetMatrix();
	vtkTransform *transformation = vtkTransform::New();
	transformation->SetMatrix(matrix);
	m_ToumoOriginActor->SetUserTransform(transformation);
	getXYZRotationAngles(matrix);
	writeOBJCase();
}
void registrationWidget::transMarker(const float x, const float y, const float z, const float rx, const float ry, const float rz){
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = setTransformation_right(x, y, z, rx, ry, rz);
	vtkTransform *transformation = vtkTransform::New();
	transformation->SetMatrix(matrix);
	m_MarkerActor->SetUserTransform(transformation);
}
void registrationWidget::setTransformation_left(const int flag, const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	transToumo(x, y, -z, -rx, -ry, rz);
	transMarker(x, y, -z, -rx, -ry, rz);
}

void registrationWidget::setCurrentTransformation()
{
	vtkMatrix4x4 *ToumoInitial = vtkMatrix4x4::New();
	transToumo(100, 200, 300, 60, 0, 0);
	transMarker(.1, .2, .3, 10, 20, 30);
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
		CTtrans->Print(cout);
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
}

void registrationWidget::readCase(std::string fileName)
{
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(fileName);
	cout << "the extension of the filename: ";
	cout << extension << endl;
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
