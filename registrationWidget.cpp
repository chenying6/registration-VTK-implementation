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
#include <vtkCylinderSource.h>
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkPolygon.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "vtkSTLWriter.h"
#define LINE_LEN 1.
#define PI 3.14159265
std::string s_dicomName = "E:\\VTK_Project\\registration_test\\data\\DICOM\\20190408\\10470000";
std::string s_modelName = "E:\\VTK_Project\\registration_test\\build\\zidingyi.obj";
std::string s_stlName = "C:\\Users\\29477\\Desktop\\registrationTest.stl";
registrationWidget::registrationWidget(QWidget *parent)
	: QWidget(parent)
{
	ui = new Ui_Form;
	ui->setupUi(this);
	connect(ui->present, SIGNAL(clicked()), this, SLOT(setPresentStates()));
	connect(ui->transform, SIGNAL(clicked()), this, SLOT(mapCT2Toumo()));
	connect(ui->ct2axis, SIGNAL(clicked()), this, SLOT(mapCT2Marker()));
	this->readCase(s_modelName);
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
#pragma region test_back_transformation
	vtkMatrix4x4 *initialTransformation = vtkMatrix4x4::New();
	initialTransformation = setTransformation_right(10, 10, 10, 10, 10, 10);
	m_ToumoOriginActor->SetUserMatrix(initialTransformation);
	m_Toumomatrix = m_ToumoOriginActor->GetUserMatrix();
	vtkMatrix4x4 *initialInvertMatrix = vtkMatrix4x4::New();
	m_Toumomatrix->Print(cout);
	initialInvertMatrix->Invert(m_Toumomatrix, initialInvertMatrix);
	initialInvertMatrix->Print(cout);
	getXYZRotationAngles(initialInvertMatrix);
#pragma endregion test_back_transformation
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
	prepareExportCTModel();
	renderer->AddActor(m_worldActor);
	renderer->AddActor(m_CTActor);
	//renderer->AddActor(m_ToumoOriginActor);
	renderer->AddActor(m_MarkerActor);
	renderer->SetBackground(.3, .3, .5);
	m_renderWindowInteractor->SetRenderWindow(ui->qvtkWidget->GetRenderWindow());
	ui->qvtkWidget->SetRenderWindow(m_renderWindow);
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
	m_renderWindow->Delete();
	m_renderWindowInteractor->Delete();
}

vtkMatrix4x4 * registrationWidget::setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	std::ofstream outTXT("C:\\Users\\29477\\Desktop\\matrix.txt", ios::out);
	vtkTransform *transformation = vtkTransform::New();
	vtkMatrix4x4 *transMatrix = vtkMatrix4x4::New();
	transformation->Identity();
	transformation->Translate(x, y, z);
	transformation->RotateY(ry);
	transformation->RotateX(rx);
	transformation->RotateZ(rz);
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
	m_Toumomatrix = m_ToumoOriginActor->GetUserMatrix();
	matrix->Multiply4x4(matrix, m_Toumomatrix, matrix);
	m_ToumoOriginActor->SetUserMatrix(matrix);
	ui->display->insertPlainText("the present transformation of Toumo:\n");
	getXYZRotationAngles(matrix);
}
void registrationWidget::transCT(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix = setTransformation_left(x, y, z, rx, ry, rz);
	m_CTActor->SetUserMatrix(matrix);
	ui->display->insertPlainText("the present transformation of CT:\n");
	getXYZRotationAngles(matrix);
	//exportCompositeModel();
	writeOBJCase();
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
	ui->display->insertPlainText("the present transformation of marker:\n");
	getXYZRotationAngles(matrix);
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
	std::ofstream outTXT("C:\\Users\\29477\\Desktop\\matrix.txt", ios::out);
	m_Markermatrix = m_MarkerActor->GetMatrix();
	m_CToriginmatrix = m_CTActor->GetMatrix();
	m_Toumomatrix = m_ToumoOriginActor->GetMatrix();
#pragma region test
	vtkMatrix4x4 *world2MarkerMatrix = m_Markermatrix;
	outTXT<<"The marker pistion and rotation:";
	world2MarkerMatrix->Print(outTXT);
	ui->display->insertPlainText("The marker pistion and rotation:\n");
	getXYZRotationAngles(world2MarkerMatrix);
	vtkMatrix4x4 *Marker2ToumoMatrix = vtkMatrix4x4::New();
	vtkMatrix4x4 *Marker2worldMatrix = vtkMatrix4x4::New();
	Marker2worldMatrix->Invert(world2MarkerMatrix, Marker2worldMatrix);
	vtkMatrix4x4 *Marker2CTMatrix = vtkMatrix4x4::New();
	Marker2CTMatrix->Multiply4x4(m_CToriginmatrix, Marker2worldMatrix, Marker2CTMatrix);
	outTXT<< "The CT in marker coordinates:";
	Marker2CTMatrix->Print(outTXT);
	ui->display->insertPlainText("the transformation of CT in marker coordinates:\n");
	getXYZRotationAngles(Marker2CTMatrix);
	Marker2ToumoMatrix->Multiply4x4(m_Toumomatrix, Marker2worldMatrix, Marker2ToumoMatrix);
	outTXT<<"The Toumo in marker coordinates:";
	Marker2ToumoMatrix->Print(outTXT); 
	ui->display->insertPlainText("the transformation of Toumo in marker coordinates:\n");
	getXYZRotationAngles(Marker2ToumoMatrix);
	vtkMatrix4x4 *CT2MarkerMatrix = vtkMatrix4x4::New();
	CT2MarkerMatrix->Invert(Marker2CTMatrix, CT2MarkerMatrix);
	vtkMatrix4x4 *CT2ToumoMatrix = vtkMatrix4x4::New();
	CT2ToumoMatrix->Multiply4x4(Marker2ToumoMatrix, CT2MarkerMatrix, CT2ToumoMatrix);
	outTXT<<"the CT2Toumo matrix is:";
	CT2ToumoMatrix->Print(outTXT);
	ui->display->insertPlainText("the CT2Toumo matrix is:\n");
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
	std::ofstream outTXT("C:\\Users\\29477\\Desktop\\matrix.txt", ios::out);
	if (firstClick2Marker == 0)
	{
		m_CToriginmatrix = m_CTActor->GetMatrix();
		printf("%lf,%lf,%lf", m_CTcenter[0], m_CTcenter[1], m_CTcenter[2]);
		vtkMatrix4x4* marker2CToriginMatrix = vtkMatrix4x4::New();
		marker2CToriginMatrix = getmarker2CToriginMatrix();
		//CTtrans->Multiply4x4(marker2CToriginMatrix, m_CToriginmatrix, CTtrans);
		//下一个语句，用于将CT模型的中心点恢复到坐标系原点
		//CTtrans->Multiply4x4(CTtrans, m_CTcenter2origin, CTtrans);
		ui->display->insertPlainText("the angles of CT in marker coordiates:\n");
		vtkMatrix4x4 *CTtrans = vtkMatrix4x4::New();
		CTtrans->Multiply4x4(marker2CToriginMatrix, m_CToriginmatrix, CTtrans);
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
	stringstream ss;
	ss << x1;
	ss << ',';
	ss << y1;
	ss << ',';
	ss << z1;
	string t;
	ss >> t;
	ui->display->insertPlainText("the rotation angles: ");
	ui->display->insertPlainText(QString::fromStdString(t));
	ui->display->insertPlainText("\n");
}

vtkMatrix4x4* registrationWidget::getmarker2CToriginMatrix() {
	std::ofstream outTXT("C:\\Users\\29477\\Desktop\\matrix.txt", ios::out);
	m_Markermatrix = m_MarkerActor->GetMatrix();
	m_CToriginmatrix = m_CTActor->GetMatrix();
	vtkMatrix4x4 *inverseMarkerMatrix = vtkMatrix4x4::New();//标记板坐标的逆矩阵
	inverseMarkerMatrix->Invert(m_Markermatrix, inverseMarkerMatrix);
	vtkMatrix4x4 *CTmarker2originMatrix = vtkMatrix4x4::New();
	CTmarker2originMatrix->Multiply4x4(m_CToriginmatrix, inverseMarkerMatrix, CTmarker2originMatrix);
	outTXT << "marker to CT origin matrix:";
	CTmarker2originMatrix->Print(outTXT);
	outTXT.close();
	return CTmarker2originMatrix;
}
void registrationWidget::constructCompositeModel() {
	vtkPoints *points = vtkPoints::New();
	points->InsertNextPoint(0.0, 0.0, 0.0);
	points->InsertNextPoint(3.0, 0.0, 0.0);
	points->InsertNextPoint(3.0, 0.0, -10.0);
	points->InsertNextPoint(1.0, 2.0, 0.0);
	points->InsertNextPoint(2.0, 2.0, 0.0);	
	vtkPolygon *polygon1 = vtkPolygon::New();
	polygon1->GetPointIds()->SetNumberOfIds(4);
	polygon1->GetPointIds()->SetId(0, 0);
	polygon1->GetPointIds()->SetId(1, 1);
	polygon1->GetPointIds()->SetId(2, 4);
	polygon1->GetPointIds()->SetId(3, 3);
	vtkPolygon *polygon2 = vtkPolygon::New();
	polygon2->GetPointIds()->SetNumberOfIds(3);
	polygon2->GetPointIds()->SetId(0, 1);
	polygon2->GetPointIds()->SetId(1, 2);
	polygon2->GetPointIds()->SetId(2, 4);
	vtkPolygon *polygon3 = vtkPolygon::New();
	polygon3->GetPointIds()->SetNumberOfIds(3);
	polygon3->GetPointIds()->SetId(0, 0);
	polygon3->GetPointIds()->SetId(1, 3);
	polygon3->GetPointIds()->SetId(2, 2);
	vtkPolygon *polygon4 = vtkPolygon::New();
	polygon4->GetPointIds()->SetNumberOfIds(3);
	polygon4->GetPointIds()->SetId(0, 2);
	polygon4->GetPointIds()->SetId(1, 3);
	polygon4->GetPointIds()->SetId(2, 4);
	vtkPolygon *polygon5 = vtkPolygon::New();
	polygon5->GetPointIds()->SetNumberOfIds(3);
	polygon5->GetPointIds()->SetId(0, 0);
	polygon5->GetPointIds()->SetId(1, 2);
	polygon5->GetPointIds()->SetId(2, 1);

	vtkCellArray *cell = vtkCellArray::New();
	cell->InsertNextCell(polygon1);
	cell->InsertNextCell(polygon2);
	cell->InsertNextCell(polygon3);
	cell->InsertNextCell(polygon4);
	cell->InsertNextCell(polygon5);
	compositepolydata->SetPoints(points);
	compositepolydata->SetPolys(cell);	
}
void registrationWidget::exportCompositeModel() {
	//用于最开始生成自定义组合模型
	vtkPolyDataMapper *compositeMapper = vtkPolyDataMapper::New();
	constructCompositeModel();
	compositeMapper->SetInputData(compositepolydata);
	vtkActor *compositeActor = vtkActor::New();
	compositeActor->SetMapper(compositeMapper);
	vtkRenderer *exportrenderer = vtkRenderer::New();
	exportrenderer->AddActor(compositeActor);
	m_exportWindow->AddRenderer(exportrenderer);
	writeOBJCase();
}
void  registrationWidget::prepareExportCTModel(){
	vtkRenderer *exportrenderer = vtkRenderer::New();
	exportrenderer->AddActor(m_CTActor);
	m_exportWindow->AddRenderer(exportrenderer);
}
int main(int argc, char *argv[]) {
	QApplication app(argc, argv);
	registrationWidget registrationWidget;
	registrationWidget.show();
	return app.exec();
}