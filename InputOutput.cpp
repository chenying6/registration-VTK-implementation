#include "InputOutput.h"
#include "vtksys/SystemTools.hxx"
#include "vtkDICOMImageReader.h"
#include "vtkOBJExporter.h"
#include "vtkDecimatePro.h"
#include "vtkOBJReader.h"
#include "vtkSTLReader.h"
#include "vtkNIFTIImageReader.h"
#include "vtkImageShrink3D.h"
#include "vtkSTLWriter.h"
#include "vtkPolyDataMapper.h"
#include "vtkDataSetMapper.h"
#include "vtkImageData.h"
#include "vtkPolygon.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkAlgorithmOutput.h"
vtkMapper * InputOutput::readCase(std::string fileName)
{
	std::string extension = vtksys::SystemTools::GetFilenameLastExtension(fileName);
	vtkPolyData *origindata = vtkPolyData::New();
	vtkPolyData *CTdata = vtkPolyData::New();
	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	vtkPolyDataMapper* polyMapper = vtkPolyDataMapper::New();
	vtkDataSetMapper* imageMapper = vtkDataSetMapper::New();
	if (extension == ".obj")
	{
		vtkSmartPointer<vtkOBJReader> obj = vtkSmartPointer<vtkOBJReader>::New();
		CTdata = this->PolyDataReader(obj, fileName.data());
		polyMapper->SetInputData(CTdata);
		return polyMapper;
	}
	else if (extension == ".stl")
	{
		vtkSmartPointer<vtkSTLReader> stl = vtkSmartPointer<vtkSTLReader>::New();
		CTdata = this->PolyDataReader(stl, fileName.data());
		polyMapper->SetInputData(CTdata);
		return polyMapper;
	}
	else if (extension == ".nii")
	{
		vtkSmartPointer<vtkNIFTIImageReader> nifti = vtkSmartPointer<vtkNIFTIImageReader>::New();
		imageData = this->ImageDataReader(nifti,fileName.data());
		imageMapper->SetInputData(imageData);
		return imageMapper;
	}
	else
	{	
		vtkSmartPointer<vtkImageShrink3D> shrink = vtkSmartPointer<vtkImageShrink3D>::New();
		vtkSmartPointer<vtkMarchingCubes> CTdata = vtkSmartPointer<vtkMarchingCubes>::New();
		vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
		shrink->SetInputConnection(DICOMReader(fileName.data()));
		shrink->SetShrinkFactors(1, 1, 1);
		shrink->AveragingOn();
		shrink->Update();
		CTdata->SetInputConnection(shrink->GetOutputPort());
		CTdata->ComputeNormalsOn();
		CTdata->SetValue(0, -150);
		CTdata->Update();
		decimate->SetInputConnection(CTdata->GetOutputPort());
		decimate->SetTargetReduction(0.8);
		decimate->Update();
		vtkAlgorithmOutput* toumoConnectivity = vtkAlgorithmOutput::New();
		toumoConnectivity = ExtractToumoConnectivity(decimate->GetOutputPort());
		//polyMapper->SetInputConnection();
		polyMapper->SetInputConnection(toumoConnectivity);
		polyMapper->Update();
		return polyMapper;
	}
}

vtkAlgorithmOutput* InputOutput::ExtractToumoConnectivity(vtkAlgorithmOutput* input) {
	vtkPolyDataConnectivityFilter* connectivityFilter = vtkPolyDataConnectivityFilter::New();
	connectivityFilter->SetInputConnection(input);
	connectivityFilter->SetExtractionModeToLargestRegion();
	connectivityFilter->Update();
	return connectivityFilter->GetOutputPort();
}
template<class T>
vtkPolyData* InputOutput::PolyDataReader(vtkSmartPointer<T> type, const char* fileName) {
	type->SetFileName(fileName);
	type->Update();
	return type->GetOutput();
}
template<class T>
vtkImageData* InputOutput::ImageDataReader(vtkSmartPointer<T>  type, const char* fileName) {
	type->SetFileName(fileName);
	type->Update();
	return type->GetOutput();
}
vtkAlgorithmOutput* InputOutput::DICOMReader(const char* fileName) {
	vtkDICOMImageReader *dicomReader = vtkDICOMImageReader::New();
	dicomReader->SetDataByteOrderToLittleEndian();
	dicomReader->SetDirectoryName(fileName);
	dicomReader->ReleaseDataFlagOn();
	dicomReader->SetDataScalarTypeToUnsignedShort();
	dicomReader->Update();
	dicomReader->GlobalWarningDisplayOff();
	return dicomReader->GetOutputPort();
}
vtkPolyData * InputOutput::constructCompositeModel() {
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(0.0, 0.0, 0.0);
	points->InsertNextPoint(3.0, 0.0, 0.0);
	points->InsertNextPoint(3.0, 0.0, -10.0);
	points->InsertNextPoint(1.0, 2.0, 0.0);
	points->InsertNextPoint(2.0, 2.0, 0.0);
	vtkSmartPointer<vtkPolygon> polygon1 = vtkSmartPointer<vtkPolygon>::New();
	polygon1->GetPointIds()->SetNumberOfIds(4);
	polygon1->GetPointIds()->SetId(0, 0);
	polygon1->GetPointIds()->SetId(1, 1);
	polygon1->GetPointIds()->SetId(2, 4);
	polygon1->GetPointIds()->SetId(3, 3);
	vtkSmartPointer<vtkPolygon> polygon2 = vtkSmartPointer<vtkPolygon>::New();
	polygon2->GetPointIds()->SetNumberOfIds(3);
	polygon2->GetPointIds()->SetId(0, 1);
	polygon2->GetPointIds()->SetId(1, 2);
	polygon2->GetPointIds()->SetId(2, 4);
	vtkSmartPointer<vtkPolygon> polygon3 = vtkSmartPointer<vtkPolygon>::New();
	polygon3->GetPointIds()->SetNumberOfIds(3);
	polygon3->GetPointIds()->SetId(0, 0);
	polygon3->GetPointIds()->SetId(1, 3);
	polygon3->GetPointIds()->SetId(2, 2);
	vtkSmartPointer<vtkPolygon> polygon4 = vtkSmartPointer<vtkPolygon>::New();
	polygon4->GetPointIds()->SetNumberOfIds(3);
	polygon4->GetPointIds()->SetId(0, 2);
	polygon4->GetPointIds()->SetId(1, 3);
	polygon4->GetPointIds()->SetId(2, 4);
	vtkSmartPointer<vtkPolygon> polygon5 = vtkSmartPointer<vtkPolygon>::New();
	polygon5->GetPointIds()->SetNumberOfIds(3);
	polygon5->GetPointIds()->SetId(0, 0);
	polygon5->GetPointIds()->SetId(1, 2);
	polygon5->GetPointIds()->SetId(2, 1);

	vtkSmartPointer<vtkCellArray> cell = vtkSmartPointer<vtkCellArray>::New();
	cell->InsertNextCell(polygon1);
	cell->InsertNextCell(polygon2);
	cell->InsertNextCell(polygon3);
	cell->InsertNextCell(polygon4);
	cell->InsertNextCell(polygon5);

	vtkPolyData *compositepolydata = vtkPolyData::New();
	compositepolydata->SetPoints(points);
	compositepolydata->SetPolys(cell);
	return compositepolydata;
}
void InputOutput::exportCompositeModel(vtkPolyData *compositepolydata, vtkRenderWindow *exportWin) {
	//用于最开始生成自定义组合模型
	vtkPolyDataMapper *compositeMapper = vtkPolyDataMapper::New();
	constructCompositeModel();
	compositeMapper->SetInputData(compositepolydata);
	vtkActor *compositeActor = vtkActor::New();
	compositeActor->SetMapper(compositeMapper);
	vtkRenderer *exportrenderer = vtkRenderer::New();
	exportrenderer->AddActor(compositeActor);
	exportWin->AddRenderer(exportrenderer);
	writeOBJCase(exportWin);
}
void  InputOutput::prepareExportModel(vtkActor* actor, vtkRenderWindow *exportWin) {
	vtkRenderer *exportrenderer = vtkRenderer::New();
	exportrenderer->AddActor(actor);
	exportWin->AddRenderer(exportrenderer);
}
void InputOutput::writeOBJCase(vtkRenderWindow *exportWin)
{
	vtkSmartPointer<vtkOBJExporter> objExporter = vtkSmartPointer<vtkOBJExporter>::New();
	objExporter->SetFilePrefix("registrationTest");
	objExporter->SetInput(exportWin);
	objExporter->Write();
}
void InputOutput::writeSTLCase(vtkMarchingCubes *data)
{
	vtkSmartPointer<vtkSTLWriter> stlExporter = vtkSmartPointer<vtkSTLWriter>::New();
	//stlExporter->SetInputData(data);
	stlExporter->SetInputConnection(data->GetOutputPort());
	stlExporter->SetFileTypeToBinary();
	stlExporter->SetFileName("C:\\Users\\29477\\Desktop\\registrationTest.stl");
	stlExporter->Update();
	stlExporter->Write();
}