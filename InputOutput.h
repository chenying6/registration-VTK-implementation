#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H
/*******************************************************************************
*Author: Chen Ying
*Content: It has the following functions:
		1. import different types of medical images;
		2. export STL OBJ models
		3. customly construct composite models with lines, vertices
		4. extract single model using the connectivity
*Date: 2019-4-17
********************************************************************************/

#include "vtkMarchingCubes.h"
#include "vtkMapper.h"
#include "vtkRenderWindow.h"
class InputOutput {
public:
	//根据输入的文件的类型，进行文件的读取
	vtkMapper * readCase(std::string fileName);
	//以OBJ格式导出模型
	void writeOBJCase(vtkRenderWindow* exportWin);
	//以STL格式导出模型
	void writeSTLCase(vtkMarchingCubes *data);
	//自定义绘制组合模型
	vtkPolyData * constructCompositeModel();
	//将组合模型输出
	void exportCompositeModel(vtkPolyData *compositepolydata, vtkRenderWindow *exportWin);
	//为输出CT模型做准备，在最开始初始化，而不是每一次发生调用时做准备（会报错）
	void prepareExportModel(vtkActor* actor, vtkRenderWindow *exportWin);
	//从TXT文件中读取
	void readFromTXT(std::string fileName, float*& matrixArray);
	vtkAlgorithmOutput* ExtractToumoConnectivity(vtkAlgorithmOutput* input);
	template<class T>
	vtkPolyData* PolyDataReader(vtkSmartPointer<T> type, const char* fileName);
	template<class T>
	vtkImageData* ImageDataReader(vtkSmartPointer<T> type,const char* fileName);
	vtkAlgorithmOutput* DICOMReader(const char* fileName);
};
#endif // !INPUTOUTPUT_H
