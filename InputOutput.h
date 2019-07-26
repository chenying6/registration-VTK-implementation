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
	//����������ļ������ͣ������ļ��Ķ�ȡ
	vtkMapper * readCase(std::string fileName);
	//��OBJ��ʽ����ģ��
	void writeOBJCase(vtkRenderWindow* exportWin);
	//��STL��ʽ����ģ��
	void writeSTLCase(vtkMarchingCubes *data);
	//�Զ���������ģ��
	vtkPolyData * constructCompositeModel();
	//�����ģ�����
	void exportCompositeModel(vtkPolyData *compositepolydata, vtkRenderWindow *exportWin);
	//Ϊ���CTģ����׼�������ʼ��ʼ����������ÿһ�η�������ʱ��׼�����ᱨ��
	void prepareExportModel(vtkActor* actor, vtkRenderWindow *exportWin);
	//��TXT�ļ��ж�ȡ
	void readFromTXT(std::string fileName, float*& matrixArray);
	vtkAlgorithmOutput* ExtractToumoConnectivity(vtkAlgorithmOutput* input);
	template<class T>
	vtkPolyData* PolyDataReader(vtkSmartPointer<T> type, const char* fileName);
	template<class T>
	vtkImageData* ImageDataReader(vtkSmartPointer<T> type,const char* fileName);
	vtkAlgorithmOutput* DICOMReader(const char* fileName);
};
#endif // !INPUTOUTPUT_H
