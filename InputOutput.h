#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H
using namespace std;

#include "vtkMarchingCubes.h"
#include "vtkMapper.h"
#include "vtkRenderWindow.h"
class InputOutput {
public:
	//����������ļ������ͣ������ļ��Ķ�ȡ
	vtkMapper * readCase(std::string fileName);
	//��obj��ʽ����ģ��
	void writeOBJCase(vtkRenderWindow *exportWin);
	//��STL��ʽ����ģ��
	void writeSTLCase(vtkMarchingCubes *data);
	//�Զ���������ģ��
	vtkPolyData * constructCompositeModel();
	//�����ģ�����
	void exportCompositeModel(vtkPolyData *compositepolydata, vtkRenderWindow *exportWin);
	//Ϊ���CTģ����׼�������ʼ��ʼ����������ÿһ�η�������ʱ��׼�����ᱨ��
	void prepareExportModel(vtkActor* actor, vtkRenderWindow *exportWin);
	vtkAlgorithmOutput* ExtractToumoConnectivity(vtkAlgorithmOutput* input);
	template<class T>
	vtkPolyData* PolyDataReader(vtkSmartPointer<T> type, const char* fileName);
	template<class T>
	vtkImageData* ImageDataReader(vtkSmartPointer<T> type,const char* fileName);
	vtkAlgorithmOutput* DICOMReader(const char* fileName);
};
#endif // !INPUTOUTPUT_H
