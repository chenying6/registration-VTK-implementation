#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H
#include "vtkMatrix4x4.h"
#include <vector>
class Transformation {
public:
	//����������ϵ�ķ�ʽ����λ��
	vtkMatrix4x4 * setTransformation_right(std::vector<int> array);
	//����������ϵ�ķ�ʽ����λ�ˣ�ʹ����Unity�еı�ʾ��������һ��
	vtkMatrix4x4 * setTransformation_left(std::vector<int> array);
	//VTK�ж��ھ�ϵ��Ĭ�ϵı任˳�������ƽ�ƣ���Z,X,Y����ת˳�򣬼����Ƹ������ת�Ƕȣ�ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	std::string getZXYRotationAngles(vtkMatrix4x4 *m);
	//Unity�ж��ھ�ϵ��Ĭ�ϵı任˳�򣬸�����ƽ�ƣ���Y,X,Z����ת˳�򣬼����Ƹ������ת�Ƕȣ�ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	std::string getYXZRotationAngles(vtkMatrix4x4 *m);
	vtkMatrix4x4* setCurrentMatrix(float* matrixArray);
	vtkMatrix4x4* setCurrentMatrix(double matrix[3][4]);
	vtkMatrix4x4* setCurrentMatrix(float* r1, float* r2, float *r3);
};
#endif // !TRANSFORMATION_H

