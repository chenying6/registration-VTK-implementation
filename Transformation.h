#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H
/*******************************************************************************
*Author: Chen Ying
*Content: It has the following functions:
		1. set the transformation according to left or right coodinate routine
		2. get the Euler rotation angles according to the rotation routine of Unity
*Date: 2019-4-17
********************************************************************************/
#include "vtkMatrix4x4.h"
#include <vector>
class Transformation {
public:
	//����������ϵ�ķ�ʽ����λ��
	vtkMatrix4x4 * setTransformation_right(std::vector<int> array);
	//����������ϵ�ķ�ʽ����λ�ˣ�ʹ����Unity�еı�ʾ��������һ��
	vtkMatrix4x4 * setTransformation_left(std::vector<int> array);
	//VTK�ж��ھ�ϵ��Ĭ�ϵı任˳�������ƽ�ƣ���Z,X,Y����ת˳�򣬼����Ƹ������ת�Ƕ�
	//ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	std::string getZXYRotationAngles(vtkMatrix4x4 *m);
	//Unity�ж��ھ�ϵ��Ĭ�ϵı任˳�򣬸�����ƽ�ƣ���Y,X,Z����ת˳�򣬼����Ƹ������ת�Ƕ�
	std::string getYXZRotationAngles(vtkMatrix4x4 *m);
	vtkMatrix4x4* setCurrentMatrix(float* matrixArray);
	vtkMatrix4x4* setCurrentMatrix(double matrix[3][4]);
	vtkMatrix4x4* setCurrentMatrix(float* r1, float* r2, float *r3);
};
#endif // !TRANSFORMATION_H

