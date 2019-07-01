#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H
#endif // !TRANSFORMATION_H
#include "vtkMatrix4x4.h"
class Transformation {
public:
	//����������ϵ�ķ�ʽ����λ��
	vtkMatrix4x4 * setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	//����������ϵ�ķ�ʽ����λ�ˣ�ʹ����Unity�еı�ʾ��������һ��
	vtkMatrix4x4 * setTransformation_left(const float x, const float y, const float z, const float rx, const float ry, const float rz);
	//VTK�ж��ھ�ϵ��Ĭ�ϵı任˳�������ƽ�ƣ���Z,X,Y����ת˳�򣬼����Ƹ������ת�Ƕȣ�ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	std::string getZXYRotationAngles(vtkMatrix4x4 *m);
	//Unity�ж��ھ�ϵ��Ĭ�ϵı任˳�򣬸�����ƽ�ƣ���Y,X,Z����ת˳�򣬼����Ƹ������ת�Ƕȣ�ע�⣬ʹ����һ��������Ӧ��setTransformation_right��Ӧ
	std::string getYXZRotationAngles(vtkMatrix4x4 *m);
	vtkMatrix4x4* setCurrentMatrix(double matrix[4][4]);
	vtkMatrix4x4* setCurrentMatrix(float* r1, float* r2, float *r3, float* r4);
	
};

