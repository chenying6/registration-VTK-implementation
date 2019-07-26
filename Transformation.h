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
	//以右手坐标系的方式设置位姿
	vtkMatrix4x4 * setTransformation_right(std::vector<int> array);
	//以左手坐标系的方式设置位姿，使其在Unity中的表示与期望的一样
	vtkMatrix4x4 * setTransformation_left(std::vector<int> array);
	//VTK中对于静系，默认的变换顺序根据先平移，再Z,X,Y的旋转顺序，计算绕各轴的旋转角度
	//注意，使用哪一个函数，应与setTransformation_right对应
	std::string getZXYRotationAngles(vtkMatrix4x4 *m);
	//Unity中对于静系，默认的变换顺序，根据先平移，再Y,X,Z的旋转顺序，计算绕各轴的旋转角度
	std::string getYXZRotationAngles(vtkMatrix4x4 *m);
	vtkMatrix4x4* setCurrentMatrix(float* matrixArray);
	vtkMatrix4x4* setCurrentMatrix(double matrix[3][4]);
	vtkMatrix4x4* setCurrentMatrix(float* r1, float* r2, float *r3);
};
#endif // !TRANSFORMATION_H

