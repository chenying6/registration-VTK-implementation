#include "Transformation.h"
#include "vtkTransform.h"
#include <sstream>
#include <qstring.h>
#define PI 3.14159265
vtkMatrix4x4 * Transformation::setTransformation_right(std::vector<int> array)
{
	vtkTransform *transformation = vtkTransform::New();
	vtkMatrix4x4 *transMatrix = vtkMatrix4x4::New();
	transformation->Identity();
	transformation->Translate(array[0], array[1], array[2]);
	transformation->RotateY(array[4]);
	transformation->RotateX(array[3]);
	transformation->RotateZ(array[5]);
	transMatrix = transformation->GetMatrix();
	return transMatrix;
}
vtkMatrix4x4 * Transformation::setTransformation_left(std::vector<int> array)
{
	//(-x, y, z, rx, -ry, -rz);
	std::vector<int> rightVector;
	rightVector.push_back(-array[0]);
	rightVector.push_back(array[1]);
	rightVector.push_back(array[2]);
	rightVector.push_back(array[3]);
	rightVector.push_back(-array[4]);
	rightVector.push_back(-array[5]);
	vtkMatrix4x4 *matrix = setTransformation_right(rightVector);
	return matrix;
}
vtkMatrix4x4 * Transformation::setCurrentMatrix(double m[3][4]) {
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix->Identity();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			matrix->SetElement(i, j, m[i][j]);
		}
	}
	return matrix;
}
vtkMatrix4x4* Transformation::setCurrentMatrix(float* matrixArray){
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix->Identity();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			matrix->SetElement(i, j, matrixArray[4*i+j]);
		}
	}
	return matrix;
}
vtkMatrix4x4* Transformation::setCurrentMatrix(float* r1, float* r2, float *r3) {
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix->SetElement(0, 0, r1[0]);
	matrix->SetElement(0, 1, r1[1]);
	matrix->SetElement(0, 2, r1[2]);
	matrix->SetElement(0, 3, r1[3]);
	matrix->SetElement(1, 0, r2[0]);
	matrix->SetElement(1, 1, r2[1]);
	matrix->SetElement(1, 2, r2[2]);
	matrix->SetElement(1, 3, r2[3]);
	matrix->SetElement(2, 0, r3[0]);
	matrix->SetElement(2, 1, r3[1]);
	matrix->SetElement(2, 2, r3[2]);
	matrix->SetElement(2, 3, r3[3]);
	return matrix;
}
std::string Transformation::getZXYRotationAngles(vtkMatrix4x4 *matrix)
{
	double y1 = 0, x1 = 0, z1 = 0;
	y1 = atan2(-matrix->GetElement(2, 0), matrix->GetElement(2, 2))*180.0 / PI;
	x1 = atan2(matrix->GetElement(2, 1), sqrt(matrix->GetElement(2, 0)*matrix->GetElement(2, 0) + matrix->GetElement(2, 2)*matrix->GetElement(2, 2)))*180.0 / PI;
	z1 = atan2(-matrix->GetElement(0, 1), matrix->GetElement(1, 1))*180.0 / PI;
	std::stringstream ss;
	ss << x1<< ','<< y1<< ','<< z1;
	std::string t;
	ss >> t;
	return t;	
}
std::string Transformation::getYXZRotationAngles(vtkMatrix4x4 *matrix) {
	double y1 = 0, x1 = 0, z1 = 0;
	x1 = atan2(-matrix->GetElement(1, 2), sqrt(matrix->GetElement(0, 2)*matrix->GetElement(0, 2) + matrix->GetElement(2, 2)*matrix->GetElement(2, 2)))*180.0 / PI;
	y1 = atan2(matrix->GetElement(0, 2), matrix->GetElement(2, 2)) * 180.0 / PI;
	z1 = atan2(matrix->GetElement(1, 0), matrix->GetElement(1, 1)) * 180.0 / PI;
	std::stringstream ss;
	ss << x1<< ','<< y1<< ','<< z1;
	std::string t;
	ss >> t;
	return t;
}