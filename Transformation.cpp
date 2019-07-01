#include "Transformation.h"
#include "vtkTransform.h"
#include <sstream>
#include <qstring.h>
#define PI 3.14159265
vtkMatrix4x4 * Transformation::setTransformation_right(const float x, const float y, const float z, const float rx, const float ry, const float rz)
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
vtkMatrix4x4 * Transformation::setTransformation_left(const float x, const float y, const float z, const float rx, const float ry, const float rz)
{
	vtkMatrix4x4 *matrix = setTransformation_right(-x, y, z, rx, -ry, -rz);
	return matrix;
}
vtkMatrix4x4 * Transformation::setCurrentMatrix(double m[4][4]) {
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
vtkMatrix4x4* Transformation::setCurrentMatrix(float* r1, float* r2, float *r3, float* r4) {
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
	matrix->SetElement(3, 0, r4[0]);
	matrix->SetElement(3, 1, r4[1]);
	matrix->SetElement(3, 2, r4[2]);
	matrix->SetElement(3, 3, r4[3]);
	return matrix;
}
std::string Transformation::getZXYRotationAngles(vtkMatrix4x4 *matrix)
{
	double y1 = 0, x1 = 0, z1 = 0;
	y1 = atan2(-matrix->GetElement(2, 0), matrix->GetElement(2, 2))*180.0 / PI;
	x1 = atan2(matrix->GetElement(2, 1), sqrt(matrix->GetElement(2, 0)*matrix->GetElement(2, 0) + matrix->GetElement(2, 2)*matrix->GetElement(2, 2)))*180.0 / PI;
	z1 = atan2(-matrix->GetElement(0, 1), matrix->GetElement(1, 1))*180.0 / PI;
	std::stringstream ss;
	ss << x1;
	ss << ',';
	ss << y1;
	ss << ',';
	ss << z1;
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
	ss << x1;
	ss << ',';
	ss << y1;
	ss << ',';
	ss << z1;
	std::string t;
	ss >> t;
	return t;
}
