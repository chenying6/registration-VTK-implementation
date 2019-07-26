#ifndef POINTPICKERINTERACTORSTYLE_H
#define POINTPICKERINTERACTORSTYLE_H
/*******************************************************************************
*Author: Chen Ying
*Content: It has the following functions:
		1. pick points and show their coordinates
*Date: 2019-4-17
********************************************************************************/
#include "vtkInteractorStyleTrackballCamera.h"
class PointPickerInteractorStyle :public vtkInteractorStyleTrackballCamera {
public:
	static PointPickerInteractorStyle *New();
	vtkTypeMacro(PointPickerInteractorStyle, vtkInteractorStyleTrackballCamera);
	PointPickerInteractorStyle();
	~PointPickerInteractorStyle();

	virtual void OnRightButtonDown();
};
#endif // !POINTPICKERINTERACTORSTYLE_H
