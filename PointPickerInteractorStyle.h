#ifndef POINTPICKERINTERACTORSTYLE_H
#define POINTPICKERINTERACTORSTYLE_H
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
