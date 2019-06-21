#include "PointPickerInteractorStyle.h"
#include "QVTKInteractor.h"
#include "vtkPointPicker.h"
#include "vtkRenderWindow.h"
#include "vtkRendererCollection.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
vtkStandardNewMacro(PointPickerInteractorStyle);
PointPickerInteractorStyle::PointPickerInteractorStyle() {

}
PointPickerInteractorStyle::~PointPickerInteractorStyle() {

}

void PointPickerInteractorStyle::OnRightButtonDown() {
		std::cout << "Picking pixel:" << this->Interactor->GetEventPosition()[0] << " "
			<< this->Interactor->GetEventPosition()[1] <<std::endl;

		this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1],
			0,
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		std::cout << "Picked value:" << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;
		vtkSphereSource* sphere = vtkSphereSource::New();
		sphere->Update();
		vtkPolyDataMapper* sphereMapper = vtkPolyDataMapper::New();
		sphereMapper->SetInputConnection(sphere->GetOutputPort());
		vtkActor* sphereActor = vtkActor::New();
		sphereActor->SetMapper(sphereMapper);
		sphereActor->SetPosition(picked);
		sphereActor->SetScale(10);
		sphereActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(sphereActor);

		vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	}
