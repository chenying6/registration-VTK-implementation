#ifndef STENCIL_H
#define STENCIL_H

#include <vtkImageViewer.h>

class transformation
{

public:
	transformation();
	~transformation();
	void StartTransformation();
	vtkImageViewer * viewer;
};


#endif