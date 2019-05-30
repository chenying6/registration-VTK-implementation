/*
*Author: Chen Ying
*Version:1.4
*Date: 2019-4-22
*History change description:
add an axisActor class to generate axis model for Unity use;
*/

#ifndef AxisActor_h
#define AxisActor_h

#include "vtkRenderingAnnotationModule.h" // For export macro
#include "vtkProp3D.h"

class vtkActor;
class vtkCaptionActor2D;
class vtkConeSource;
class vtkCylinderSource;
class vtkLineSource;
class vtkPolyData;
class vtkPropCollection;
class vtkProperty;
class vtkRenderer;
class vtkSphereSource;

class AxisActor : public vtkProp3D
{
public:
	static AxisActor *New();
	void PrintSelf(ostream& os, vtkIndent indent) override;

	/**
	* For some exporters and other other operations we must be
	* able to collect all the actors or volumes. These methods
	* are used in that process.
	*/
	void GetActors(vtkPropCollection *) override;

	//@{
	/**
	* Support the standard render methods.
	*/
	int RenderOpaqueGeometry(vtkViewport *viewport) override;
	int RenderTranslucentPolygonalGeometry(vtkViewport *viewport) override;
	int RenderOverlay(vtkViewport *viewport) override;
	//@}

	/**
	* Does this prop have some translucent polygonal geometry?
	*/
	int HasTranslucentPolygonalGeometry() override;

	///**
	//* Shallow copy of an axes actor. Overloads the virtual vtkProp method.
	//*/
	//void ShallowCopy(vtkProp *prop) override;

	/**
	* Release any graphics resources that are being consumed by this actor.
	* The parameter window could be used to determine which graphic
	* resources to release.
	*/
	void ReleaseGraphicsResources(vtkWindow *) override;

	//@{
	/**
	* Get the bounds for this Actor as (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax). (The
	* method GetBounds(double bounds[6]) is available from the superclass.)
	*/
	void GetBounds(double bounds[6]);
	double *GetBounds() override;
	//@}

	/**
	* Get the actors mtime plus consider its properties and texture if set.
	*/
	vtkMTimeType GetMTime() override;

	/**
	* Return the mtime of anything that would cause the rendered image to
	* appear differently. Usually this involves checking the mtime of the
	* prop plus anything else it depends on such as properties, textures
	* etc.
	*/
	vtkMTimeType GetRedrawMTime() override;

	//@{
	/**
	* Set the total length of the axes in 3 dimensions.
	*/
	void SetTotalLength(double v[3])
	{
		this->SetTotalLength(v[0], v[1], v[2]);
	}
	void SetTotalLength(double x, double y, double z);
	vtkGetVectorMacro(TotalLength, double, 3);
	//@}

	//@{
	/**
	* Set the normalized (0-1) length of the shaft.
	*/
	void SetNormalizedShaftLength(double v[3])
	{
		this->SetNormalizedShaftLength(v[0], v[1], v[2]);
	}
	void SetNormalizedShaftLength(double x, double y, double z);
	vtkGetVectorMacro(NormalizedShaftLength, double, 3);
	//@}

	//@{
	/**
	* Set the normalized (0-1) length of the tip.  Normally, this would be
	* 1 - the normalized length of the shaft.
	*/
	void SetNormalizedTipLength(double v[3])
	{
		this->SetNormalizedTipLength(v[0], v[1], v[2]);
	}
	void SetNormalizedTipLength(double x, double y, double z);
	vtkGetVectorMacro(NormalizedTipLength, double, 3);
	//@}

	//@{
	/**
	* Set the normalized (0-1) position of the label along the length of
	* the shaft.  A value > 1 is permissible.
	*/
	void SetNormalizedLabelPosition(double v[3])
	{
		this->SetNormalizedLabelPosition(v[0], v[1], v[2]);
	}
	void SetNormalizedLabelPosition(double x, double y, double z);
	vtkGetVectorMacro(NormalizedLabelPosition, double, 3);
	//@}

	//@{
	/**
	* Set/get the resolution of the pieces of the axes actor.
	*/
	vtkSetClampMacro(ConeResolution, int, 3, 128);
	vtkGetMacro(ConeResolution, int);
	vtkSetClampMacro(SphereResolution, int, 3, 128);
	vtkGetMacro(SphereResolution, int);
	vtkSetClampMacro(CylinderResolution, int, 3, 128);
	vtkGetMacro(CylinderResolution, int);
	//@}

	//@{
	/**
	* Set/get the radius of the pieces of the axes actor.
	*/
	vtkSetClampMacro(ConeRadius, double, 0, VTK_FLOAT_MAX);
	vtkGetMacro(ConeRadius, double);
	vtkSetClampMacro(SphereRadius, double, 0, VTK_FLOAT_MAX);
	vtkGetMacro(SphereRadius, double);
	vtkSetClampMacro(CylinderRadius, double, 0, VTK_FLOAT_MAX);
	vtkGetMacro(CylinderRadius, double);
	//@}

	//@{
	/**
	* Set the type of the shaft to a cylinder, line, or user defined geometry.
	*/
	void SetShaftType(int type);
	void SetShaftTypeToCylinder()
	{
		this->SetShaftType(AxisActor::CYLINDER_SHAFT);
	}
	void SetShaftTypeToLine()
	{
		this->SetShaftType(AxisActor::LINE_SHAFT);
	}
	void SetShaftTypeToUserDefined()
	{
		this->SetShaftType(AxisActor::USER_DEFINED_SHAFT);
	}
	vtkGetMacro(ShaftType, int);
	//@}

	//@{
	/**
	* Set the type of the tip to a cone, sphere, or user defined geometry.
	*/
	void SetTipType(int type);
	void SetTipTypeToCone()
	{
		this->SetTipType(AxisActor::CONE_TIP);
	}
	void SetTipTypeToSphere()
	{
		this->SetTipType(AxisActor::SPHERE_TIP);
	}
	void SetTipTypeToUserDefined()
	{
		this->SetTipType(AxisActor::USER_DEFINED_TIP);
	}
	vtkGetMacro(TipType, int);
	//@}

	//@{
	/**
	* Set the user defined tip polydata.
	*/
	void SetUserDefinedTip(vtkPolyData *);
	vtkGetObjectMacro(UserDefinedTip, vtkPolyData);
	//@}

	//@{
	/**
	* Set the user defined shaft polydata.
	*/
	void SetUserDefinedShaft(vtkPolyData *);
	vtkGetObjectMacro(UserDefinedShaft, vtkPolyData);
	//@}

	//@{
	/**
	* Get the tip properties.
	*/
	vtkProperty *GetXAxisTipProperty();
	vtkProperty *GetYAxisTipProperty();
	vtkProperty *GetZAxisTipProperty();
	//@}

	//@{
	/**
	* Get the shaft properties.
	*/
	vtkProperty *GetXAxisShaftProperty();
	vtkProperty *GetYAxisShaftProperty();
	vtkProperty *GetZAxisShaftProperty();
	//@}

	/**
	* Retrieve handles to the X, Y and Z axis (so that you can set their text
	* properties for example)
	*/
	vtkCaptionActor2D *GetXAxisCaptionActor2D()
	{
		return this->XAxisLabel;
	}
	vtkCaptionActor2D *GetYAxisCaptionActor2D()
	{
		return this->YAxisLabel;
	}
	vtkCaptionActor2D *GetZAxisCaptionActor2D()
	{
		return this->ZAxisLabel;
	}

	//@{
	/**
	* Set/get the label text.
	*/
	vtkSetStringMacro(XAxisLabelText);
	vtkGetStringMacro(XAxisLabelText);
	vtkSetStringMacro(YAxisLabelText);
	vtkGetStringMacro(YAxisLabelText);
	vtkSetStringMacro(ZAxisLabelText);
	vtkGetStringMacro(ZAxisLabelText);
	//@}

	//@{
	/**
	* Enable/disable drawing the axis labels.
	*/
	vtkSetMacro(AxisLabels, int);
	vtkGetMacro(AxisLabels, int);
	vtkBooleanMacro(AxisLabels, int);
	//@}

	enum
	{
		CYLINDER_SHAFT,
		LINE_SHAFT,
		USER_DEFINED_SHAFT
	};

	enum
	{
		CONE_TIP,
		SPHERE_TIP,
		USER_DEFINED_TIP
	};

protected:
	AxisActor();
	~AxisActor() override;

	vtkCylinderSource *CylinderSource;
	vtkLineSource     *LineSource;
	vtkConeSource     *ConeSource;
	vtkSphereSource   *SphereSource;

	vtkActor          *XAxisShaft;
	vtkActor          *YAxisShaft;
	vtkActor          *ZAxisShaft;

	vtkActor          *XAxisTip;
	vtkActor          *YAxisTip;
	vtkActor          *ZAxisTip;

	void               UpdateProps();

	double             TotalLength[3];
	double             NormalizedShaftLength[3];
	double             NormalizedTipLength[3];
	double             NormalizedLabelPosition[3];

	int                ShaftType;
	int                TipType;

	vtkPolyData       *UserDefinedTip;
	vtkPolyData       *UserDefinedShaft;

	char              *XAxisLabelText;
	char              *YAxisLabelText;
	char              *ZAxisLabelText;

	vtkCaptionActor2D *XAxisLabel;
	vtkCaptionActor2D *YAxisLabel;
	vtkCaptionActor2D *ZAxisLabel;

	int                AxisLabels;


	int                ConeResolution;
	int                SphereResolution;
	int                CylinderResolution;

	double             ConeRadius;
	double             SphereRadius;
	double             CylinderRadius;

private:
	AxisActor(const AxisActor&) = delete;
	void operator=(const AxisActor&) = delete;
};

#endif

