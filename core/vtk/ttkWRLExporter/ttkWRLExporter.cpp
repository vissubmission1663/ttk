// VTK includes
#include <vtkActor.h>
#include <vtkActorCollection.h>
#include <vtkAssemblyPath.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkGeometryFilter.h>
#include <vtkLightCollection.h>
#include <vtkMapper.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVRMLExporter.h>

// base code includes
#include <Debug.h>

#include <ttkWRLExporter.h>

/// @cond

TTKWRLEXPORTER_EXPORT vtkPolyData *ttkWRLExporterPolyData_ = nullptr;

// Over-ride the appropriate functions of the vtkVRMLExporter class.
TTKWRLEXPORTER_EXPORT void vtkVRMLExporter::WriteAnActor(vtkActor *anActor,
                                                         FILE *fp) {

  ttk::Debug debugInfo;
  debugInfo.setDebugMsgPrefix("WRLExporter");
  debugInfo.printWrn("Using TTK fix for VRML export...");

  vtkSmartPointer<vtkPolyData> pd;
  vtkPointData *pntData;
  vtkPoints *points;
  vtkDataArray *normals = nullptr;
  vtkDataArray *tcoords = nullptr;
  int i, i1, i2;
  double *tempd;
  vtkCellArray *cells;
  vtkIdType npts = 0;
#ifdef VTK_CELL_ARRAY_V2
  vtkIdType const *indx = nullptr;
#else
  vtkIdType *indx = nullptr;
#endif
  int pointDataWritten = 0;
  vtkPolyDataMapper *pm;
  vtkUnsignedCharArray *colors;
  double *p;
  unsigned char *c;
  vtkTransform *trans;

  // see if the actor has a mapper. it could be an assembly
  if(anActor->GetMapper() == nullptr) {
    return;
  }

  if(anActor->GetVisibility() == 0) {
    return;
  }

  // first stuff out the transform
  trans = vtkTransform::New();
  trans->SetMatrix(anActor->vtkProp3D::GetMatrix());

  fprintf(fp, "    Transform {\n");
  tempd = trans->GetPosition();
  fprintf(fp, "      translation %g %g %g\n", tempd[0], tempd[1], tempd[2]);
  tempd = trans->GetOrientationWXYZ();
  fprintf(fp, "      rotation %g %g %g %g\n", tempd[1], tempd[2], tempd[3],
          tempd[0] * vtkMath::Pi() / 180.0);
  tempd = trans->GetScale();
  fprintf(fp, "      scale %g %g %g\n", tempd[0], tempd[1], tempd[2]);
  fprintf(fp, "      children [\n");
  trans->Delete();

  vtkDataObject *inputDO = anActor->GetMapper()->GetInputDataObject(0, 0);

  if(inputDO == nullptr) {
    return;
  }

  // we really want polydata
  if(inputDO->IsA("vtkCompositeDataSet")) {
    vtkCompositeDataGeometryFilter *gf = vtkCompositeDataGeometryFilter::New();
    gf->SetInputConnection(anActor->GetMapper()->GetInputConnection(0, 0));
    gf->Update();
    pd = gf->GetOutput();
    gf->Delete();
  } else if(inputDO->GetDataObjectType() != VTK_POLY_DATA) {
    vtkGeometryFilter *gf = vtkGeometryFilter::New();
    gf->SetInputConnection(anActor->GetMapper()->GetInputConnection(0, 0));
    gf->Update();
    pd = gf->GetOutput();
    gf->Delete();
  } else {
    anActor->GetMapper()->Update();
    pd = static_cast<vtkPolyData *>(inputDO);
  }

  // BUG fix here
  ttkWRLExporterPolyData_ = static_cast<vtkPolyData *>(pd);
  // end of BUG fix here

  pm = vtkPolyDataMapper::New();
  pm->SetInputData(pd);
  pm->SetScalarRange(anActor->GetMapper()->GetScalarRange());
  pm->SetScalarVisibility(anActor->GetMapper()->GetScalarVisibility());
  pm->SetLookupTable(anActor->GetMapper()->GetLookupTable());
  pm->SetScalarMode(anActor->GetMapper()->GetScalarMode());

  if(pm->GetScalarMode() == VTK_SCALAR_MODE_USE_POINT_FIELD_DATA
     || pm->GetScalarMode() == VTK_SCALAR_MODE_USE_CELL_FIELD_DATA) {
    if(anActor->GetMapper()->GetArrayAccessMode() == VTK_GET_ARRAY_BY_ID) {
      pm->ColorByArrayComponent(anActor->GetMapper()->GetArrayId(),
                                anActor->GetMapper()->GetArrayComponent());
    } else {
      pm->ColorByArrayComponent(anActor->GetMapper()->GetArrayName(),
                                anActor->GetMapper()->GetArrayComponent());
    }
  }

  points = pd->GetPoints();
  pntData = pd->GetPointData();
  normals = pntData->GetNormals();
  tcoords = pntData->GetTCoords();
  colors = pm->MapScalars(1.0);

  // write out polys if any
  if(pd->GetNumberOfPolys() > 0) {
    WriteShapeBegin(anActor, fp, pd, pntData, colors);
    fprintf(fp, "          geometry IndexedFaceSet {\n");
    // two sided lighting ? for now assume it is on
    fprintf(fp, "            solid FALSE\n");

    if(!pointDataWritten) {
      this->WritePointData(points, normals, tcoords, colors, fp);
      pointDataWritten = 1;
    } else {
      fprintf(fp, "            coord  USE VTKcoordinates\n");

      if(normals) {
        fprintf(fp, "            normal  USE VTKnormals\n");
      }

      if(tcoords) {
        fprintf(fp, "            texCoord  USE VTKtcoords\n");
      }

      if(colors) {
        fprintf(fp, "            color  USE VTKcolors\n");
      }
    }

    fprintf(fp, "            coordIndex  [\n");

    cells = pd->GetPolys();

    for(cells->InitTraversal(); cells->GetNextCell(npts, indx);) {
      fprintf(fp, "              ");

      for(i = 0; i < npts; i++) {
        // treating vtkIdType as int
        fprintf(fp, "%i, ", static_cast<int>(indx[i]));
      }

      fprintf(fp, "-1,\n");
    }

    fprintf(fp, "            ]\n");
    fprintf(fp, "          }\n");
    WriteShapeEnd(fp);
  }

  // write out tstrips if any
  if(pd->GetNumberOfStrips() > 0) {
    WriteShapeBegin(anActor, fp, pd, pntData, colors);
    fprintf(fp, "          geometry IndexedFaceSet {\n");

    if(!pointDataWritten) {
      this->WritePointData(points, normals, tcoords, colors, fp);
      pointDataWritten = 1;
    } else {
      fprintf(fp, "            coord  USE VTKcoordinates\n");

      if(normals) {
        fprintf(fp, "            normal  USE VTKnormals\n");
      }

      if(tcoords) {
        fprintf(fp, "            texCoord  USE VTKtcoords\n");
      }

      if(colors) {
        fprintf(fp, "            color  USE VTKcolors\n");
      }
    }

    fprintf(fp, "            coordIndex  [\n");
    cells = pd->GetStrips();

    for(cells->InitTraversal(); cells->GetNextCell(npts, indx);) {
      for(i = 2; i < npts; i++) {
        if(i % 2) {
          i1 = i - 1;
          i2 = i - 2;
        } else {
          i1 = i - 2;
          i2 = i - 1;
        }

        // treating vtkIdType as int
        fprintf(fp, "              %i, %i, %i, -1,\n",
                static_cast<int>(indx[i1]), static_cast<int>(indx[i2]),
                static_cast<int>(indx[i]));
      }
    }

    fprintf(fp, "            ]\n");
    fprintf(fp, "          }\n");
    WriteShapeEnd(fp);
  }

  // write out lines if any
  if(pd->GetNumberOfLines() > 0) {
    WriteShapeBegin(anActor, fp, pd, pntData, colors);
    fprintf(fp, "          geometry IndexedLineSet {\n");

    if(!pointDataWritten) {
      this->WritePointData(points, nullptr, nullptr, colors, fp);
    } else {
      fprintf(fp, "            coord  USE VTKcoordinates\n");

      if(colors) {
        fprintf(fp, "            color  USE VTKcolors\n");
      }
    }

    fprintf(fp, "            coordIndex  [\n");

    cells = pd->GetLines();

    for(cells->InitTraversal(); cells->GetNextCell(npts, indx);) {
      fprintf(fp, "              ");

      for(i = 0; i < npts; i++) {
        // treating vtkIdType as int
        fprintf(fp, "%i, ", static_cast<int>(indx[i]));
      }

      fprintf(fp, "-1,\n");
    }

    fprintf(fp, "            ]\n");
    fprintf(fp, "          }\n");
    WriteShapeEnd(fp);
  }

  // write out verts if any
  if(pd->GetNumberOfVerts() > 0) {
    WriteShapeBegin(anActor, fp, pd, pntData, colors);
    fprintf(fp, "          geometry PointSet {\n");
    cells = pd->GetVerts();
    fprintf(fp, "            coord Coordinate {");
    fprintf(fp, "              point [");

    for(cells->InitTraversal(); cells->GetNextCell(npts, indx);) {
      fprintf(fp, "              ");

      for(i = 0; i < npts; i++) {
        p = points->GetPoint(indx[i]);
        fprintf(fp, "              %g %g %g,\n", p[0], p[1], p[2]);
      }
    }

    fprintf(fp, "              ]\n");
    fprintf(fp, "            }\n");

    if(colors) {
      fprintf(fp, "            color Color {");
      fprintf(fp, "              color [");

      for(cells->InitTraversal(); cells->GetNextCell(npts, indx);) {
        fprintf(fp, "              ");

        for(i = 0; i < npts; i++) {
          c = colors->GetPointer(4 * indx[i]);
          fprintf(fp, "           %g %g %g,\n", c[0] / 255.0, c[1] / 255.0,
                  c[2] / 255.0);
        }
      }

      fprintf(fp, "              ]\n");
      fprintf(fp, "            }\n");
    }

    fprintf(fp, "          }\n");
    WriteShapeEnd(fp);
  }

  fprintf(fp, "      ]\n"); // close the original transforms children
  fprintf(fp, "    }\n"); // close the original transform

  pm->Delete();
}

TTKWRLEXPORTER_EXPORT void vtkVRMLExporter::WriteData() {

  vtkRenderer *ren;
  vtkActorCollection *ac;
  vtkActor *anActor, *aPart;
  vtkLightCollection *lc;
  vtkLight *aLight;
  // //   vtkCamera *cam;
  //   double *tempd;
  FILE *fp;

  // make sure the user specified a FileName or FilePointer
  if(!this->FilePointer && (this->FileName == nullptr)) {
    vtkErrorMacro(<< "Please specify FileName to use");
    return;
  }

  // Always pick the first renderer
  // first make sure there is only one renderer in this rendering window
  // if (this->RenderWindow->GetRenderers()->GetNumberOfItems() > 1)
  //  {
  //  vtkErrorMacro(<< "VRML files only support one renderer per window.");
  //  return;
  //  }

  // get the renderer
  ren = this->RenderWindow->GetRenderers()->GetFirstRenderer();

  // make sure it has at least one actor
  if(ren->GetActors()->GetNumberOfItems() < 1) {
    vtkErrorMacro(<< "no actors found for writing VRML file.");
    return;
  }

  // try opening the files
  if(!this->FilePointer) {
    fp = fopen(this->FileName, "w");

    if(!fp) {
      vtkErrorMacro(<< "unable to open VRML file " << this->FileName);
      return;
    }
  } else {
    fp = this->FilePointer;
  }

  //
  //  Write header
  //
  vtkDebugMacro("Writing VRML file");
  fprintf(fp, "#VRML V2.0 utf8\n");
  fprintf(fp, "# VRML file written by the visualization toolkit\n\n");

  // Start write the Background
  double background[3];
  ren->GetBackground(background);
  fprintf(fp, "    Background {\n ");
  fprintf(fp, "   skyColor [%f %f %f, ]\n", background[0], background[1],
          background[2]);
  fprintf(fp, "    }\n ");
  // End of Background

  // BUG fix
  // do the camera
  //   cam = ren->GetActiveCamera();
  //   fprintf(fp, "    Viewpoint\n      {\n      fieldOfView %f\n",
  //           cam->GetViewAngle()*vtkMath::Pi() / 180.0);
  //   fprintf(fp, "      position %f %f %f\n", cam->GetPosition()[0],
  //           cam->GetPosition()[1], cam->GetPosition()[2]);
  //   fprintf(fp, "      description \"Default View\"\n");
  //   tempd = cam->GetOrientationWXYZ();
  //   fprintf(fp, "      orientation %g %g %g %g\n      }\n", tempd[1],
  //   tempd[2],
  //           tempd[3], tempd[0]*vtkMath::Pi() / 180.0);
  //
  //   // do the lights first the ambient then the others
  //   fprintf(fp,
  //           "    NavigationInfo {\n      type [\"EXAMINE\",\"FLY\"]\n speed
  //           %f\n",
  //           this->Speed);
  //
  //   if (ren->GetLights()->GetNumberOfItems() == 0){
  //     fprintf(fp, "      headlight TRUE}\n\n");
  //   }
  //   else{
  //     fprintf(fp, "      headlight FALSE}\n\n");
  //   }
  // end of BUG fix

  fprintf(
    fp, "DirectionalLight { ambientIntensity 1 intensity 0 # ambient light\n");
  fprintf(fp, "      color %f %f %f }\n\n", ren->GetAmbient()[0],
          ren->GetAmbient()[1], ren->GetAmbient()[2]);

  // make sure we have a default light
  // if we dont then use a headlight
  lc = ren->GetLights();
  vtkCollectionSimpleIterator lsit;

  for(lc->InitTraversal(lsit); (aLight = lc->GetNextLight(lsit));) {
    this->WriteALight(aLight, fp);
  }

  // do the actors now
  ac = ren->GetActors();
  vtkAssemblyPath *apath;
  vtkCollectionSimpleIterator ait;

  for(ac->InitTraversal(ait); (anActor = ac->GetNextActor(ait));) {
    for(anActor->InitPathTraversal(); (apath = anActor->GetNextPath());) {
      aPart = static_cast<vtkActor *>(apath->GetLastNode()->GetViewProp());
      this->WriteAnActor(aPart, fp);
    }
  }

  if(!this->FilePointer) {
    fclose(fp);
  }
}

TTKWRLEXPORTER_EXPORT void
  vtkVRMLExporter::WritePointData(vtkPoints *points,
                                  vtkDataArray *normals,
                                  vtkDataArray *tcoords,
                                  vtkUnsignedCharArray *colors,
                                  FILE *fp) {

  double *p;
  unsigned char *c;

  // write out the points
  fprintf(fp, "            coord DEF VTKcoordinates Coordinate {\n");
  fprintf(fp, "              point [\n");
  for(int i = 0; i < points->GetNumberOfPoints(); i++) {
    p = points->GetPoint(i);
    fprintf(fp, "              %g %g %g,\n", p[0], p[1], p[2]);
  }
  fprintf(fp, "              ]\n");
  fprintf(fp, "            }\n");

  // write out the point data
  if(normals) {

    fprintf(fp, "            normal DEF VTKnormals Normal {\n");
    fprintf(fp, "              vector [\n");
    for(int i = 0; i < normals->GetNumberOfTuples(); i++) {
      p = normals->GetTuple(i);
      fprintf(fp, "           %g %g %g,\n", p[0], p[1], p[2]);
    }
    fprintf(fp, "            ]\n");
    fprintf(fp, "          }\n");
  }

  // write out the point data
  if(tcoords) {
    fprintf(fp, "            texCoord DEF VTKtcoords TextureCoordinate {\n");
    fprintf(fp, "              point [\n");
    for(int i = 0; i < tcoords->GetNumberOfTuples(); i++) {
      p = tcoords->GetTuple(i);
      fprintf(fp, "           %g %g,\n", p[0], p[1]);
    }
    fprintf(fp, "            ]\n");
    fprintf(fp, "          }\n");

    // BUG fix here.
    if(ttkWRLExporterPolyData_) {
      fprintf(fp, "          texCoordIndex[\n");
      vtkCellArray *cells = ttkWRLExporterPolyData_->GetPolys();
      vtkIdType npts = 0;
#ifdef VTK_CELL_ARRAY_V2
      vtkIdType const *indx = nullptr;
#else
      vtkIdType *indx = nullptr;
#endif
      for(cells->InitTraversal(); cells->GetNextCell(npts, indx);) {
        fprintf(fp, "            ");
        for(int i = 0; i < npts; i++) {
          fprintf(fp, "%i, ", static_cast<int>(indx[i]));
        }
        fprintf(fp, "-1,\n");
      }
      fprintf(fp, "          ]\n");
    }
    // end of BUG fix here.
  }

  // write out the point data
  if(colors) {
    fprintf(fp, "            color DEF VTKcolors Color {\n");
    fprintf(fp, "              color [\n");
    for(int i = 0; i < colors->GetNumberOfTuples(); i++) {
      c = colors->GetPointer(4 * i);
      fprintf(
        fp, "           %g %g %g,\n", c[0] / 255.0, c[1] / 255.0, c[2] / 255.0);
    }
    fprintf(fp, "            ]\n");
    fprintf(fp, "          }\n");
  }
}

/// @endcond
