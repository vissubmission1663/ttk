#include <ttkBottleneckDistance.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkBottleneckDistance);

ttkBottleneckDistance::ttkBottleneckDistance() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

int ttkBottleneckDistance::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkBottleneckDistance::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

// Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
int augmentDiagrams(const ttk::DiagramType &diag0,
                    const ttk::DiagramType &diag1,
                    const std::vector<matchingTuple> &matchings,
                    vtkUnstructuredGrid *const vtu0,
                    vtkUnstructuredGrid *const vtu1) {

  if(matchings.empty()) {
    return 0;
  }

  vtkNew<vtkIntArray> matchingIds0{};
  matchingIds0->SetName("MatchingIdentifier");
  matchingIds0->SetNumberOfComponents(1);
  matchingIds0->SetNumberOfTuples(vtu0->GetNumberOfCells());
  vtu0->GetCellData()->AddArray(matchingIds0);

  vtkNew<vtkIntArray> matchingIds1{};
  matchingIds1->SetName("MatchingIdentifier");
  matchingIds1->SetNumberOfComponents(1);
  matchingIds1->SetNumberOfTuples(vtu1->GetNumberOfCells());
  vtu1->GetCellData()->AddArray(matchingIds1);

  // Unaffected by default
  matchingIds0->Fill(-1);
  matchingIds1->Fill(-1);

  // Affect bottleneck matchings
  for(size_t i = 0; i < matchings.size(); ++i) {
    const auto &t = matchings[i];
    matchingIds0->SetTuple1(std::get<0>(t), i);
    matchingIds1->SetTuple1(std::get<1>(t), i);
  }

  return 1;
}

int translateDiagram(vtkUnstructuredGrid *output,
                     vtkUnstructuredGrid *input,
                     const double spacing) {
  vtkNew<vtkTransform> tr{};
  tr->Translate(0, 0, spacing);
  vtkNew<vtkTransformFilter> trf{};
  trf->SetTransform(tr);
  trf->SetInputData(input);
  trf->Update();
  output->ShallowCopy(trf->GetOutputDataObject(0));

  return 1;
}

int generateMatchings(vtkUnstructuredGrid *const outputCT3,
                      const ttk::DiagramType &diagram1,
                      const ttk::DiagramType &diagram2,
                      const std::vector<matchingTuple> &matchings,
                      const double spacing,
                      const bool is2D0,
                      const bool is2D1) {

  vtkNew<vtkUnstructuredGrid> vtu{};

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * matchings.size());
  vtu->SetPoints(points);

  vtkNew<vtkDoubleArray> costs{};
  costs->SetName("Cost");
  costs->SetNumberOfComponents(1);
  costs->SetNumberOfTuples(matchings.size());
  vtu->GetCellData()->AddArray(costs);

  vtkNew<vtkIntArray> matchingIds{};
  matchingIds->SetName("MatchingIdentifier");
  matchingIds->SetNumberOfComponents(1);
  matchingIds->SetNumberOfTuples(matchings.size());
  vtu->GetCellData()->AddArray(matchingIds);

  // Build matchings.
  for(size_t i = 0; i < matchings.size(); ++i) {
    const auto &t = matchings[i];
    const auto n1 = std::get<0>(t);
    const auto n2 = std::get<1>(t);

    const auto &pair0 = diagram1[n1];
    const auto &pair1 = diagram2[n2];

    const auto pairPoint = [](const ttk::PairTuple &pair, const bool is2D,
                              const double zval) -> std::array<double, 3> {
      if(is2D) {
        return {std::get<6>(pair), std::get<10>(pair), zval};
      } else {
        return {std::get<11>(pair), std::get<12>(pair), std::get<13>(pair)};
      }
    };

    const auto p0 = pairPoint(pair0, is2D0, -spacing / 2.0);
    points->SetPoint(2 * i + 0, p0.data());
    const auto p1 = pairPoint(pair1, is2D1, spacing / 2.0);
    points->SetPoint(2 * i + 1, p1.data());

    std::array<vtkIdType, 2> ids{
      2 * static_cast<vtkIdType>(i) + 0,
      2 * static_cast<vtkIdType>(i) + 1,
    };
    vtu->InsertNextCell(VTK_LINE, 2, ids.data());

    costs->SetTuple1(i, std::get<2>(t));
    matchingIds->SetTuple1(i, i);
  }

  outputCT3->ShallowCopy(vtu);

  return 1;
}

int ttkBottleneckDistance::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto outputDiagrams = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto outputMatchings = vtkUnstructuredGrid::GetData(outputVector, 1);

  // Wrap
  // this->setWrapper(this);
  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0]);
  std::vector<vtkUnstructuredGrid *> inputDiags{};

  if(blocks == nullptr) {
    this->printErr("No input diagrams");
    return 0;
  }

  for(size_t i = 0; i < blocks->GetNumberOfBlocks(); ++i) {
    const auto diag = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
    if(diag != nullptr) {
      inputDiags.emplace_back(diag);
    }
  }

  if(inputDiags.size() < 2) {
    this->printErr("Less than two input diagrams");
    return 0;
  }
  if(inputDiags.size() > 2) {
    this->printWrn("More than two input diagrams: "
                   + std::to_string(inputDiags.size()));
  }

  const auto coords0 = vtkFloatArray::SafeDownCast(
    inputDiags[0]->GetPointData()->GetArray("Coordinates"));
  const auto coords1 = vtkFloatArray::SafeDownCast(
    inputDiags[1]->GetPointData()->GetArray("Coordinates"));

  const bool is2D0 = coords0 != nullptr;
  const bool is2D1 = coords1 != nullptr;

  // Call package
  int status = 0;

  ttk::DiagramType diagram0{}, diagram1{};

  status = VTUToDiagram(diagram0, inputDiags[0], *this);
  if(status < 0) {
    this->printErr("Could not extract diagram from first input data-set");
    return 0;
  }

  status = VTUToDiagram(diagram1, inputDiags[1], *this);
  if(status < 0) {
    this->printErr("Could not extract diagram from second input data-set");
    return 0;
  }

  this->setCTDiagram1(&diagram0);
  this->setCTDiagram2(&diagram1);

  this->setWasserstein(this->WassersteinMetric);
  this->setAlgorithm(this->DistanceAlgorithm);
  this->setPVAlgorithm(this->PVAlgorithm);

  // Empty matchings.
  std::vector<matchingTuple> matchings;
  this->setOutputMatchings(&matchings);

  // Exec.
  status = this->execute<double>(this->UsePersistenceMetric);
  if(status != 0) {
    this->printErr("Base layer failed with error status "
                   + std::to_string(status));
    return 0;
  }

  // Generate matchings
  if(this->UseOutputMatching) {
    status = generateMatchings(outputMatchings, diagram0, diagram1, matchings,
                               this->Spacing, is2D0, is2D1);

    if(status != 1) {
      this->printErr("Could not compute matchings");
      return 0;
    }
  }

  // Translate diagrams
  vtkNew<vtkUnstructuredGrid> vtu0{}, vtu1{};
  if(this->UseGeometricSpacing) {
    translateDiagram(vtu0, inputDiags[0], -this->Spacing / 2.0);
    translateDiagram(vtu1, inputDiags[1], this->Spacing / 2.0);
  } else {
    vtu0->ShallowCopy(inputDiags[0]);
    vtu1->ShallowCopy(inputDiags[1]);
  }

  // Add matchings infos on diagrams
  status = augmentDiagrams(diagram0, diagram1, matchings, vtu0, vtu1);
  if(status != 1) {
    this->printErr("Could not augment diagrams");
    return 0;
  }

  // Set output.
  outputDiagrams->SetNumberOfBlocks(2);
  outputDiagrams->SetBlock(0, vtu0);
  outputDiagrams->SetBlock(1, vtu1);

  return 1;
}
