
// ./anisofracture -test 8

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string helper_output;
helper_output = "output/RIG_fleece_fibres";
std::cout << "Attempt 1" << std::endl;
sim.output_dir.path = helper_output;
sim.end_frame = 7;
T frameRate = 1;
sim.step.frame_dt = (T)1 / frameRate;
sim.gravity = -9.81 * TV::Unit(1);
sim.step.max_dt = 1e-3;
sim.symplectic = true;
sim.verbose = false;
sim.transfer_scheme = MpmSimulationBase<T, dim>::FLIP_blend_PIC;

sim.dump_F_for_meshing = true;
T particle_per_cell = 50;
sim.rpic_damping_iteration = 0;

// ****************************************************************************
// Loop for all files
// ****************************************************************************

std::ifstream file;
file.open("../../Data/TetMesh/mesh_files/fibre_direction.csv");
std::vector<std::vector<double>> directions;
std::string line;

if (!file.is_open())
  throw std::runtime_error("Could not open file");

while (std::getline(file, line)) {
  std::stringstream s(line);
  std::vector<double> row;
  std::string val;
  while (std::getline(s, val, ',')) {
    row.push_back(std::stod(val));
  }

  directions.push_back(row);
}

file.close();

for (size_t file_index = 0; file_index < 98; ++file_index) {

  // Construct the file path based on the index
  std::string meshFilePath =
      "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/file_" +
      std::to_string(file_index + 1) + ".mesh";

  FILE *file = std::fopen(meshFilePath.c_str(), "r");

  if (file != nullptr) {
    std::fclose(file);

    // Get the direction corresponding to this file
    std::vector<double> helper_fiber_vector = directions[file_index];

    bool helper_anisotropic =
        std::all_of(helper_fiber_vector.begin(), helper_fiber_vector.end(),
                    [](int i) { return i == 0; });

    TV helper_fiber(helper_fiber_vector[0], helper_fiber_vector[1],
                    helper_fiber_vector[2]);

    // Print the values to verify
    // std::cout << "File index: " << file_index << " Helper Isotropic: " <<
    // helper_isotropic
    // << " Helper Fiber: [" << helper_fiber[0] << ", " << helper_fiber[1] << ",
    // " << helper_fiber[2] << "]" << std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Parameters that can be modified ////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::map<std::string, double> parameters;

    // Open the parameters file for reading
    std::ifstream paramFile("../cmaes/parameters.txt");

    if (!paramFile.is_open()) {
      std::cerr << "Error: Could not open parameters file." << std::endl;
    }

    std::string line;
    while (std::getline(paramFile, line)) {
      size_t separatorPos = line.find(":");
      if (separatorPos != std::string::npos) {
        std::string paramName = line.substr(0, separatorPos);
        double paramValue = std::stod(line.substr(separatorPos + 1));
        parameters[paramName] = paramValue;
      }
    }

    // Close the parameters file
    paramFile.close();

    // Access the parameters by their names
    T Youngs = parameters["Youngs"];
    T nu = parameters["nu"];
    T rho = parameters["rho"];
    T percentage = parameters["percentage"];
    T eta = parameters["eta"];
    double fiberstifness = parameters["fiber"];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// End ////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    sim.cfl = 0.4;
    T helper_alpha = 0;
    T residual_stress = 0.01;

    // std::string filename = "TetMesh/RIG_test_holes.mesh";
    MpmParticleHandleBase<T, dim> particles_handle =
        init_helper.sampleFromTetWildFile(meshFilePath, rho);
    T total_volume = particles_handle.total_volume;
    T particle_count = particles_handle.particle_range.length();
    T per_particle_volume = total_volume / particle_count;
    sim.dx = std::pow(particle_per_cell * per_particle_volume, (T)1 / (T)3);

    /*
    if (1) {
      StdVector<TV> samples;
      StdVector<Vector<int, 4>> indices;
      std::string absolute_path = DataDir().absolutePath(meshFilePath);
      readTetMeshTetWild(absolute_path, samples, indices);
      sim.output_dir.createPath();
      std::string vtk_path = DataDir().path + "/../Projects/anisofracture/" +
                             sim.output_dir.path + "/tet2.vtk";
      writeTetmeshVtk(vtk_path, samples, indices);
    }
    */

    QRAnisotropic<T, dim> model(Youngs, nu, helper_anisotropic);
    StdVector<TV> a_0;
    StdVector<T> alphas;
    TV a_1, a_2;
    a_1 = helper_fiber;
    a_1.normalize();
    a_0.push_back(a_1);
    alphas.push_back(helper_alpha);
    T l0 = 0.5 * sim.dx;

    T zeta = 1;
    bool allow_damage = true;
    // model.scaleFiberStiffness(0, 2);
    model.scaleFiberStiffness(0, fiberstifness);
    particles_handle.addFBasedMpmForceWithAnisotropicPhaseField(
        a_0, alphas, percentage, l0, model, eta, zeta, allow_damage,
        residual_stress);
  } else {
    std::cout << "File does not exist: " << meshFilePath << std::endl;
    continue;
  }
}

// ****************************************************************************
// Ground Plane
// ****************************************************************************
TV ground_origin = TV(2, 1.817, 2);
TV ground_normal(0, 1, 0);
HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
AnalyticCollisionObject<T, dim>
    ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
ground_object.setFriction(1);
init_helper.addAnalyticCollisionObject(ground_object);

double sphere_radius = 0.0078;

// ****************************************************************************
// Left Hand
// ****************************************************************************

auto left1Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
  T t = time;
  TV translation = TV(0, 0, -.006 * time);
  TV translation_velocity(0, 0, -.006);
  object.setTranslation(translation, translation_velocity);
};

//////////
// 1
/////////
{
  CappedCylinder<T, dim> cylinder11(sphere_radius, 1, Vector<T, 4>(1, 0, 0, 0),
                                    TV(1.975, 1.825, 1.975)); // 1.907
  HalfSpace<T, dim> board1(TV(0, 1, 0), TV(0, 1, 0));
  DifferenceLevelSet<T, dim> cutsphere11;
  cutsphere11.add(cylinder11, board1);
  AnalyticCollisionObject<T, dim> leftObject1(
      left1Transform, cutsphere11, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject1);
}

//////////
// 2
/////////
{

  CappedCylinder<T, dim> cylinder12(
      sphere_radius, 0.5, Vector<T, 4>(1, 0, 0, 0), TV(2.017, 1.825, 1.975));
  HalfSpace<T, dim> board1(TV(0, 1, 0), TV(0, 1, 0));
  DifferenceLevelSet<T, dim> cutsphere12;
  cutsphere12.add(cylinder12, board1);
  AnalyticCollisionObject<T, dim> leftObject2(
      left1Transform, cutsphere12, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject2);
}

//////////
// 3
/////////
{

  CappedCylinder<T, dim> cylinder13(
      sphere_radius, 0.5, Vector<T, 4>(1, 0, 0, 0), TV(2.065, 1.825, 1.983));
  HalfSpace<T, dim> board1(TV(0, 1, 0), TV(0, 1, 0));
  DifferenceLevelSet<T, dim> cutsphere13;
  cutsphere13.add(cylinder13, board1);
  AnalyticCollisionObject<T, dim> leftObject3(
      left1Transform, cutsphere13, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject3);
}
// ****************************************************************************
// Rigth Hand
// ****************************************************************************

auto right1Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
  TV translation = TV(0, 0, .006 * time);
  TV translation_velocity(0, 0, .006);
  object.setTranslation(translation, translation_velocity);
};
//////////
// 1
/////////
{

  CappedCylinder<T, dim> cylinder21(
      sphere_radius, 0.5, Vector<T, 4>(1, 0, 0, 0), TV(1.972, 1.825, 2.01));
  HalfSpace<T, dim> board2(TV(0, 1, 0), TV(0, 1, 0));
  DifferenceLevelSet<T, dim> cutsphere21;
  cutsphere21.add(cylinder21, board2);
  AnalyticCollisionObject<T, dim> rightObject1(
      right1Transform, cutsphere21, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject1);
}
//////////
// 2
/////////
{

  CappedCylinder<T, dim> cylinder22(
      sphere_radius, 0.5, Vector<T, 4>(1, 0, 0, 0), TV(2.015, 1.825, 2.01));
  HalfSpace<T, dim> board2(TV(0, 1, 0), TV(0, 1, 0));
  DifferenceLevelSet<T, dim> cutsphere22;
  cutsphere22.add(cylinder22, board2);
  AnalyticCollisionObject<T, dim> rightObject2(
      right1Transform, cutsphere22, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject2);
}
//////////
// 3
/////////
{

  CappedCylinder<T, dim> cylinder23(sphere_radius, 0.5,
                                    Vector<T, 4>(1, 0, 0, 0),
                                    TV(2.065, 1.825, 2.015)); // End frame 2.019
  HalfSpace<T, dim> board2(TV(0, 1, 0), TV(0, 1, 0));
  DifferenceLevelSet<T, dim> cutsphere23;
  cutsphere23.add(cylinder23, board2);
  AnalyticCollisionObject<T, dim> rightObject3(
      right1Transform, cutsphere23, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject3);
}
// init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx,
// AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for
// SPGrid.
