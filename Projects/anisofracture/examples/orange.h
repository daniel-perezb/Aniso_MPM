// Test with an anisotropic top part with fibres pointing straigth up
// Borrom has an isotropic face because of contaminants forcing fibres together
// ./anisofracture -test 8
// #include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cstdio>

std::string helper_output;
helper_output = "output/RIG_fleece_fibres";
std::cout << "Attempt 1" << std::endl;
sim.output_dir.path = helper_output;
sim.end_frame = 150;
T frameRate = 10;
sim.step.frame_dt = (T)1 / frameRate;
sim.gravity = .0003 * TV::Unit(1);
sim.step.max_dt = 1e-3;
sim.symplectic = true;
sim.verbose = false;
sim.cfl = 0.4;
sim.transfer_scheme = MpmSimulationBase<T, dim>::FLIP_blend_PIC;
sim.flip_pic_ratio = 0; // FULL PIC for damping
sim.dump_F_for_meshing = true;
T particle_per_cell = 50;
sim.rpic_damping_iteration = 0;

// ****************************************************************************
// Top
// ****************************************************************************
// if (1) {
//   T Youngs = 10000;
//   T nu = 0.4;
//   T rho = 500;
//   T helper_isotropic = false;
//   TV helper_fiber = TV(0, 1, 0); // fibres going straigth up
//   T helper_alpha = -1;

//   std::string filename = "TetMesh/RIG_test1.mesh";
//   MpmParticleHandleBase<T, dim> particles_handle =
//       init_helper.sampleFromTetWildFile(filename, rho);
//   T total_volume = particles_handle.total_volume;
//   T particle_count = particles_handle.particle_range.length();
//   T per_particle_volume = total_volume / particle_count;
//   sim.dx = std::pow(particle_per_cell * per_particle_volume, (T)1 / (T)3);
//   if (1) {
//     StdVector<TV> samples;
//     StdVector<Vector<int, 4>> indices;
//     std::string absolute_path = DataDir().absolutePath(filename);
//     readTetMeshTetWild(absolute_path, samples, indices);
//     sim.output_dir.createPath();
//     std::string vtk_path = DataDir().path + "/../Projects/anisofracture/" +
//                            sim.output_dir.path + "/tet1.vtk";
//     writeTetmeshVtk(vtk_path, samples, indices);
//   }

//   QRAnisotropic<T, dim> model(Youngs, nu, helper_isotropic);
//   StdVector<TV> a_0;
//   StdVector<T> alphas;
//   TV a_1, a_2;
//   a_1 = helper_fiber;
//   a_1.normalize();
//   a_0.push_back(a_1);
//   alphas.push_back(helper_alpha);
//   T percentage = 0.15;
//   T l0 = 0.5 * sim.dx;
//   T eta = 0.01;
//   T zeta = 1;
//   bool allow_damage = true;
//   T residual_stress = 0.001;
//   // model.scaleFiberStiffness(0, 2);
//   particles_handle.addFBasedMpmForceWithAnisotropicPhaseField(
//       a_0, alphas, percentage, l0, model, eta, zeta, allow_damage,
//       residual_stress);

//   // SnowPlasticity<T> p(0, 1, 0.5);
//   // particles_handle.addPlasticity(model, p, "F");
//   std::cout << "Particle count: " << sim.particles.count << std::endl;
// }

// ****************************************************************************
// Bottom
// ****************************************************************************
// if (1) {
//   T Youngs = 20000;
//   T nu = 0.4;
//   T rho = 500;
//   T helper_isotropic = true;
//   TV helper_fiber = TV(1, 1, 1); // No fibre direction (Isotropic model)
//   T helper_alpha = 0;

//   std::string filename = "TetMesh/RIG_test_holes.mesh";
//   MpmParticleHandleBase<T, dim> particles_handle =
//       init_helper.sampleFromTetWildFile(filename, rho);
//   T total_volume = particles_handle.total_volume;
//   T particle_count = particles_handle.particle_range.length();
//   T per_particle_volume = total_volume / particle_count;
//   sim.dx = std::pow(particle_per_cell * per_particle_volume, (T)1 / (T)3);
//   if (1) {
//     StdVector<TV> samples;
//     StdVector<Vector<int, 4>> indices;
//     std::string absolute_path = DataDir().absolutePath(filename);
//     readTetMeshTetWild(absolute_path, samples, indices);
//     sim.output_dir.createPath();
//     std::string vtk_path = DataDir().path + "/../Projects/anisofracture/" +
//                            sim.output_dir.path + "/tet2.vtk";
//     writeTetmeshVtk(vtk_path, samples, indices);
//   }

//   QRAnisotropic<T, dim> model(Youngs, nu, helper_isotropic);
//   StdVector<TV> a_0;
//   StdVector<T> alphas;
//   TV a_1, a_2;
//   a_1 = helper_fiber;
//   a_1.normalize();
//   a_0.push_back(a_1);
//   alphas.push_back(helper_alpha);
//   T percentage = 999;
//   T l0 = 0.5 * sim.dx;
//   T eta = 0.1;
//   T zeta = 1;
//   bool allow_damage = true;
//   T residual_stress = 0.005;
//   // model.scaleFiberStiffness(0, 2);
//   particles_handle.addFBasedMpmForceWithAnisotropicPhaseField(
//       a_0, alphas, percentage, l0, model, eta, zeta, allow_damage,
//       residual_stress);
// }

// ****************************************************************************
// Loop for all files
// ****************************************************************************

std::ifstream file;
file.open("/ziran2020/Data/TetMesh/cutfiles/fibre_direction.csv");
std::vector<std::vector<double>> directions;
std::string line;

if(!file.is_open()) throw std::runtime_error("Could not open file");

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
for (int rowd =0 ; rowd < 12; rowd++){
  for (double num : directions[rowd]) {
      std::cout << num << " " ;
    }
    std::cout << std::endl;
}
   
for (size_t file_index = 0; file_index < 12; ++file_index) {
    
    // Construct the file path based on the index
    std::string meshFilePath = "/ziran2020/Data/TetMesh/mesh_files/file_" + std::to_string(file_index + 1) + ".mesh";

    FILE* file = std::fopen(meshFilePath.c_str(), "r");

    if (file != nullptr) {
        std::fclose(file);
    
    // Get the direction corresponding to this file
    std::vector<double> helper_fiber_vector = directions[file_index];

    bool helper_isotropic = std::all_of(helper_fiber_vector.begin(), helper_fiber_vector.end(), [](int i) { return i == 0; });

    TV helper_fiber(helper_fiber_vector[0], helper_fiber_vector[1], helper_fiber_vector[2]);
    
    // Print the values to verify
    std::cout << "File index: " << file_index << " Helper Isotropic: " << helper_isotropic 
    << " Helper Fiber: [" << helper_fiber[0] << ", " << helper_fiber[1] << ", " << helper_fiber[2] << "]" << std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Parameters that can be modified //// 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::map<std::string, double> parameters;

    // Open the parameters file for reading
    std::ifstream paramFile("parameters.txt");
    
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
    T residual_stress = parameters["residual_stress"];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// End //// 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 

    T helper_alpha = 0;

    // std::string filename = "TetMesh/RIG_test_holes.mesh";
    MpmParticleHandleBase<T, dim> particles_handle =
      init_helper.sampleFromTetWildFile(meshFilePath, rho);
    T total_volume = particles_handle.total_volume;
    T particle_count = particles_handle.particle_range.length();
    T per_particle_volume = total_volume / particle_count;
    sim.dx = std::pow(particle_per_cell * per_particle_volume, (T)1 / (T)3);
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

  QRAnisotropic<T, dim> model(Youngs, nu, helper_isotropic);
  StdVector<TV> a_0;
  StdVector<T> alphas;
  TV a_1, a_2;
  a_1 = helper_fiber;
  a_1.normalize();
  a_0.push_back(a_1);
  alphas.push_back(helper_alpha);
  T percentage = 999;
  T l0 = 0.5 * sim.dx;
  T eta = 0.1;
  T zeta = 1;
  bool allow_damage = true;
  // model.scaleFiberStiffness(0, 2);
  particles_handle.addFBasedMpmForceWithAnisotropicPhaseField(
    a_0, alphas, percentage, l0, model, eta, zeta, allow_damage,
    residual_stress);
} 
else {
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

double sphere_radius = 0.01;

// ****************************************************************************
// Left Hand
// ****************************************************************************

{
  auto left1Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    T t = time;
    TV translation = TV(0, 0, -.05 * time);
    TV translation_velocity(0, 0, -.05);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> leftLS1(TV(1.998, 1.83, 2.02), sphere_radius);
  AnalyticCollisionObject<T, dim> leftObject(
      left1Transform, leftLS1, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject);
}

{
  auto left2Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    T t = time;
    TV translation = TV(0, 0, -.05 * time);
    TV translation_velocity(0, 0, -.05);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> leftLS2(TV(2.028, 1.83, 2.02), sphere_radius);
  AnalyticCollisionObject<T, dim> leftObject(
      left2Transform, leftLS2, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject);
}

{
  auto left3Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    T t = time;
    TV translation = TV(0, 0, -.05 * time);
    TV translation_velocity(0, 0, -.05);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> leftLS3(TV(2.058, 1.83, 2.02), sphere_radius);
  AnalyticCollisionObject<T, dim> leftObject(
      left3Transform, leftLS3, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject);
}

// ****************************************************************************
// Rigth Hand
// ****************************************************************************

{
  auto right1Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    TV translation = TV(0, 0, .05 * time);
    TV translation_velocity(0, 0, .05);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> sphere1(TV(1.998, 1.83, 2.047), sphere_radius);
  AnalyticCollisionObject<T, dim> rightObject(
      right1Transform, sphere1, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject);
}

{
  auto right2Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    TV translation = TV(0, 0, .05 * time);
    TV translation_velocity(0, 0, .05);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> sphere2(TV(2.028, 1.83, 2.047), sphere_radius);
  AnalyticCollisionObject<T, dim> rightObject(
      right2Transform, sphere2, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject);
}

{
  auto right3Transform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    TV translation = TV(0, 0, .05 * time);
    TV translation_velocity(0, 0, .05);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> sphere3(TV(2.058, 1.83, 2.047), sphere_radius);
  AnalyticCollisionObject<T, dim> rightObject(
      right3Transform, sphere3, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject);
}

// init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx,
// AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for
// SPGrid.