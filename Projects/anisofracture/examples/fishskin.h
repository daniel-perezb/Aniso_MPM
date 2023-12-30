
std::string helper_output;
helper_output = "output/dual_fleece";
std::cout << "Attempt 1" << std::endl;
sim.output_dir.path = helper_output;
sim.end_frame = 160;
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
if (1) {
  T Youngs = 10000;
  T nu = 0.4;
  T rho = 500;
  T helper_isotropic = false;
  TV helper_fiber = TV(0.7071, 0, 0.7071); // 45 degree angle fibres
  T helper_alpha = -1;

  std::string filename = "TetMesh/split/sample1_bottom.mesh";
  MpmParticleHandleBase<T, dim> particles_handle =
      init_helper.sampleFromTetWildFile(filename, rho);
  T total_volume = particles_handle.total_volume;
  T particle_count = particles_handle.particle_range.length();
  T per_particle_volume = total_volume / particle_count;
  sim.dx = std::pow(particle_per_cell * per_particle_volume, (T)1 / (T)3);
  if (1) {
    StdVector<TV> samples;
    StdVector<Vector<int, 4>> indices;
    std::string absolute_path = DataDir().absolutePath(filename);
    readTetMeshTetWild(absolute_path, samples, indices);
    sim.output_dir.createPath();
    std::string vtk_path = DataDir().path + "/../Projects/anisofracture/" +
                           sim.output_dir.path + "/tet1.vtk";
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
  T percentage = 0.15;
  T l0 = 0.5 * sim.dx;
  T eta = 0.01;
  T zeta = 1;
  bool allow_damage = true;
  T residual_stress = 0.001;
  // model.scaleFiberStiffness(0, 2);
  particles_handle.addFBasedMpmForceWithAnisotropicPhaseField(
      a_0, alphas, percentage, l0, model, eta, zeta, allow_damage,
      residual_stress);
}

// ****************************************************************************
// Bottom
// ****************************************************************************
if (1) {
  T Youngs = 50000;
  T nu = 0.4;
  T rho = 500;
  T helper_isotropic = true;
  TV helper_fiber = TV(1, 1, 1); // No fibre direction (Isotropic model)
  T helper_alpha = 0;

  std::string filename = "TetMesh/split/sample1_top.mesh";
  MpmParticleHandleBase<T, dim> particles_handle =
      init_helper.sampleFromTetWildFile(filename, rho);
  T total_volume = particles_handle.total_volume;
  T particle_count = particles_handle.particle_range.length();
  T per_particle_volume = total_volume / particle_count;
  sim.dx = std::pow(particle_per_cell * per_particle_volume, (T)1 / (T)3);
  if (1) {
    StdVector<TV> samples;
    StdVector<Vector<int, 4>> indices;
    std::string absolute_path = DataDir().absolutePath(filename);
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
  T residual_stress = 0.005;
  // model.scaleFiberStiffness(0, 2);
  particles_handle.addFBasedMpmForceWithAnisotropicPhaseField(
      a_0, alphas, percentage, l0, model, eta, zeta, allow_damage,
      residual_stress);
}

// ****************************************************************************
// Collision objects
// ****************************************************************************
double sphere_radius = 0.01;

{
  auto leftTransform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    T t = time;
    TV translation = TV(-0.05 * time, 0, 0);
    TV translation_velocity(-0.05, 0, 0);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> leftLS(TV(2.02, 1.83, 1.9), sphere_radius);
  AnalyticCollisionObject<T, dim> leftObject(
      leftTransform, leftLS, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(leftObject);
}

{
  auto rightTransform = [](T time, AnalyticCollisionObject<T, dim> &object) {
    TV translation = TV(0.05 * time, 0, 0);
    TV translation_velocity(.05, 0, 0);
    object.setTranslation(translation, translation_velocity);
  };
  Sphere<T, dim> sphere2(TV(2.05, 1.83, 1.9), sphere_radius);
  AnalyticCollisionObject<T, dim> rightObject(
      rightTransform, sphere2, AnalyticCollisionObject<T, dim>::STICKY);
  init_helper.addAnalyticCollisionObject(rightObject);
}
