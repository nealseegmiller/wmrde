#ifndef _WMRDE_WMRMODEL_H_
#define _WMRDE_WMRMODEL_H_

#include <memory>

#include <wmrde/algebra/spatial.h>
#include <wmrde/contactgeom.h>

namespace wmrde
{

//TODO, move this?
class MassProperties
{
/*!
 * \class MassProperties
 * A class for storing mass properties that precomputes the conversion
 * to/from a spatial inertia matrix.
 */
 public:
  //constructors
  MassProperties() {}
  MassProperties(
      const Real scalar_mass,
      const Vec3& center_of_mass,
      const Mat3& moment_of_inertia)
 {
    m = scalar_mass;
    cm = center_of_mass;
    I = moment_of_inertia;
    Is = toSpatialInertia(m,cm,I);
 }

  MassProperties(const Mat6& spatial_inertia)
  {
    Is = spatial_inertia;
    fromSpatialInertia(Is,m,cm,I);
  }

  //get methods
  const Real getScalarMass() const { return m; }
  const Vec3& getCenterOfMass() const { return cm; }
  const Mat3& getMomentOfInertia() const { return I; }
  const Mat6& getSpatialInertia() const { return Is; }

 private:

  Real m;   //scalar mass
  Vec3 cm;  //center of mass
  Mat3 I;   //moment of inertia about the center of mass
  Mat6 Is; //spatial inertia
};

Mat3 calcMomentOfInertiaBox(
    const Real mass,
    const Real lx,
    const Real ly,
    const Real lz);

Mat3 calcMomentOfInertiaCylinder(
    const Real mass,
    const Real radius,
    const Real height,
    const int axis); //x=0, y=1, z=2

//Frame is a sub class of WmrModel
class Frame
{
 public:
  enum DofTypes {
    FREE = -1, //6 Dof
    REV_X = 0, //revolute about x
    REV_Y = 1, //revolute about y
    REV_Z = 2, //revolute about z
    LIN_X = 3, //linear along x
    LIN_Y = 4, //linear along y
    LIN_Z = 5, //linear along z
  };

  Frame() :
    name("UNINITIALIZED"),
    dof_type(-1),
    parent_idx(-1),
    is_actuated(false),
    is_fixed(false)
  {
//    HT_parent_jd0.setIdentity();
  }

  //set every member in constructor
  Frame(const std::string& frame_name,
        const int degree_of_freedom_type,
        const int parent_frame_index,
        const HTransform& HT_parent_joint_disp0,
        const bool frame_is_actuated,
        const bool frame_is_fixed,
        WheelGeom* wheel_geom_ptr = nullptr, //Frame will take ownership of the object pointed to
        const MassProperties& mass_properties = MassProperties())
  :
    name(frame_name),
    dof_type(degree_of_freedom_type),
    parent_idx(parent_frame_index),
    HT_parent_jd0(HT_parent_joint_disp0),
    is_actuated(frame_is_actuated),
    is_fixed(frame_is_fixed),
    wheel_geom(wheel_geom_ptr),
    mass(mass_properties)
  {}

  //default copy constructor is deleted because of unique_ptr member
  //TODO, a better way?
  Frame deepcopy() const
  {
    return Frame(name, dof_type, parent_idx, HT_parent_jd0, is_actuated, is_fixed,
        wheel_geom ? new WheelGeom(*wheel_geom) : nullptr, //deepcopy, handle nullptr
        mass);
  }

  std::string name;
  int dof_type; //degree of freedom type, set to one of the DofTypes values
  int parent_idx; //index of parent frame
  HTransform HT_parent_jd0; //Homogeneous Transform of frame wrt parent if joint displacement=0
  bool is_actuated;
  bool is_fixed; //if true, keep fixed when initializing terrain contact
  std::unique_ptr<WheelGeom> wheel_geom; //nullptr for non-wheel frames
  MassProperties mass; //Only necessary for dynamic sim.
};


class WmrModel
{
/*!
 * \class WmrModel
 * A class to represent the kinematics and mass properties of a wheels mobile robot
 */
 public:
	//constructor
	WmrModel() :
	  num_wheel_frames_(0),
	  num_actuated_frames_(0)
  {}
	
	/*!
	 * Adds named body frame to the WmrModel. Must call this before
	 * calling addFrame().
	 * \param name the body frame name
	 */
	bool addBodyFrame(const std::string& name);

	/*!
	 * Adds a frame to the WmrModel. must call addBodyFrame() first
	 * Checks the validity of frame member variables:
	 * *name must not match any existing frame.
	 * *dof_type must be one of the DofType enum values.
	 * *parent_idx must be valid index of existing frame.
	 * \frame the frame to add
	 * \return true if frame added successfully
	 */
	bool addFrame(const Frame& frame);
	
	void setFrameMass(const int frame_idx, const MassProperties& mass) { frames_[frame_idx].mass = mass; }

	//get methods
  /*!
   * returns index of frame with the specified name
   * returns -1 if name doesn't match any frame
   */
  int frameNameToIndex(const std::string& name) const ;

	/*!
	 * Return const reference to a frame.
	 */
	const Frame& getFrame(int idx) const { return frames_[idx]; }
	
	int numFrames() const { return frames_.size(); }
	int numJoints() const { return std::max(0, numFrames()-1); } //body frame has no joint
	int numWheelFrames() const { return num_wheel_frames_; }
	int numActuatedFrames() const { return num_actuated_frames_; }
	const std::vector<int>& getWheelFrameIndices() const { return wheel_frame_inds_; }
	const std::vector<int>& getActuatedFrameIndices() const { return actuated_frame_inds_; }
	
 private:
	std::vector<Frame> frames_;
	
	//these values are updated automatically in calls to addFrame()
	//This makes it faster to iterate over frame subsets
	int num_wheel_frames_;
	int num_actuated_frames_;
	std::vector<int> wheel_frame_inds_;
	std::vector<int> actuated_frame_inds_;
};

} //namespace

/*
  //functions to construct model
  void addBodyFrame(const std::string name);
  void addJointFrame(const std::string name, const std::string parent_name, const std::string dof_string, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0);
  void addWheelFrame(const std::string name, const std::string parent_name, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0,
    const Real radius);
  void addSprocketFrame(const std::string name, const std::string parent_name, const bool isactuated, const HomogeneousTransform HT_parent_jointdisp0,
    const Real radius, const Real radius2, const Real len);

  void setFrameMass(const int fi, const Real mass, const Vec3 center_of_mass, const Mat3 moment_of_inertia);
  void set_isfixed(const int fi, const bool val) { frames[fi].isfixed = val; }
  void set_njc( const int num_joint_constraints ) { njc = num_joint_constraints; }

  int nameToInd(const std::string Name) const; //returns index of frame by name

  //get methods
  //const correctness, const before { means read only access to WmrModel object
  const int get_nf() const { return nf; }


  const Frame* get_frames() const { return frames; }
  const Frame* get_frameptr(const int fi) const { return frames+fi; }
  //const Frame get_frame(const int fi) const { return frames[fi]; } //don't use this, slow!

  const int get_njc() const { return njc; }

  //
  const int get_nw() const { return nw; }
  const int* get_wheelframeinds() const { return wheelframeinds; }
  const int get_nt() const { return nt; }
  const int* get_sprocketframeinds() const { return sprocketframeinds; }
  const int get_na() const { return na; }
  const int* get_actframeinds() const { return actframeinds; }

  //TODO, make these private?
  //FOR KINEMATIC MODEL
  int min_npic;     //minimum number of points in contact
  Real dz_target[MAXNW];  //target Delta z (contact height error)
  //time constants for Baumgarte's stabilization method:
  Real tc_z[MAXNW];   //wheel-ground contact constraints, z dir
  Real tc_j[MAXNJC];    //holonomic joint constraints

  //FOR DYNAMIC MODEL
  Real grav;        //scalar acceleration of gravity
  bool use_constraints;


  //user must set these function pointers

  void (*controller) (const WmrModel& mdl, const Real time, const Real state[], //inputs
    Real u[], Real qvel_cmd[]); //outputs

  void (*wheelGroundContactModel) ( const int wheelno, const Real params[], const Vec3 vc, const Real Rw, const Real dz, //inputs
    Vec3 fw, Real J[]); //outputs
  Real wgc_p[MAXNPAR]; //wheel-ground contact model params
  //wheelno:  wheel number (in a size nw array)
  //J:    3x5 Jacobian of fw wrt [vx vy vz Rw dz]


  //one function for all actuators, may be coupled
  void (*actuatorModel) ( const Real params[], const Real ucmd[], const Real u[], const Real interr[], //inputs
    Real f[], Real err[], Real* dfdu); //outputs
  Real act_p[MAXNPAR]; //actuator model params
  //INPUTS
  //params: size np
  //ucmd:   size na, commanded joint velocities
  //u:    size na, actual joint velocities
  //interr: size na, integrated error
  //OUTPUTS
  //f:    size na, forces
  //err:    size na, error
  //dfdu:   na x na Jacobian, df/du

  //must also set njc > 0
  void (*holonomicJointConstraints) ( const WmrModel& mdl, const Real jd[], const Real jr[], //inputs
    Real c[], Real Jc[], Real f[], Real df_djd[], Real df_djr[]); //outputs
  //nj = number of joints
  //nc = number of joint constraints
  //INPUTS
  //jd:   size nj, joint displacement
  //jr:   size nj, joint rates (d/dt jd)
  //OUTPUTS
  //c:    size nc, evaluation of joint constraint function, zeros if satisfied
  //Jc:   nc x nj, Jacobian of c wrt jd
  //f:    size nc, force
  //df_djd: nc x nj, Jacobian of f wrt jd
  //df_djr: nc x nj, Jacobian of f wrt jr



  //the following are only used if no wheel-ground contact model is provided

  //for wheel-ground contact constraints
  Real erp_z[MAXNW];    //error reduction parameter, z dir (0-1)
  Real cfm_z[MAXNW];    //constraint force mixing, z dir (>0)
  Real fds_x[MAXNW];    //force dependent slip, x dir
  Real fds_y[MAXNW];    //force dependent slip, y dir

  //for holonomic joint constraints
  Real erp_j[MAXNJC];   //error reduction parameter
  Real cfm_j[MAXNJC];   //constraint force mixing
 */

#endif //_WMRDE_WMRMODEL_H_
