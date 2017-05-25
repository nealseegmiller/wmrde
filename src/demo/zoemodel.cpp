#include <wmrde/demo/zoemodel.h>
#include <wmrde/wmrde_ros_interface.h> //for makeColorRGBA

namespace vm = visualization_msgs; //abbreviate

//TODO, move this to wmrde_ros_interface?
vm::Marker makeMeshMarker(
    const wmrde::WmrModel& mdl,
    const int frame_idx,
    const std::string& resource,
    const std_msgs::ColorRGBA color,
    const double scale)
{
  vm::Marker marker;
  marker.type = vm::Marker::MESH_RESOURCE;
  marker.id = frame_idx;
  marker.header.frame_id = mdl.getFrame(frame_idx).name;
  marker.mesh_resource = resource;
  marker.color = color;
  marker.scale.x = scale;
  marker.scale.y = scale;
  marker.scale.z = scale;
  return marker;
}

namespace wmrde
{

void makeZoeModel(WmrModel& mdl)
{
  //Zoe rover model

  //dimensions for kinematics (meters)
  const Real length = 1.91; //distance from front to rear axle
  const Real width = 1.64; //distance from left to right wheel
  const Real axle_drop = 0.119; //vertical offset between axle and wheel joints
  const Real wheel_rad = 0.325; //wheel radius

  const int num_wheels = 4; //number of wheels

  //masses (kg)
  const Real total_mass = 198.0; //kg
  const Real mass_body = 0.8*total_mass;
  const Real mass_axle = 0.1*total_mass;
  const Real mass_wheel = 0.1*total_mass/Real(num_wheels);

  //compute mass properties
  Vec3 cm_body = (Vec3() << 0, 0, 0.25).finished(); //body center of mass

  //approximate moments of inertia
  Mat3 I_body = calcMomentOfInertiaBox(mass_body, length, 1.0, 0.5);
  Mat3 I_axle = calcMomentOfInertiaBox(mass_axle, 0.1, width, 0.1);
  Mat3 I_wheel = calcMomentOfInertiaCylinder(mass_wheel, wheel_rad, 0.07, 1);

  MassProperties mp_body(mass_body, cm_body, I_body);
  MassProperties mp_axle(mass_axle, Vec3::Zero(), I_axle);
  MassProperties mp_wheel(mass_wheel, Vec3::Zero(), I_wheel);

  mdl = WmrModel();
  //add the body frame
  {
    mdl.addBodyFrame("body");
    mdl.setFrameMass(0, mp_body);
  }

  //add the front axle steer frame
  {
    //lock front axle roll
    HTransform HT_parent_jd0(Mat3::Identity(), (Vec3() << length/2.0, 0, 0).finished());
    Frame frame("front_steer", Frame::REV_Z, mdl.frameNameToIndex("body"), HT_parent_jd0, false, true, nullptr, mp_axle);
    mdl.addFrame(frame);
  }

  //add the rear axle roll frame
  {
    HTransform HT_parent_jd0(Mat3::Identity(), (Vec3() << -length/2.0, 0, 0).finished());
    Frame frame("rear_roll", Frame::REV_X, mdl.frameNameToIndex("body"), HT_parent_jd0, false, false, nullptr);
    mdl.addFrame(frame); //no mass
  }

  //add the rear axle steer frame
  {
    Frame frame("rear_steer", Frame::REV_Z, mdl.frameNameToIndex("rear_roll"), HTransform::Identity(), false, true, nullptr, mp_axle);
    mdl.addFrame(frame);
  }

  //add the wheel frames
  {
    //front wheels
    Vec3 t_left = (Vec3() << 0, width/2.0, -axle_drop).finished();
    Vec3 t_right = (Vec3() << 0, -width/2.0, -axle_drop).finished();

    HTransform HT_parent_jd0(Mat3::Identity(), t_left);
    Frame frame("front_left_wheel", Frame::REV_Y, mdl.frameNameToIndex("front_steer"), HT_parent_jd0, true, true, new WheelGeom(), mp_wheel);
    frame.wheel_geom->radius = wheel_rad;
    mdl.addFrame(frame);

    frame.name = "front_right_wheel";
    frame.HT_parent_jd0.t = t_right;
    mdl.addFrame(frame);

    //rear wheels
    frame.name = "rear_left_wheel";
    frame.parent_idx = mdl.frameNameToIndex("rear_steer");
    frame.HT_parent_jd0.t = t_left;
    mdl.addFrame(frame);

    frame.name = "rear_right_wheel";
    frame.HT_parent_jd0.t = t_right;
    mdl.addFrame(frame);
  }
}

void makeZoeModelMarkers(
    const WmrModel& mdl, //model set using makeZoeModel
    vm::MarkerArray& markers)
{
  markers = vm::MarkerArray();

  std::string prefix = "package://wmrde/CAD/Zoe/"; //url to mesh resources
  //colors
  std_msgs::ColorRGBA body_color = makeColorRGBA(0.5, 0, 0.5);
  std_msgs::ColorRGBA axle_color = makeColorRGBA(0.5, 0.5, 0.5);
  std_msgs::ColorRGBA wheel_color = makeColorRGBA(0.25, 0.25, 0.25);
  double scale = 1.0; //.stl files in meters

  markers.markers.push_back(makeMeshMarker(
      mdl, mdl.frameNameToIndex("body"), prefix + "ZoeBody.stl", makeColorRGBA(0.5, 0, 0.5), scale));
  markers.markers.push_back(makeMeshMarker(
      mdl, mdl.frameNameToIndex("front_steer"), prefix + "ZoeFrontAxle.stl", axle_color, scale));
  markers.markers.push_back(makeMeshMarker(
      mdl, mdl.frameNameToIndex("rear_steer"), prefix + "ZoeRearAxle.stl", axle_color, scale));
  for (int frame_idx : mdl.wheelFrameIndices())
  {
    markers.markers.push_back(makeMeshMarker(
        mdl, frame_idx, prefix + "ZoeWheel.stl", wheel_color, scale));
  }
}

} //namespace

/*
void zoe(WmrModel& mdl, Real state[], Real qvel[]) {
	//Zoe rover model

	const bool FIX_FRONT_AXLE_ROLL = true;
	
	//dimensions
	const Real L = 1.91; //front to rear
	const Real B = 1.64; //axle width
	const Real D = .119; //vertical offset between wheel centers and axle joints
	const Real rad = .325; //wheel radius

	const int nw=4; //number of wheels

	//masses
	const Real TotalMass=198; //kg
	const Real Mb = .9*TotalMass; //mass of body
	const Real Mw = .1*TotalMass/Real(nw); //mass of wheel

	const Real Hb = .5; //height of body
	const Real Wb = 1.0; //width of body
	const Real Ww = .07; //width of wheel
	Vec3 cmb = {0,0,Hb/2}; //body center of mass (in body coords)
	Vec3 cmw = {0,0,0}; //wheel center of mass

	mdl.addBodyFrame("Body");

	HomogeneousTransform HT;
	VecEuler euler = {0,0,0};
	VecOrient orient;
#if WMRSIM_USE_QUATERNION
	eulerToQuat(euler,orient);
#else
	copyEuler(euler,orient);
#endif
	Vec3 pos;
	
	//front axle
	setVec3(L/2,0,0,pos);
	poseToHT(orient,pos,HT);
	if (FIX_FRONT_AXLE_ROLL) {
		mdl.addJointFrame("SteerFront","Body","RZ",false,HT);
	} else {
		mdl.addJointFrame("RollFront","Body","RX",false,HT);
		setVec3(0,0,0,pos);
		poseToHT(orient,pos,HT);
		mdl.addJointFrame("SteerFront","RollFront","RZ",false,HT);
	}

	setVec3(0,B/2,-D,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("WheelFrontLeft","SteerFront",true,HT,rad); //1
	setVec3(0,-B/2,-D,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("WheelFrontRight","SteerFront",true,HT,rad); //2

	//rear axle
	setVec3(-L/2,0,0,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("RollRear","Body","RX",false,HT);
	setVec3(0,0,0,pos);
	poseToHT(orient,pos,HT);
	mdl.addJointFrame("SteerRear","RollRear","RZ",false,HT);
	
	setVec3(0,B/2,-D,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("WheelRearLeft","SteerRear",true,HT,rad); //3
	setVec3(0,-B/2,-D,pos);
	poseToHT(orient,pos,HT);
	mdl.addWheelFrame("WheelRearRight","SteerRear",true,HT,rad); //4


	//set mass properties
	int nf = mdl.get_nf();
	const int* wheelframeinds = mdl.get_wheelframeinds();

	Mat3 I; //Inertia
	inertiaBox(Mb,L,Wb,Hb,I);
	mdl.setFrameMass(0,Mb,cmb,I);

	inertiaCylinder(Mw,rad,Ww,1,I);
	for (int i=0; i<nw; i++)
		mdl.setFrameMass(wheelframeinds[i],Mw,cmw,I);


	//FOR KINEMATIC SIM
	mdl.min_npic = nw;

	//contact height error
	setVec(nw,-.02,mdl.dz_target);
	setVec(nw,.1,mdl.tc_z);
	if (!FIX_FRONT_AXLE_ROLL)
		mdl.tc_j[0] = .1;

	//FOR DYNAMIC SIM
	mdl.wheelGroundContactModel = uniformWgc;

	Real Kp = TotalMass*mdl.grav/(nw*-mdl.dz_target[0]);
	setWgcParams(Kp,mdl.wgc_p);

	mdl.actuatorModel = zoeAct;
	mdl.act_p[0] = 2e3; //Kp
	mdl.act_p[1] = 0.0; //Ki
	mdl.act_p[2] = REALMAX; //max

	//ODE contact model parameters
	Real erp,cfm;
	KpKdToErpCfm(Kp, Kp/20, .04, erp, cfm);
	Real fds = 1.0/(Kp*1e-1);

	setVec(nw,erp,mdl.erp_z);
	setVec(nw,cfm,mdl.cfm_z);
	setVec(nw,fds,mdl.fds_x);
	setVec(nw,fds,mdl.fds_y);
	if (!FIX_FRONT_AXLE_ROLL) {
		mdl.erp_j[0] = .2;
		mdl.cfm_j[0] = 1e-6;
	}

	//FOR BOTH
	//set function pointers
	mdl.controller = zoeController;

	if (!FIX_FRONT_AXLE_ROLL) {
		mdl.holonomicJointConstraints = zoeConstraints;
		mdl.set_njc(1);
	}

	//initialize the state vector
	setEuler(DEGTORAD(0),DEGTORAD(0),DEGTORAD(0),euler);
	setVec3(0.0, 0.0, rad + D-.01, pos);

#if WMRSIM_USE_QUATERNION
	eulerToQuat(euler,orient);
#else
	copyEuler(euler,orient);
#endif
	int ns = NUMSTATE(nf); //number of elements in state
	setVec(ns,0.0,state);
	copyOrient(orient,state+SI_ORIENT);
	copyVec3(pos,state+SI_POS);

	//keep steer angles fixed when initializing terrain contact
	int steerfront_fi, steerrear_fi;
	steerfront_fi = mdl.nameToInd("SteerFront");
	steerrear_fi = mdl.nameToInd("SteerRear");

	mdl.set_isfixed(steerfront_fi, true);
	mdl.set_isfixed(steerrear_fi, true);

	//DEBUGGING, INIT STEER ANGLES
	if (0) {
		Real sa = DEGTORAD(15);
		state[TOSTATEI(steerfront_fi)] = sa; 
		state[TOSTATEI(steerrear_fi)] = -sa; 
	}


	if (qvel != 0) { //not null
		//initialize qvel
		//setVec(mdl.get_na(), 0.0, qvel); //to zeros

		//to cmd
		Real u[WmrModel::MAXNA];
		mdl.controller(mdl, 0.0, state, u, qvel);
	}

}


void zoeController(const WmrModel& mdl, const Real time, const Real state[], //inputs
					Real u[], Real qvel_cmd[]) { //outputs


	static int steerfront_fi, steerrear_fi; //steer frame indices, in WmrModel
	static bool init = false; //flag, static variables initialized
	const int nw = 4;
	const int front = 0;
	const int rear = 1;

	if (!init) {
		steerfront_fi = mdl.nameToInd("SteerFront");
		steerrear_fi = mdl.nameToInd("SteerRear");
	}

	Real speed,turnrad; //commanded speed, turn radius

	speed = .5;
	turnrad = 1000;

	
	//if (time < .5) {
	//	speed = .5;
	//	turnrad = 1000;
	//} else if (time >= .5 && time < 1) {
	//	speed = 1.0;
	//	turnrad = 1000;
	//} else if (time >= 1 && time < 5) {
	//	speed = 1.0;
	//	turnrad = 5.0;
	//} else if (time >= 5) {
	//	speed = 1.0;
	//	turnrad = -5.0;
	//}
	

	//measured steer angles
	Real steer[2];
	steer[front] = state[TOSTATEI(steerfront_fi)];
	steer[rear] = state[TOSTATEI(steerrear_fi)];

	//allocate outputs
	Real steer_cmd[2];
	Real yawrate = 0;
	Real steer_rate[2];

	zoeArcController( speed, turnrad, steer, //inputs
		steer_cmd, yawrate, u, steer_rate); //outputs

	
	if (qvel_cmd != 0) { //not null

		//get from WmrModel
		int nv = NUMQVEL(mdl.get_nf());
		const int* wheelframeinds = mdl.get_wheelframeinds();

		setVec(nv, 0.0, qvel_cmd);
		qvel_cmd[VI_ANG + 2] = yawrate;
		qvel_cmd[VI_LIN + 0] = speed;

		//wheel velocities
		for (int i=0; i<nw; i++)
			qvel_cmd[TOQVELI(wheelframeinds[i])] = u[i];
		
		//steer rates
		qvel_cmd[TOQVELI(steerfront_fi)] = steer_rate[front];
		qvel_cmd[TOQVELI(steerrear_fi)] = steer_rate[rear];
	}
}

void zoeArcController( const Real speed, const Real turnrad, const Real* steer, //inputs
	Real* steer_cmd, Real yawrate, Real* wheelvel, Real* steer_rate) { //outputs
	//turn arc commands into wheel velocity commands

	const int front = 0;
	const int rear = 1;

	Real L,B,rad; //dimensions
	L = 1.91;
	B = 1.64;
	rad = .325;

	//options
	const Real min_rad = 3.0;	//min turn radius
	const Real max_dvel = 1.5;	//max difference btwn left & right wheel vel due to feedback
	const Real Kp = 1.0;		//proportional gain

	//compute commanded steer angles and yaw rate

	

	if (fabs(turnrad) >= 1000) {
		//drive straight
		steer_cmd[0] = 0;
		steer_cmd[1] = 0;
		yawrate = 0;
	} else {
		//enforce radius limit
		Real turnrad_new = turnrad;
		if (fabs(turnrad) < min_rad) 
			turnrad_new = REALSIGN(turnrad)*min_rad;

		steer_cmd[front] = atan(L/2/turnrad_new);
		steer_cmd[rear] = -steer_cmd[front];

		yawrate = speed/turnrad_new;
	}

	//compute commanded wheel velocities

	//feedforward
	const int nw=4;
	Real tmp;

	tmp = 1.0/cos(steer_cmd[front]);
	wheelvel[0] = tmp*speed - B/2*yawrate; //front left
	wheelvel[1] = tmp*speed + B/2*yawrate; //front right
	tmp = 1.0/cos(steer_cmd[rear]);
	wheelvel[2] = tmp*speed - B/2*yawrate; //rear left
	wheelvel[3] = tmp*speed + B/2*yawrate; //rear right

	//feedback
	Real dvel_front,dvel_rear;
	//front wheels
	dvel_front = Kp*(steer_cmd[front] - steer[front]);
	if ( fabs(dvel_front) > max_dvel ) 
		dvel_front = REALSIGN(dvel_front)*max_dvel;
	wheelvel[0] -= dvel_front; //left
	wheelvel[1] += dvel_front; //right
	//rear wheels
	dvel_rear = Kp*(steer_cmd[rear] - steer[rear]);
	if ( fabs(dvel_rear) > max_dvel ) 
		dvel_rear = REALSIGN(dvel_rear)*max_dvel;
	wheelvel[2] -= dvel_rear;
	wheelvel[3] += dvel_rear;

	//convert from m/s to rad/s
	for (int i=0; i<nw; i++) 
		wheelvel[i] = wheelvel[i]/rad;

	//steer angle rates in rad/s
	steer_rate[front] = 2*dvel_front/B;
	steer_rate[rear] = 2*dvel_rear/B;
}


void zoeConstraints( const WmrModel& mdl, const Real jd[], const Real jr[], //inputs
					 Real c[], Real Jc[], Real f[], Real df_djd[], Real df_djr[]) { //outputs
	static int rollfront_ji, rollrear_ji; //roll angle indices (in jd)
	static bool init = false; //flag, initialized

	if (!init) {
		rollfront_ji = mdl.nameToInd("RollFront")-1;
		rollrear_ji = mdl.nameToInd("RollRear")-1;
	}

	c[0] = jd[rollfront_ji] + jd[rollrear_ji]; //angles should be equal magnitude, opposite sign
	
	
	int nj = mdl.get_nf()-1; //number of joints
	

	//init to zeros
	setVec(nj,0.0,Jc);

	Jc[rollfront_ji] = 1.0;
	Jc[rollrear_ji] = 1.0;

	if ( f != 0 ) { //if not null
		//for dynamic sim
		Real Kp = 22500;
		Real Kd = Kp/20;
		Real cd = jr[rollfront_ji] + jr[rollrear_ji];

		f[0] = -(Kp*c[0] + Kd*cd);
		if ( df_djd != 0 ) {
			copyVec(nj,Jc,df_djd);
			mulcVec(nj,-Kp,df_djd);

			copyVec(nj,Jc,df_djr);
			mulcVec(nj,-Kd,df_djr);
		}
	}
}
*/
