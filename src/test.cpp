#include <wmrde/test.h>

inline double tosec(timeval tim)
{
  return tim.tv_sec + (tim.tv_usec/1000000.0);
}

void test_common() {

	if (1) {
		//test angle macros in common.h
		Real a=181;
		Real b=-181;
		Real c=b-360;
		std::cout << "DIFFDEG(" << a << "," << b << ") = " << DIFFDEG(a,b) << std::endl;
		std::cout << "DIFFDEG(" << a << "," << c << ") = " << DIFFDEG(a,c) << std::endl;
		std::cout << "WRAPDEG(" << c << "," << 180 << ") = " << WRAPDEG(c,180) << std::endl;
		std::cout << std::endl;
	}

	if (1) {
		//test sortIndex()

		const int SIZE = 5;
		Real val[SIZE] = {-1.5, 0, 2.5, -2.5, 1.5};
		int idx[SIZE];
		Real val_sort[SIZE];

		sortIndex(SIZE, val, idx, val_sort);

		std::cout << "val=\n"; printMatReal(1,SIZE,val,-1,-1);
		std::cout << "idx=\n"; printMatInt(1,SIZE,idx,-1);
		std::cout << "val_sort=\n"; printMatReal(1,SIZE,val_sort,-1,-1);
	}

}

void test_linalg3() {
	//compare Eigen vs. linalg3.h for 3x3 matrix addition, multiplication

	//for timing
	int n = (int) 1e9;
	clock_t t;

	Mat3 A,B,C;

	setMat3(1,2,3,4,5,6,7,8,9,A);
	addcMat3(A,1,B);

	//computation time
	t=clock();
	for (int i=0; i<n; i++) {
		multMatMat3(A,B,C);
	}
	t=clock()-t;

	//must print or compiler will optimize out loop
	std::cout << "(using linalg3) C=\n";
	printMat3(C,-1,-1);
	std::cout << "iterations: " << (Real) n << std::endl;
	std::cout << "clock (ms): " << t << std::endl;
	std::cout << std::endl;


	if (1) {
		//test Eigen
		Eigen::Matrix<Real,3,3> A_;
		Eigen::Matrix<Real,3,3> B_;
		Eigen::Matrix<Real,3,3> C_;

		//.data() returns a pointer to the data array of the matrix
		copyVec3ToArray(3, (Vec3*) A, A_.data());
		copyVec3ToArray(3, (Vec3*) B, B_.data());

		//std::cout << "A=\n" << A_ << std::endl;
		//std::cout << "B=\n" << B_ << std::endl;
	
		t=clock();
		for (int i=0; i<n; i++) { 
			//C_ = A_+B_;
			C_ = A_*B_;
		}
		t=clock()-t;

		std::cout << "(using Eigen) C=\n" << C_ << std::endl;
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
		std::cout << std::endl;
	}
	
}


void test_transform() {

	VecEuler euler;
	Vec3 translation;
	HomogeneousTransform HT,HT2;

	setEuler(DEGTORAD(10),DEGTORAD(20),DEGTORAD(30),euler);
	setVec3(1,2,3,translation);
	poseToHT(euler,translation,HT);

	std::cout << "HT=\n"; printHT(HT,-1,-1);
	

	//test rotation matrix
	Mat3 Rx,Ry,Rz,Rot,Tmp;

	Rotx(euler[0],Rx);
	Roty(euler[1],Ry);
	Rotz(euler[2],Rz);
	multMatMat3(Ry,Rx,Tmp);
	multMatMat3(Rz,Tmp,Rot);
	std::cout << "Rotz(yaw)*Roty(pit)*Rotx(rol)=\n"; printMat3(Rot,-1,-1);

	//test invert transform
	invertHT(HT,HT2);
	std::cout << "inv(HT)=\n"; printHT(HT2,-1,-1);

	//test compose HT
	if (1) {
		//time it

		int n = (int) 1e8;
		clock_t t;

		//compose Homogeneous Transforms
		t=clock();
		for (int i=0; i<n; i++) { 
			composeHT(HT,HT,HT2);
		}
		t=clock()-t;

		std::cout << "HT*HT=\n"; printHT(HT2,-1,-1);
	
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl << std::endl;
	} else {
		composeHT(HT,HT,HT2);
		std::cout << "HT*HT=\n"; printHT(HT2,-1,-1);
	}

	composeInvHT(HT,HT,HT2);
	std::cout << "inv(HT)*HT=\n"; printHT(HT2,-1,-1);


	//test transform point
	Vec3 p,q;
	setVec3(2,3,4,p);

	std::cout << "p=\n"; printVec3(p,-1,-1);
	applyHT(HT,p,q);
	std::cout << "q=HT*p=\n";	printVec3(q,-1,-1);

	applyInvHT(HT,p,q);
	std::cout << "q=inv(HT)*p=\n"; printVec3(q,-1,-1);

	
	//test converting between angular velocity and Euler rates
	Vec3 vel={1,2,3};
	VecEuler eulerrate;
	MatEuler T;
	std::cout << "angular velocity=\n";	printVec3(vel,-1,-1);

	velToEulerrate(euler,vel,eulerrate,T);
	std::cout << "eulerrate=\n"; printEuler(eulerrate,-1,-1);
	std::cout << "T_vel_to_eulerrate=\n"; printMatReal(3,3,T,-1,-1);

	std::cout << "convert back:\n";
	eulerrateToVel(euler,eulerrate,vel,0);
	std::cout << "angular velocity=\n";	printVec3(vel,-1,-1);
	std::cout << std::endl;

	
	//test at singularity
	setEuler(0,M_PI/2,0,euler);
	setVec3(0,0,1,vel);
	velToEulerrate(euler,vel,eulerrate,T);
	std::cout << "test at singularity:\n";
	std::cout << "euler=\n"; printEuler(euler,-1,-1);
	std::cout << "angular velocity=\n"; printVec3(vel,-1,-1);
	std::cout << "eulerrate=\n"; printEuler(eulerrate,-1,-1);
	
}



void test_spatial() {
	//for timing
	int n;
	clock_t t;

	//construct Homogeneous Transform
	VecEuler euler;
	Vec3 translation;
	HomogeneousTransform HT;

	setEuler(DEGTORAD(10),DEGTORAD(20),DEGTORAD(30),euler);
	setVec3(1,2,3,translation);
	poseToHT(euler,translation,HT);


	std::cout << "HT=\n";
	printHT(HT,-1,-1);

	//convert to Plucker transform
	Mat6b P;

	HTToPlucker(HT,P);
	std::cout << "P=\n"; printMat6b(P,-1,-1);

	//inverse Plucker transform
	Mat6b P2;
	
	if (1) {
		//time it
		n= (int) 1e8;
		
		t=clock();
		for (int i=0; i<n; i++) {
			invHTToPlucker(HT,P2);
		}
		t=clock()-t;
		std::cout << "inv(P)=\n"; printMat6b(P2,-1,-1);
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
		std::cout << std::endl;
	} else {
		invHTToPlucker(HT,P2);
		std::cout << "inv(P)=\n"; printMat6b(P2,-1,-1);
	}

	Vec6b m,f,v,v2;
	setVec6b(1,2,3,4,5,6,m);
	setVec6b(2,3,4,5,6,7,f);

	std::cout << "m=\n"; printVec6b(m,-1,-1);
	std::cout << "f=\n"; printVec6b(f,-1,-1);

	//cross products
	crossVec6bMotion(m,f,v);
	std::cout << "crossMotion(m)*f=\n"; printVec6b(v,-1,-1);

	crossVec6bForce(m,f,v);
	std::cout << "crossForce(m)*f=\n"; printVec6b(v,-1,-1);
	
	//multiply vector by Plucker transform
	std::cout << "using spatial.h\n";
	if (1) {
		//time it
		n= (int) 1e8;
		
		t=clock();
		for (int i=0; i<n; i++) {
			multPluckerVec6b(P,m,v);
			multPluckerTVec6b(P,m,v2);
		}
		t=clock()-t;
		std::cout << "P*m=\n"; printVec6b(v,-1,-1);
		std::cout << "P'*m=\n"; printVec6b(v2,-1,-1);
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
		std::cout << std::endl;


	} else {
		multPluckerVec6b(P,m,v);
		std::cout << "P*m=\n"; printVec6b(v,-1,-1);

		multPluckerTVec6b(P,m,v2);
		std::cout << "P'*m=\n"; printVec6b(v2,-1,-1);
	}

	if (1) {
		//using Eigen
		Eigen::Matrix<Real,6,6> P_;
		Eigen::Matrix<Real,6,1> m_, v_, v2_;
		
		//copy to Eigen matrices
		copyMat6bToArray(P,P_.data());
		copyVec6bToArray(m,m_.data());

		t=clock();
		for (int i=0; i<n; i++) {
			v_ = P_*m_;
			v2_ = P_.transpose()*m_;
		}
		t=clock()-t;

		std::cout << "using Eigen\n";
		std::cout << "P*m=\n" << v_ << std::endl;
		std::cout << "P'*m=\n" << v2_ << std::endl;
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}

	//compute P'*I*P
	//copy tricycle body frame
	Real mass = 10;
	Real L = 1;
	Real W = .8;
	Real H = .25;
	Vec3 cm;
	Mat3 I;
	setVec3(L/2,0,H/2,cm);

	inertiaBox(mass,L,W,H,I);

	//using spatial.h
	Mat6b Is; //spatial inertia
	toSpatialInertia(mass,cm,I,Is);

	std::cout << "Is=\n";
	printMat6b(Is,-1,-1);

	Mat6b Is2;

	std::cout << "using spatial.h\n";
	if (1) {
		//time it
		n = (int) 1e7;

		t=clock();
	
		for (int i=0; i<n; i++) {
			//multPluckerTMat6bPlucker(P,Is,Is2); //deprecated
			multPluckerTInertiaPlucker(P,Is,Is2);
		}

		t=clock()-t;
		std::cout << "P'*Is*P=\n";	printMat6b(Is2,-1,-1);
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
		std::cout << std::endl;
	} else {
		//multPluckerTMat6bPlucker(P,Is,Is2);
		multPluckerTInertiaPlucker(P,Is,Is2);

		std::cout << "P'*Is*P=\n";	printMat6b(Is2,-1,-1);
	}
	
	if (1) {
		//using Eigen
		Eigen::Matrix<Real,6,6> P_, Is_, Is2_;
		
		//copy to Eigen matrices
		copyMat6bToArray(P,P_.data());
		copyMat6bToArray(Is,Is_.data());

		t=clock();
		for (int i=0; i<n; i++) {
			Is2_ = P_.transpose()*Is_*P_;
		}
		t=clock()-t;

		std::cout << "using Eigen\n";
		std::cout << "P'*Is*P=\n" << Is2_ << std::endl;
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
	
}



void test_matrix() {
	//matrix.h vs. Eigen

	//for timing
	int n = (int) 1e6;
	clock_t t;

	const int SIZE = 10;

	//matrix.h matrices
	Real A[SIZE*SIZE];
	Real B[SIZE*SIZE];
	Real B0[SIZE*SIZE]; //backup, to initialize Eigen
	Real b[SIZE];

	for (int i=0; i<SIZE*SIZE; i++) {
		A[i] = (Real) i+1;
		B[i] = 1;
		B0[i] = B[i];
	}

	for (int i=0; i<SIZE; i++) {
		b[i] = 1;
	}


	//std::cout << "A=\n"; printMatReal(SIZE,SIZE,A,-1,-1);
	//std::cout << "B=\n"; printMatReal(SIZE,SIZE,B,-1,-1);
	//std::cout << "b=\n"; printMatReal(SIZE,1,b,-1,-1);


	//matrix.h
	int ri = 2;
	int ci = 3;
	int brows = SIZE/2;
	int bcols = SIZE/2;
	Real val = 2.0;
	//Real m = 1.0 + 1e-6;
	Real m = 1.0;

	t=clock();

	for (int i=0; i<n; i++) {
		//setMat(SIZE,SIZE,val,B);
		//setMatRow(SIZE,SIZE,ri,val,B);
		//setMatCol(SIZE,ci,val,B);
		//setMatBlock(SIZE,ri,ci,brows,bcols,val,B);
		
		//mulcMat(SIZE,SIZE,m,B);
		//mulcMatRow(SIZE,SIZE,ri,m,B);
		//mulcMatCol(SIZE,ci,m,B);
		//mulcMatBlock(SIZE,ri,ci,brows,bcols,m,B);
		
		//copyMat(SIZE,SIZE,A,B);
		//copyMatRow(SIZE,SIZE,ri,A,SIZE,ri,B);
		//copyMatCol(SIZE,ci,A,ci,B);
		//copyMatBlock(SIZE,ri,ci,brows,bcols,A, SIZE,ri,ci,B);
		
		//copyTMat(SIZE,SIZE,A,B);
		
		//addmMat(SIZE,SIZE,A,m,B);
		//addmMatRow(SIZE,SIZE,ri,A,SIZE,ri,m,B);
		//addmMatCol(SIZE,ci,A,ci,m,B);
		//addmMatBlock(SIZE,ri,ci,brows,bcols,A, SIZE,ri,ci,m,B);

		//multMatVec(SIZE,SIZE,A,b,m,B);
		//multMatTVec(SIZE,SIZE,A,b,m,B);
		//multMatMat(SIZE,SIZE,A,SIZE,A,m,B);
		multMatTMat(SIZE,SIZE,A,SIZE,A,m,B);

	}
	t=clock()-t;

	std::cout << "(matrix.h) B=\n"; printMatReal(SIZE,SIZE,B,0,-1);
	std::cout << "iterations: " << (Real) n << std::endl;
	std::cout << "clock (ms): " << t << std::endl;
	std::cout << std::endl;


	
	if (1) {
		//Eigen

		//dynamic Eigen matrices
		Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> A_(SIZE,SIZE);
		Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> B_(SIZE,SIZE);
		Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> b_(SIZE,1);

		//fixed Eigen matrices
		//Eigen::Matrix<Real,SIZE,SIZE> A_;
		//Eigen::Matrix<Real,SIZE,SIZE> B_;
		//Eigen::Matrix<Real,SIZE,1> b_;

		memcpy(A_.data(), A, sizeof(Real)*SIZE*SIZE);
		memcpy(B_.data(), B0, sizeof(Real)*SIZE*SIZE);
		memcpy(b_.data(), b, sizeof(Real)*SIZE);

		//std::cout << "Eigen:\n";
		//std::cout << "A_=\n" << A_ << std::endl;
		//std::cout << "B_=\n" << B_ << std::endl;
		//std::cout << "b_=\n" << b_ << std::endl;
	
		t=clock();
		for (int i=0; i<n; i++) {
			//B_.setConstant(val);
			//B_.row(ri).setConstant(val);
			//B_.col(ci).setConstant(val);
			//B_.block(ri,ci,brows,bcols).setConstant(val);
		
			//B_ *= m;
			//B_.row(ri) *= m;
			//B_.col(ci) *= m;
			//B_.block(ri,ci,brows,bcols) *= m;
		
			//B_ = A_;
			//B_.row(ri) = A_.row(ri);
			//B_.col(ci) = A_.col(ci);
			//B_.block(ri,ci,brows,bcols) = A_.block(ri,ci,brows,bcols);
		
			//B_ = A_.transpose();

			//B_ += (A_ * m);
			//B_.row(ri) += (A_.row(ri) * m);
			//B_.col(ci) += (A_.col(ci) * m);
			//B_.block(ri,ci,brows,bcols) += (A_.block(ri,ci,brows,bcols) * m);
		
			//B_.col(0) = A_*b_;
			//B_.col(0) = A_.transpose()*b_;
			//B_ = A_*A_;
			B_ = A_.transpose()*A_;
		}
		t=clock()-t;

		std::cout << "(Eigen) B=\n" << B_ << std::endl;
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
}



void test_surface() {
	//vector of pointers to base surface class

	SurfaceVector surfs;

	//flat(surfs);
	//ramp(surfs);
	grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	const int MAXNP=1000;

	Vec3 pts[MAXNP];
	Real height[MAXNP];
	Real dz[MAXNP];
	int surfidx[MAXNP]; //surface index
	Vec3 normal[MAXNP];
	
	
	Real lo = 2.1;
	Real hi = 4.1;
	int np = 5;
	
	Real d = (hi-lo)/Real(std::max(np-1,1));
	for (int i=0; i<np; i++) {
		pts[i][0] = lo + Real(i)*d;
		pts[i][1] = -1.0;
		pts[i][2] = .25;
	}

	if (1) {
		for (int i=0; i<np; i++) {
			surfidx[i] = surfacesHeight(surfs, pts[i], height[i]);
			surfidx[i] = surfacesDz(surfs, pts[i], dz[i], normal[i] );
			surfs[surfidx[i]]->surfaceNormal(pts[i], -1, normal[i]);
		}

		//print results
		std::cout << "points =\n";
		printnVec3(np,pts,-1,-1);
		std::cout << std::endl;

		std::cout << "height =\n";
		printMatReal(1,np,height,-1,-1);
		std::cout << std::endl;

		std::cout << "dz =\n";
		printMatReal(1,np,dz,-1,-1);
		std::cout << std::endl;

		std::cout << "surface index =\n";
		printMatInt(1,np,surfidx,-1);
		std::cout << std::endl << std::endl;

		std::cout << "surface normals =\n";
		printnVec3(np,normal,-1,-1);
	}

	if (1) {

		//time it

		int n = (int) 1e5;
		clock_t t;
		
		t=clock();
		for (int iter=0; iter<n; iter++) {
			for (int i=0; i<np; i++) { //loop over points
				surfidx[i] = surfacesDz(surfs, pts[i], dz[i], normal[i] );
			}
		}
		t=clock()-t;

		std::cout << "number of points: " << np << std::endl;
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;

	}
	

}


void test_updateWheelContactGeom() {

	const Real rad = .325;
	WheelContactGeom contact;

	//terrain
	SurfaceVector surfs;

	//flat(surfs);
	//ramp(surfs);
	grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );
	
	
	//wheel pose
	HomogeneousTransform HT_wheel_to_world;
	VecOrient orient = {0,0,0}; //assume WMRSIM_USE_QUATERNION 0
	Vec3 pos = {2.5, -1.0, rad};

	//adjust height of wheel
	Real zsurf;
	surfacesHeight(surfs,pos,zsurf);
	pos[2] += zsurf;

	poseToHT(orient,pos,HT_wheel_to_world); //assumes WMRSIM_USE_QUATERNION 0

	updateWheelContactGeom(surfs, HT_wheel_to_world, rad, contact );

	//print
	std::cout << "wheel pos = \n";
	printVec3(pos,-1,-1);
	std::cout << "dz = " << contact.dz << std::endl;
	std::cout << "contact angle (deg)= " << RADTODEG(contact.angle) << std::endl;
	std::cout << "\nHT_contact_to_wheel=\n";
	printHT(contact.HT_wheel,-1,-1);
	std::cout << "\nHT_contact_to_world=\n";
	printHT(contact.HT_world,-1,-1);
	
	if (1) {
		//time it
		int nw = 4;
		clock_t t;

		int n= (int) 1e5;

		t=clock();
		for (int i=0; i<n*nw; i++) {
			updateWheelContactGeom(surfs, HT_wheel_to_world, rad, contact );
		}
		t=clock()-t;
		std::cout << "num wheels: " << nw << std::endl;
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
}

void test_updateTrackContactGeom() {


	const Real rad = .2;
	const Real rad2 = .2;
	const Real L = .5;

	TrackContactGeom contact;
	initTrackContactGeom(rad,rad2,L,contact);
	const int np = contact.get_np();
	

	//terrain
	SurfaceVector surfs;

	//flat(surfs);
	//ramp(surfs);
	grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//track pose
	HomogeneousTransform HT_track_to_world;
	VecOrient orient = {0,0,0}; //assume WMRSIM_USE_QUATERNION 0
	Vec3 pos = {3.0, 0, rad-.01};

	//adjust height of track
	Real zsurf;
	surfacesHeight(surfs,pos,zsurf);
	pos[2] += zsurf;

	poseToHT(orient,pos,HT_track_to_world); 

	updateTrackContactGeom(surfs, HT_track_to_world, contact );

	//print
	std::cout << "track pos = \n";
	printVec3(pos,-1,-1);
	std::cout << "dz=\n"; printMatReal(1,np,contact.dz,-1,-1); std::cout << std::endl;

}


void test_stateToHT() {
	
	const int MAXNS = NUMSTATE(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	//rocky(mdl,state,0);
	zoe(mdl,state,0);

	//get from WmrModel
	int nf = mdl.get_nf();
	int ns = NUMSTATE(nf);
	
	//convert state to homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl, state, HT_parent, HT_world);


	//DEBUGGING, PRINT
	std::cout << "state=\n";
	printMatReal(ns,1,state,-1,-1);
	std::cout << std::endl;

	for (int i=0; i<nf; i++) { 
		std::cout << "HT_parent[" << i << "]=\n";
		printHT(HT_parent[i],-1,-1);
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int i=0; i<nf; i++) { 
		std::cout << "HT_world[" << i << "]=\n";
		printHT(HT_world[i],-1,-1);
		std::cout << std::endl;
	}

	if (1) {
		//time it
		clock_t t;
		int n= (int) 1e6;

		t=clock();
		for (int i=0; i<n; i++) {
			stateToHT(mdl, state, HT_parent, HT_world);
		}
		t=clock()-t;
		std::cout << "stateToHT()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl << std::endl;
	}


}


void test_qvelToQdot() {

	//constants
	const int nf = 5;
	const int ns = NUMSTATE(nf);
	const int nv = NUMQVEL(nf);

	Real statedot[ns];
	Real qvel[nv],qvel_cvtback[nv];
	
	VecEuler euler = {DEGTORAD(10),DEGTORAD(20),DEGTORAD(30)};
	VecOrient orient;
	Mat3 R_body_to_world;
	
	if (WMRSIM_USE_QUATERNION) 
		eulerToQuat(euler,orient);
	else 
		copyEuler(euler,orient);

	orientToRot(orient,R_body_to_world);
	for (int i=0; i<nv; i++) {
		qvel[i] = i;
	}

	qvelToQdot(nf,qvel,orient,R_body_to_world,statedot);
	qdotToQvel(nf,statedot,orient,R_body_to_world,qvel_cvtback); //convert back to check

	//PRINT
	std::cout << "orient =\n";
	printOrient(orient,-1,-1);
	std::cout << std::endl;

	std::cout << "R_body_to_world =\n";
	printMat3(R_body_to_world,-1,-1);
	std::cout << std::endl;

	std::cout << "qvel =\n";
	printMatReal(nv,1,qvel,-1,-1);
	std::cout << std::endl;

	std::cout << "statedot =\n";
	printMatReal(ns,1,statedot,-1,-1);
	std::cout << std::endl;

	std::cout << "qvel (converted back) =\n";
	printMatReal(nv,1,qvel_cvtback,-1,-1);
	std::cout << std::endl;

}


//helper subfunctions
//initialize TrackContactGeom objects
void sub_initTrackContactGeom(const WmrModel& mdl, TrackContactGeom* contacts) {
	
	//get from WmrModel
	const int nt = mdl.get_nt();
	const int* sprocketframeinds = mdl.get_sprocketframeinds();
	const Frame* frames = mdl.get_frames();
		
	for (int tno = 0; tno < nt; tno++) {
		int fi = sprocketframeinds[tno];
		initTrackContactGeom(frames[fi].rad, frames[fi].rad2, frames[fi].L, contacts[tno]);
	}
}

void test_wheelJacobians() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	zoe(mdl,state,0);
	//rocky(mdl,state,0);

	//get from WmrModel
	const int nw = mdl.get_nw();

	//terrain
	SurfaceVector surfs;

	flat(surfs);
	//ramp(surfs);
	//grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//convert state to homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF + WmrModel::MAXNW];
	HomogeneousTransform HT_world[WmrModel::MAXNF + WmrModel::MAXNW];
	stateToHT(mdl,state,HT_parent,HT_world);

	//contact geometry
	WheelContactGeom wcontacts[WmrModel::MAXNW];
	ContactGeom* contacts = static_cast<ContactGeom*>(wcontacts);

	updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);

	//number of points in contact
	int npic = 0;
	for (int wno = 0; wno < nw; wno++)
		npic += wcontacts[wno].incontact;

	//allocate matrix
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	const int MAXNP = WmrModel::MAXNW; //max number of contact points
	Real A[3*MAXNP * MAXNV];

	wheelJacobians(mdl, HT_world, wcontacts, A);

	//DEBUGGING, print
	int ncc = 3*npic; //number of contact constraints
	int nv = NUMQVEL(mdl.get_nf());
	
	std::cout << "A=\n";
	printMatReal(ncc,nv,A,-1,-1);
	std::cout << std::endl;
	
	
	if (1) {
		//time it
		int n= (int) 1e6;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			wheelJacobians(mdl, HT_world, wcontacts, A);
		}
		t=clock()-t;
		std::cout << "wheelJacobians()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
	

}


void test_trackJacobians() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	talon(mdl,state,0);

	//get from WmrModel
	const int nt = mdl.get_nt();

	//terrain
	SurfaceVector surfs;

	flat(surfs);
	//ramp(surfs);
	//grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//convert state to homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF + WmrModel::MAXNW];
	HomogeneousTransform HT_world[WmrModel::MAXNF + WmrModel::MAXNW];
	stateToHT(mdl,state,HT_parent,HT_world);

	//contact geometry
	TrackContactGeom tcontacts[WmrModel::MAXNW];
	sub_initTrackContactGeom(mdl, tcontacts);

	ContactGeom* contacts = static_cast<ContactGeom*>(tcontacts);

	updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);

	//number of points in contact
	int npic = 0;
	for (int tno = 0; tno < nt; tno++) {
		npic += logicalCount(tcontacts[tno].get_np(), tcontacts[tno].incontact);
	}

	//allocate matrix
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	const int MAXNP = WmrModel::MAXNT * ContactGeom::MAXNP; //max number of contact points
	Real A[3*MAXNP * MAXNV];

	trackJacobians(mdl, HT_world, tcontacts, A);

	//DEBUGGING, print
	int ncc = 3*npic; //number of contact constraints
	int nv = NUMQVEL(mdl.get_nf());
	
	std::cout << "A=\n";
	printMatReal(ncc,nv,A,-1,-1);
	std::cout << std::endl;
	
	
	if (1) {
		//time it
		int n= (int) 1e6;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			trackJacobians(mdl, HT_world, tcontacts, A);
		}
		t=clock()-t;
		std::cout << "trackJacobians()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
	

}


void test_forwardVelKin() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	zoe(mdl,state,0);
	//rocky(mdl,state,0);
	//talon(mdl,state,0); //TODO

	//get from WmrModel
	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	const int na = mdl.get_na();
	const int nv = NUMQVEL(mdl.get_nf());

	//terrain
	SurfaceVector surfs;

	flat(surfs);
	//ramp(surfs);
	//grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//init contact geometry
	WheelContactGeom wcontacts[WmrModel::MAXNW];
	TrackContactGeom tcontacts[WmrModel::MAXNW];
	ContactGeom* contacts = 0;

	if (nw > 0) {
		contacts = static_cast<ContactGeom*>(wcontacts);
	} else if (nt > 0) {
		sub_initTrackContactGeom(mdl, tcontacts);
		contacts = static_cast<ContactGeom*>(tcontacts);
	}

	//convert state to homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,state,HT_parent,HT_world);

	//update contact geometry
	updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);

	//controller inputs
	Real u[WmrModel::MAXNA];
	Real qvel_cmd[MAXNV];
	Real time = 0;

	mdl.controller(mdl, 0, state, u, qvel_cmd);

	//allocate outputs
	Real qvel[MAXNV];
	Vec3 vc[WmrModel::MAXNW];

	//compute joint space velocity
	forwardVelKin(mdl, state, u, HT_world, contacts, qvel, vc);
	
	//PRINT
	std::cout << "u =\n";
	printMatReal(na,1,u,-1,-1);
	std::cout << std::endl;

	std::cout << "qvel_cmd =\n";
	printMatReal(nv,1,qvel_cmd,-1,-1);
	std::cout << std::endl;

	std::cout << "qvel =\n";
	printMatReal(nv,1,qvel,-1,-1);
	std::cout << std::endl;

	std::cout << "vc =\n";
	printnVec3(nw,vc,-1,-1);
	std::cout << std::endl;

	if (1) {
		//time it
		int n= (int) 1e5;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			forwardVelKin(mdl, state, u, HT_world, contacts, qvel, 0);
		}
		t=clock()-t;
		std::cout << "wheelJacobians()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
	

}


void test_initTerrainContact() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	zoe(mdl,state,0);
	//rocky(mdl,state,0);
	//talon(mdl,state,0); //TODO

	//get from WmrModel
	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	int nf = mdl.get_nf();
	int ns = NUMSTATE(nf);	

	//terrain
	SurfaceVector surfs;

	//flat(surfs);
	//ramp(surfs);
	grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//init contact geometry
	WheelContactGeom wcontacts[WmrModel::MAXNW];
	TrackContactGeom tcontacts[WmrModel::MAXNW];
	ContactGeom* contacts = 0;

	if (nw > 0) {
		contacts = static_cast<ContactGeom*>(wcontacts);
	} else if (nt > 0) {
		sub_initTrackContactGeom(mdl, tcontacts);
		contacts = static_cast<ContactGeom*>(tcontacts);
	}

	//backup
	Real state_before[MAXNS];
	copyVec(ns,state,state_before);

	std::cout << "state (before) =\n"; printMatReal(ns,1,state_before,-1,-1); std::cout << std::endl;

	initTerrainContact(mdl, surfs, contacts, state);
	
	std::cout << "state (after) =\n"; printMatReal(ns,1,state,-1,-1); std::cout << std::endl;

	if (1) {
		//time it
		int n= (int) 1e4;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			copyVec(ns,state_before,state);
			initTerrainContact(mdl, surfs, contacts, state);
		}
		t=clock()-t;
		std::cout << "initTerrainContact()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}

}


void test_subtreeInertias() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	zoe(mdl,state,0);
	//rocky(mdl,state,0);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const Frame* frames = mdl.get_frames();

	//homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	//HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,state,HT_parent,0);

	//Plucker transforms
	Mat6b Xup[WmrModel::MAXNF];
	for (int fi=0; fi < nf; fi++)
		invHTToPlucker(HT_parent[fi], Xup[fi]);

	//subtree inertias
	Mat6b Is_subt[WmrModel::MAXNF];
	subtreeInertias(mdl, Xup, Is_subt);

	//PRINT
	//for (int fi=0; fi<nf; fi++) {
	//	std::cout << "Is[" << fi << "]=\n";
	//	printMat6b(frames[fi].get_Is(),-1,-1);
	//	std::cout << std::endl;
	//}

	//for (int fi=1; fi<nf; fi++) {
	//	std::cout << "Xup[" << fi << "]=\n";
	//	printMat6b(Xup[fi],-1,-1);
	//	std::cout << std::endl;
	//}

	for (int fi=0; fi<nf; fi++) {
		std::cout << "Is (subtree) [" << fi << "]=\n";
		printMat6b(Is_subt[fi],-1,-1);
		std::cout << std::endl;
	}

	if (1) {
		//time it
		int n= (int) 1e6;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			subtreeInertias(mdl, Xup, Is_subt);
		}
		t=clock()-t;
		std::cout << "subtreeInertias()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
}


void test_jointSpaceInertia() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	
	zoe(mdl,state,0);
	//rocky(mdl,state,0);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	//homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,state,HT_parent,HT_world);

	//Plucker transforms
	Mat6b Xup[WmrModel::MAXNF];
	for (int fi=0; fi < nf; fi++)
		invHTToPlucker(HT_parent[fi], Xup[fi]);

	//subtree inertias
	Mat6b Is_subt[WmrModel::MAXNF];
	subtreeInertias(mdl, Xup, Is_subt);

	//joint space inertia
	Real H[MAXNV * MAXNV];

	jointSpaceInertia(mdl,Xup,Is_subt,H);

	//PRINT
	std::cout << "H=\n"; printMatReal(nv,nv,H,-1,-1);

	if (1) {
		//time it
		int n= (int) 1e6;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			jointSpaceInertia(mdl,Xup,Is_subt,H);
		}
		t=clock()-t;
		std::cout << "jointSpaceInertia()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
}



void test_jointSpaceBiasForce() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	Real qvel[MAXNV];
	
	zoe(mdl,state,qvel);
	//rocky(mdl,state,qvel);

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nv = NUMQVEL(nf);

	//homogeneous transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF];
	HomogeneousTransform HT_world[WmrModel::MAXNF];
	stateToHT(mdl,state,HT_parent,HT_world);

	//Plucker transforms
	Mat6b Xup[WmrModel::MAXNF];
	for (int fi=0; fi < nf; fi++)
		invHTToPlucker(HT_parent[fi], Xup[fi]);

	//joint bias force
	//DEBUGGING
	for (int vi=0; vi<nv; vi++) 
		qvel[vi]=vi+1;
	
	Real C[MAXNV];
	jointSpaceBiasForce(mdl,Xup,qvel,C);

	//PRINT
	std::cout << "qvel=\n"; printMatReal(nv,1,qvel,-1,-1);
	std::cout << "C=\n"; printMatReal(nv,1,C,-1,-1);

	if (1) {
		//time it
		int n= (int) 1e6;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			jointSpaceBiasForce(mdl,Xup,qvel,C);
		}
		t=clock()-t;
		std::cout << "jointSpaceBiasForce()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
}



void test_forwardDyn() {

	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	Real qvel[MAXNV];
	
	zoe(mdl,state,qvel);
	//rocky(mdl,state,qvel);
	//talon(mdl,state,qvel); //TODO

	//setVec(nv,0.0,qvel); //DEBUGGING
	mdl.actuatorModel=0; //ideal actuators

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int ns = NUMSTATE(nf);
	const int nv = NUMQVEL(nf);

	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();

	//terrain
	SurfaceVector surfs;

	flat(surfs);
	//ramp(surfs);
	//grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//init contact geometry
	WheelContactGeom wcontacts[WmrModel::MAXNW];
	TrackContactGeom tcontacts[WmrModel::MAXNW];
	ContactGeom* contacts = 0;

	if (nw > 0) {
		contacts = static_cast<ContactGeom*>(wcontacts);
	} else if (nt > 0) {
		sub_initTrackContactGeom(mdl, tcontacts);
		contacts = static_cast<ContactGeom*>(tcontacts);
	}

	//initTerrainContact(mdl, surfs, contacts, state); //DEBUGGING

	//convert state to Homogeneous Transforms
	HomogeneousTransform HT_parent[WmrModel::MAXNF + WmrModel::MAXNW];
	HomogeneousTransform HT_world[WmrModel::MAXNF + WmrModel::MAXNW];
	stateToHT(mdl,state,HT_parent,HT_world);

	//update contact geometry
	updateModelContactGeom(mdl, surfs, HT_world, 0, contacts);

	//controller inputs
	ControllerIO u;
	Real time = 0;
	
	mdl.controller(mdl, 0, state, u.cmd, 0);

	//additional inputs
	setVec(mdl.get_na(),0.0,u.interr);

	Real dt = .04; //time step

	//outputs
	Real qacc[MAXNV];

	forwardDyn(mdl, state, qvel, u, HT_parent, HT_world, contacts, dt, qacc);

	//DEBUGGING, print
	std::cout << "state=\n"; printMatReal(ns,1,state,-1,-1); std::cout << std::endl;
	std::cout << "qvel=\n"; printMatReal(nv,1,qvel,-1,-1); std::cout << std::endl;
	std::cout << "qacc=\n"; printMatReal(nv,1,qacc,-1,-1); std::cout << std::endl;


	if (1) {
		//time it
		int n= (int) 1e4;
		clock_t t;

		t=clock();
		for (int i=0; i<n; i++) {
			forwardDyn(mdl, state, qvel, u, HT_parent, HT_world, contacts, dt, qacc);
		}
		t=clock()-t;
		std::cout << "forwardDyn()\n";
		std::cout << "iterations: " << (Real) n << std::endl;
		std::cout << "clock (ms): " << t << std::endl;
	}
}


void test_simulate() {

	//options
	bool do_dyn = true; //do dynamic simulation, else kinematic
//	bool use_erp_cfm = false;
	bool ideal_actuators = true;
	bool do_anim = true; //do animation
	bool do_time = true;

	const Real dt = .04;
	const int nsteps = (int) floor(10.0/dt);
	Real time = 0;

	//for allocation
	const int MAXNS = NUMSTATE(WmrModel::MAXNF);
	const int MAXNV = NUMQVEL(WmrModel::MAXNF);
	const int MAXNY = MAXNS+MAXNV+WmrModel::MAXNA; //for dynamic sim

	//make WmrModel object
	WmrModel mdl;
	Real state[MAXNS];
	Real qvel[MAXNV]; //for dynamic sim

	//uncomment one of the following:
	zoe(mdl,state,qvel);
//	rocky(mdl,state,qvel);
//	talon(mdl,state,qvel);

	//also uncomment the corresponding scene function below!

	//initialize wheel-ground contact model
	mdl.wheelGroundContactModel(0, mdl.wgc_p, 0, 0, 0, //inputs
		0, 0); //outputs

//	if (use_erp_cfm)
//		mdl.wheelGroundContactModel=0;
	if (ideal_actuators)
		mdl.actuatorModel=0;

	//get from WmrModel
	const int nf = mdl.get_nf();
	const int nw = mdl.get_nw();
	const int nt = mdl.get_nt();
	const int ns = NUMSTATE(nf); //number of states
	//for dynamic sim
	const int nv = NUMQVEL(nf); //size of joint space vel
	const int na = mdl.get_na(); 
	int ny;

	//terrain
	SurfaceVector surfs;

	//flat(surfs);
	//ramp(surfs); //must also uncomment flat

	grid(surfs, ResourceDir() + std::string("gridsurfdata.txt") );

	//init contact geometry
	WheelContactGeom wcontacts[WmrModel::MAXNW];
	TrackContactGeom tcontacts[WmrModel::MAXNW];
	ContactGeom* contacts =0; //base class

	if (nw > 0) {
		contacts = static_cast<ContactGeom*>(wcontacts);
	} else if (nt > 0) {
		sub_initTrackContactGeom(mdl, tcontacts);
		contacts = static_cast<ContactGeom*>(tcontacts);
	}

	initTerrainContact(mdl, surfs, contacts, state); //DEBUGGING


	//allocate
	Real y[MAXNY];
	Real ydot[MAXNY];

	//init y
	if (do_dyn) { //dynamic sim
		copyVec(ns,state,y);
		copyVec(nv,qvel,y+ns);
		setVec(na,0.0,y+ns+nv); //interr
		ny = ns + nv + na;
	} else {
		copyVec(ns,state,y);
		ny = ns;
	}

	//backup
	Real y0[MAXNY];
	copyVec(ny,y,y0); 

	//allocate
	HomogeneousTransform HT_parent[WmrModel::MAXNF];

#if WMRSIM_ENABLE_ANIMATION
	WmrAnimation anim;
	if (do_anim) { //animate
		anim.start();

		//uncomment the scene function that corresponds to the model function above

		zoeScene(mdl, anim);
//		rockyScene(mdl, anim);
//		talonScene(mdl, tcontacts, anim);

		for (int i=0; i<surfs.size(); i++)
			anim.addEntitySurface(surfs[i].get());
	}
#endif

	std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,y,-1,-1); std::cout << std::endl;

	for (int i=0; i<nsteps; i++) {

		if (do_dyn) {
			odeDyn(time, y, mdl, surfs, contacts, dt, ydot, HT_parent);
		} else {
			odeKin(time, y, mdl, surfs, contacts, ydot, HT_parent);
		}
		addmVec(ny,ydot,dt,y);
		time += dt;

#if WMRSIM_ENABLE_ANIMATION
		if (do_anim) {
			anim.updateNodesLines(nf, HT_parent, nw + nt, contacts);

			if (!anim.updateRender())
				goto stop;

			while (anim.get_mPause()) {
				if (!anim.updateRender())
					goto stop;
			}
		}
#endif
	}

	std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,y,-1,-1); std::cout << std::endl;

#if WMRSIM_ENABLE_ANIMATION
	if (do_anim) {
		//DEBUGGING, wait for escape key before exiting
		while (1) {
			if (!anim.updateRender())
				goto stop;
		}
	}
	stop: //goto
#endif
	
	if (do_time) {
		//time it
		int n= (int) 100;
//		clock_t t;

		time = 0;

//		t=clock();

	  timeval t0,t1;
	  gettimeofday(&t0, NULL);
		for (int iter=0; iter<n; iter++) {
			copyVec(ny,y0,y); //reset to backup

			for (int i=0; i<nsteps; i++) {
				if (do_dyn) {
					odeDyn(time, y, mdl, surfs, contacts, dt, ydot, HT_parent);
				} else {
					odeKin(time, y, mdl, surfs, contacts, ydot, HT_parent);
				}
				addmVec(ny,ydot,dt,y);
				time += dt;
			}
			//std::cout << "state(" << time << ")=\n"; printMatReal(ns,1,y,-1,-1); std::cout << std::endl;
		}
    gettimeofday(&t1, NULL);

//		t=clock()-t;

		std::cout << "simulate\n";
		std::cout << "iterations: " << (Real) n << std::endl;
//		std::cout << "clock (ms): " << t << std::endl;
		std::cout << "clock (sec): " << tosec(t1)-tosec(t0) << std::endl;
	}

}
