#include <wmrde/algebra/spatial.h>

namespace wmrde
{

//convert to spatial inertia
//m:	mass
//c:	center of mass location
//I:	moment of inertia
//Is:	spatial inertia
//speed not critical, only called when making WmrModel
void toSpatialInertia(const Real m, const Vec3 c, const Mat3 I, Mat6b Is) {
	//MATLAB:
	//C = skew3(c)
	//Is = [ I + m*(C*C'), m*C
	//		 m*C',		   m*eye(3) ]
	
	Mat3 C,CT;
	skewVec3(c,C);
	copyTMat3(C,CT);

	multMatMat3(C,CT,Is);
	addmMat3(I,m,Is,Is);
	mulcMat3(CT,m,Is+BLOCK1); 
	mulcMat3(C,m,Is+BLOCK2); 
	setMat3Diagonal(m,m,m,Is+BLOCK3);
}

//convert from spatial inertia
void fromSpatialInertia(const Mat6b Is, Real m, Vec3 c, Mat3 I) {
	
	Mat3 C,CT;

	m=(Is+BLOCK3)[0];
	mulcMat3(Is+BLOCK2,1.0/m,C); 
	mulcMat3(Is+BLOCK1,1.0/m,CT);
	unskewVec3(C,c);
	multMatMat3(C,CT,I); //temporary
	mulcMat3(I,m,I);
	addmMat3(Is,-1,I,I);
}

//transform the coordinates of spatial inertia
//R = P*I*P'
//Plucker transforms & spatial inertia have special structure
//code generated using MATLAB symbolic toolbox
//sym_transform_spatial_inertia.m
void multPluckerTInertiaPlucker(const Mat6b P, const Mat6b I, Mat6b R) {

	Real I1_11 = I[0];
	Real I1_12 = I[1];
	Real I1_13 = I[2];
	Real I2_21 = I[13];
	Real I2_31 = I[14];
	Real I1_22 = I[5];
	Real I1_23 = I[6];
	Real I2_12 = I[16];
	Real I2_32 = I[18];
	Real I1_33 = I[10];
	Real I2_13 = I[20];
	Real I2_23 = I[21];
	Real I4_11 = I[36];

	Real P1_11 = P[0];
	Real P1_21 = P[1];
	Real P1_31 = P[2];
	Real P2_11 = P[12];
	Real P2_21 = P[13];
	Real P2_31 = P[14];
	Real P1_12 = P[4];
	Real P1_22 = P[5];
	Real P1_32 = P[6];
	Real P2_12 = P[16];
	Real P2_22 = P[17];
	Real P2_32 = P[18];
	Real P1_13 = P[8];
	Real P1_23 = P[9];
	Real P1_33 = P[10];
	Real P2_13 = P[20];
	Real P2_23 = P[21];
	Real P2_33 = P[22];

	Real t2,t3,t4,t5,t6,t7,t8,t9,
	t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,
	t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,
	t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,
	t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,
	t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,
	t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,
	t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,
	t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,
	t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,
	t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,
	t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,
	t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,
	t130;

	t2 = I1_11*P1_11;
	t3 = I1_12*P1_21;
	t4 = I1_13*P1_31;
	t5 = I2_21*P2_21;
	t6 = I2_31*P2_31;
	t7 = t2+t3+t4+t5+t6;
	t8 = I1_12*P1_11;
	t9 = I1_22*P1_21;
	t10 = I1_23*P1_31;
	t11 = I2_12*P2_11;
	t12 = I2_32*P2_31;
	t13 = t8+t9+t10+t11+t12;
	t14 = I1_13*P1_11;
	t15 = I1_23*P1_21;
	t16 = I1_33*P1_31;
	t17 = I2_13*P2_11;
	t18 = I2_23*P2_21;
	t19 = t14+t15+t16+t17+t18;
	t20 = I2_12*P1_21;
	t21 = I2_13*P1_31;
	t22 = I4_11*P2_11;
	t23 = t20+t21+t22;
	t24 = I2_21*P1_11;
	t25 = I2_23*P1_31;
	t26 = I4_11*P2_21;
	t27 = t24+t25+t26;
	t28 = I2_31*P1_11;
	t29 = I2_32*P1_21;
	t30 = I4_11*P2_31;
	t31 = t28+t29+t30;
	t32 = I1_11*P1_12;
	t33 = I1_12*P1_22;
	t34 = I1_13*P1_32;
	t35 = I2_21*P2_22;
	t36 = I2_31*P2_32;
	t37 = t32+t33+t34+t35+t36;
	t38 = I1_12*P1_12;
	t39 = I1_22*P1_22;
	t40 = I1_23*P1_32;
	t41 = I2_12*P2_12;
	t42 = I2_32*P2_32;
	t43 = t38+t39+t40+t41+t42;
	t44 = I1_13*P1_12;
	t45 = I1_23*P1_22;
	t46 = I1_33*P1_32;
	t47 = I2_13*P2_12;
	t48 = I2_23*P2_22;
	t49 = t44+t45+t46+t47+t48;
	t50 = I2_12*P1_22;
	t51 = I2_13*P1_32;
	t52 = I4_11*P2_12;
	t53 = t50+t51+t52;
	t54 = I2_21*P1_12;
	t55 = I2_23*P1_32;
	t56 = I4_11*P2_22;
	t57 = t54+t55+t56;
	t58 = I2_31*P1_12;
	t59 = I2_32*P1_22;
	t60 = I4_11*P2_32;
	t61 = t58+t59+t60;
	t62 = I1_11*P1_13;
	t63 = I1_12*P1_23;
	t64 = I1_13*P1_33;
	t65 = I2_21*P2_23;
	t66 = I2_31*P2_33;
	t67 = t62+t63+t64+t65+t66;
	t68 = I1_12*P1_13;
	t69 = I1_22*P1_23;
	t70 = I1_23*P1_33;
	t71 = I2_12*P2_13;
	t72 = I2_32*P2_33;
	t73 = t68+t69+t70+t71+t72;
	t74 = I1_13*P1_13;
	t75 = I1_23*P1_23;
	t76 = I1_33*P1_33;
	t77 = I2_13*P2_13;
	t78 = I2_23*P2_23;
	t79 = t74+t75+t76+t77+t78;
	t80 = I2_12*P1_23;
	t81 = I2_13*P1_33;
	t82 = I4_11*P2_13;
	t83 = t80+t81+t82;
	t84 = I2_21*P1_13;
	t85 = I2_23*P1_33;
	t86 = I4_11*P2_23;
	t87 = t84+t85+t86;
	t88 = I2_31*P1_13;
	t89 = I2_32*P1_23;
	t90 = I4_11*P2_33;
	t91 = t88+t89+t90;
	t92 = I2_13*P1_11;
	t93 = I2_23*P1_21;
	t94 = t92+t93;
	t95 = I2_12*P1_11;
	t96 = I2_32*P1_31;
	t97 = t95+t96;
	t98 = I2_21*P1_21;
	t99 = I2_31*P1_31;
	t100 = t98+t99;
	t101 = I2_13*P1_12;
	t102 = I2_23*P1_22;
	t103 = t101+t102;
	t104 = I2_12*P1_12;
	t105 = I2_32*P1_32;
	t106 = t104+t105;
	t107 = I2_21*P1_22;
	t108 = I2_31*P1_32;
	t109 = t107+t108;
	t110 = I4_11*P1_11*P1_12;
	t111 = I4_11*P1_21*P1_22;
	t112 = I4_11*P1_31*P1_32;
	t113 = t110+t111+t112;
	t114 = I2_13*P1_13;
	t115 = I2_23*P1_23;
	t116 = t114+t115;
	t117 = I2_12*P1_13;
	t118 = I2_32*P1_33;
	t119 = t117+t118;
	t120 = I2_21*P1_23;
	t121 = I2_31*P1_33;
	t122 = t120+t121;
	t123 = I4_11*P1_11*P1_13;
	t124 = I4_11*P1_21*P1_23;
	t125 = I4_11*P1_31*P1_33;
	t126 = t123+t124+t125;
	t127 = I4_11*P1_12*P1_13;
	t128 = I4_11*P1_22*P1_23;
	t129 = I4_11*P1_32*P1_33;
	t130 = t127+t128+t129;

	R[0] = P1_11*t7+P1_21*t13+P1_31*t19+P2_11*t23+P2_21*t27+P2_31*t31;
	R[4] = P1_12*t7+P1_22*t13+P1_32*t19+P2_12*t23+P2_22*t27+P2_32*t31;
	R[8] = P1_13*t7+P1_23*t13+P1_33*t19+P2_13*t23+P2_23*t27+P2_33*t31;
	R[24] = P1_11*t23+P1_21*t27+P1_31*t31;
	R[28] = P1_12*t23+P1_22*t27+P1_32*t31;
	R[32] = P1_13*t23+P1_23*t27+P1_33*t31;
	R[1] = P1_11*t37+P1_21*t43+P1_31*t49+P2_11*t53+P2_21*t57+P2_31*t61;
	R[5] = P1_12*t37+P1_22*t43+P1_32*t49+P2_12*t53+P2_22*t57+P2_32*t61;
	R[9] = P1_13*t37+P1_23*t43+P1_33*t49+P2_13*t53+P2_23*t57+P2_33*t61;
	R[25] = P1_11*t53+P1_21*t57+P1_31*t61;
	R[29] = P1_12*t53+P1_22*t57+P1_32*t61;
	R[33] = P1_13*t53+P1_23*t57+P1_33*t61;
	R[2] = P1_11*t67+P1_21*t73+P1_31*t79+P2_11*t83+P2_21*t87+P2_31*t91;
	R[6] = P1_12*t67+P1_22*t73+P1_32*t79+P2_12*t83+P2_22*t87+P2_32*t91;
	R[10] = P1_13*t67+P1_23*t73+P1_33*t79+P2_13*t83+P2_23*t87+P2_33*t91;
	R[26] = P1_11*t83+P1_21*t87+P1_31*t91;
	R[30] = P1_12*t83+P1_22*t87+P1_32*t91;
	R[34] = P1_13*t83+P1_23*t87+P1_33*t91;
	R[12] = P1_11*t100+P1_21*t97+P1_31*t94+I4_11*P1_11*P2_11+I4_11*P1_21*P2_21+I4_11*P1_31*P2_31;
	R[16] = P1_12*t100+P1_22*t97+P1_32*t94+I4_11*P1_11*P2_12+I4_11*P1_21*P2_22+I4_11*P1_31*P2_32;
	R[20] = P1_13*t100+P1_23*t97+P1_33*t94+I4_11*P1_11*P2_13+I4_11*P1_21*P2_23+I4_11*P1_31*P2_33;
	R[36] = I4_11*(P1_11*P1_11)+I4_11*(P1_21*P1_21)+I4_11*(P1_31*P1_31);
	R[40] = t113;
	R[44] = t126;
	R[13] = P1_11*t109+P1_21*t106+P1_31*t103+I4_11*P1_12*P2_11+I4_11*P1_22*P2_21+I4_11*P1_32*P2_31;
	R[17] = P1_12*t109+P1_22*t106+P1_32*t103+I4_11*P1_12*P2_12+I4_11*P1_22*P2_22+I4_11*P1_32*P2_32;
	R[21] = P1_13*t109+P1_23*t106+P1_33*t103+I4_11*P1_12*P2_13+I4_11*P1_22*P2_23+I4_11*P1_32*P2_33;
	R[37] = t113;
	R[41] = I4_11*(P1_12*P1_12)+I4_11*(P1_22*P1_22)+I4_11*(P1_32*P1_32);
	R[45] = t130;
	R[14] = P1_11*t122+P1_21*t119+P1_31*t116+I4_11*P1_13*P2_11+I4_11*P1_23*P2_21+I4_11*P1_33*P2_31;
	R[18] = P1_12*t122+P1_22*t119+P1_32*t116+I4_11*P1_13*P2_12+I4_11*P1_23*P2_22+I4_11*P1_33*P2_32;
	R[22] = P1_13*t122+P1_23*t119+P1_33*t116+I4_11*P1_13*P2_13+I4_11*P1_23*P2_23+I4_11*P1_33*P2_33;
	R[38] = t126;
	R[42] = t130;
	R[46] = I4_11*(P1_13*P1_13)+I4_11*(P1_23*P1_23)+I4_11*(P1_33*P1_33);
}

} //namespace
