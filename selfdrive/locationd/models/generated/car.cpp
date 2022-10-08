#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6892551631316742231) {
   out_6892551631316742231[0] = delta_x[0] + nom_x[0];
   out_6892551631316742231[1] = delta_x[1] + nom_x[1];
   out_6892551631316742231[2] = delta_x[2] + nom_x[2];
   out_6892551631316742231[3] = delta_x[3] + nom_x[3];
   out_6892551631316742231[4] = delta_x[4] + nom_x[4];
   out_6892551631316742231[5] = delta_x[5] + nom_x[5];
   out_6892551631316742231[6] = delta_x[6] + nom_x[6];
   out_6892551631316742231[7] = delta_x[7] + nom_x[7];
   out_6892551631316742231[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_632344986528208149) {
   out_632344986528208149[0] = -nom_x[0] + true_x[0];
   out_632344986528208149[1] = -nom_x[1] + true_x[1];
   out_632344986528208149[2] = -nom_x[2] + true_x[2];
   out_632344986528208149[3] = -nom_x[3] + true_x[3];
   out_632344986528208149[4] = -nom_x[4] + true_x[4];
   out_632344986528208149[5] = -nom_x[5] + true_x[5];
   out_632344986528208149[6] = -nom_x[6] + true_x[6];
   out_632344986528208149[7] = -nom_x[7] + true_x[7];
   out_632344986528208149[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6681531782435372252) {
   out_6681531782435372252[0] = 1.0;
   out_6681531782435372252[1] = 0;
   out_6681531782435372252[2] = 0;
   out_6681531782435372252[3] = 0;
   out_6681531782435372252[4] = 0;
   out_6681531782435372252[5] = 0;
   out_6681531782435372252[6] = 0;
   out_6681531782435372252[7] = 0;
   out_6681531782435372252[8] = 0;
   out_6681531782435372252[9] = 0;
   out_6681531782435372252[10] = 1.0;
   out_6681531782435372252[11] = 0;
   out_6681531782435372252[12] = 0;
   out_6681531782435372252[13] = 0;
   out_6681531782435372252[14] = 0;
   out_6681531782435372252[15] = 0;
   out_6681531782435372252[16] = 0;
   out_6681531782435372252[17] = 0;
   out_6681531782435372252[18] = 0;
   out_6681531782435372252[19] = 0;
   out_6681531782435372252[20] = 1.0;
   out_6681531782435372252[21] = 0;
   out_6681531782435372252[22] = 0;
   out_6681531782435372252[23] = 0;
   out_6681531782435372252[24] = 0;
   out_6681531782435372252[25] = 0;
   out_6681531782435372252[26] = 0;
   out_6681531782435372252[27] = 0;
   out_6681531782435372252[28] = 0;
   out_6681531782435372252[29] = 0;
   out_6681531782435372252[30] = 1.0;
   out_6681531782435372252[31] = 0;
   out_6681531782435372252[32] = 0;
   out_6681531782435372252[33] = 0;
   out_6681531782435372252[34] = 0;
   out_6681531782435372252[35] = 0;
   out_6681531782435372252[36] = 0;
   out_6681531782435372252[37] = 0;
   out_6681531782435372252[38] = 0;
   out_6681531782435372252[39] = 0;
   out_6681531782435372252[40] = 1.0;
   out_6681531782435372252[41] = 0;
   out_6681531782435372252[42] = 0;
   out_6681531782435372252[43] = 0;
   out_6681531782435372252[44] = 0;
   out_6681531782435372252[45] = 0;
   out_6681531782435372252[46] = 0;
   out_6681531782435372252[47] = 0;
   out_6681531782435372252[48] = 0;
   out_6681531782435372252[49] = 0;
   out_6681531782435372252[50] = 1.0;
   out_6681531782435372252[51] = 0;
   out_6681531782435372252[52] = 0;
   out_6681531782435372252[53] = 0;
   out_6681531782435372252[54] = 0;
   out_6681531782435372252[55] = 0;
   out_6681531782435372252[56] = 0;
   out_6681531782435372252[57] = 0;
   out_6681531782435372252[58] = 0;
   out_6681531782435372252[59] = 0;
   out_6681531782435372252[60] = 1.0;
   out_6681531782435372252[61] = 0;
   out_6681531782435372252[62] = 0;
   out_6681531782435372252[63] = 0;
   out_6681531782435372252[64] = 0;
   out_6681531782435372252[65] = 0;
   out_6681531782435372252[66] = 0;
   out_6681531782435372252[67] = 0;
   out_6681531782435372252[68] = 0;
   out_6681531782435372252[69] = 0;
   out_6681531782435372252[70] = 1.0;
   out_6681531782435372252[71] = 0;
   out_6681531782435372252[72] = 0;
   out_6681531782435372252[73] = 0;
   out_6681531782435372252[74] = 0;
   out_6681531782435372252[75] = 0;
   out_6681531782435372252[76] = 0;
   out_6681531782435372252[77] = 0;
   out_6681531782435372252[78] = 0;
   out_6681531782435372252[79] = 0;
   out_6681531782435372252[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5873591982632034770) {
   out_5873591982632034770[0] = state[0];
   out_5873591982632034770[1] = state[1];
   out_5873591982632034770[2] = state[2];
   out_5873591982632034770[3] = state[3];
   out_5873591982632034770[4] = state[4];
   out_5873591982632034770[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5873591982632034770[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5873591982632034770[7] = state[7];
   out_5873591982632034770[8] = state[8];
}
void F_fun(double *state, double dt, double *out_5222924237119091959) {
   out_5222924237119091959[0] = 1;
   out_5222924237119091959[1] = 0;
   out_5222924237119091959[2] = 0;
   out_5222924237119091959[3] = 0;
   out_5222924237119091959[4] = 0;
   out_5222924237119091959[5] = 0;
   out_5222924237119091959[6] = 0;
   out_5222924237119091959[7] = 0;
   out_5222924237119091959[8] = 0;
   out_5222924237119091959[9] = 0;
   out_5222924237119091959[10] = 1;
   out_5222924237119091959[11] = 0;
   out_5222924237119091959[12] = 0;
   out_5222924237119091959[13] = 0;
   out_5222924237119091959[14] = 0;
   out_5222924237119091959[15] = 0;
   out_5222924237119091959[16] = 0;
   out_5222924237119091959[17] = 0;
   out_5222924237119091959[18] = 0;
   out_5222924237119091959[19] = 0;
   out_5222924237119091959[20] = 1;
   out_5222924237119091959[21] = 0;
   out_5222924237119091959[22] = 0;
   out_5222924237119091959[23] = 0;
   out_5222924237119091959[24] = 0;
   out_5222924237119091959[25] = 0;
   out_5222924237119091959[26] = 0;
   out_5222924237119091959[27] = 0;
   out_5222924237119091959[28] = 0;
   out_5222924237119091959[29] = 0;
   out_5222924237119091959[30] = 1;
   out_5222924237119091959[31] = 0;
   out_5222924237119091959[32] = 0;
   out_5222924237119091959[33] = 0;
   out_5222924237119091959[34] = 0;
   out_5222924237119091959[35] = 0;
   out_5222924237119091959[36] = 0;
   out_5222924237119091959[37] = 0;
   out_5222924237119091959[38] = 0;
   out_5222924237119091959[39] = 0;
   out_5222924237119091959[40] = 1;
   out_5222924237119091959[41] = 0;
   out_5222924237119091959[42] = 0;
   out_5222924237119091959[43] = 0;
   out_5222924237119091959[44] = 0;
   out_5222924237119091959[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5222924237119091959[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5222924237119091959[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5222924237119091959[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5222924237119091959[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5222924237119091959[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5222924237119091959[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5222924237119091959[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5222924237119091959[53] = -9.8000000000000007*dt;
   out_5222924237119091959[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5222924237119091959[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5222924237119091959[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5222924237119091959[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5222924237119091959[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5222924237119091959[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5222924237119091959[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5222924237119091959[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5222924237119091959[62] = 0;
   out_5222924237119091959[63] = 0;
   out_5222924237119091959[64] = 0;
   out_5222924237119091959[65] = 0;
   out_5222924237119091959[66] = 0;
   out_5222924237119091959[67] = 0;
   out_5222924237119091959[68] = 0;
   out_5222924237119091959[69] = 0;
   out_5222924237119091959[70] = 1;
   out_5222924237119091959[71] = 0;
   out_5222924237119091959[72] = 0;
   out_5222924237119091959[73] = 0;
   out_5222924237119091959[74] = 0;
   out_5222924237119091959[75] = 0;
   out_5222924237119091959[76] = 0;
   out_5222924237119091959[77] = 0;
   out_5222924237119091959[78] = 0;
   out_5222924237119091959[79] = 0;
   out_5222924237119091959[80] = 1;
}
void h_25(double *state, double *unused, double *out_7830725943130014969) {
   out_7830725943130014969[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8936357706928973574) {
   out_8936357706928973574[0] = 0;
   out_8936357706928973574[1] = 0;
   out_8936357706928973574[2] = 0;
   out_8936357706928973574[3] = 0;
   out_8936357706928973574[4] = 0;
   out_8936357706928973574[5] = 0;
   out_8936357706928973574[6] = 1;
   out_8936357706928973574[7] = 0;
   out_8936357706928973574[8] = 0;
}
void h_24(double *state, double *unused, double *out_7215733804563166273) {
   out_7215733804563166273[0] = state[4];
   out_7215733804563166273[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7337736767775078476) {
   out_7337736767775078476[0] = 0;
   out_7337736767775078476[1] = 0;
   out_7337736767775078476[2] = 0;
   out_7337736767775078476[3] = 0;
   out_7337736767775078476[4] = 1;
   out_7337736767775078476[5] = 0;
   out_7337736767775078476[6] = 0;
   out_7337736767775078476[7] = 0;
   out_7337736767775078476[8] = 0;
   out_7337736767775078476[9] = 0;
   out_7337736767775078476[10] = 0;
   out_7337736767775078476[11] = 0;
   out_7337736767775078476[12] = 0;
   out_7337736767775078476[13] = 0;
   out_7337736767775078476[14] = 1;
   out_7337736767775078476[15] = 0;
   out_7337736767775078476[16] = 0;
   out_7337736767775078476[17] = 0;
}
void h_30(double *state, double *unused, double *out_4171406999558043376) {
   out_4171406999558043376[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2019667365437356819) {
   out_2019667365437356819[0] = 0;
   out_2019667365437356819[1] = 0;
   out_2019667365437356819[2] = 0;
   out_2019667365437356819[3] = 0;
   out_2019667365437356819[4] = 1;
   out_2019667365437356819[5] = 0;
   out_2019667365437356819[6] = 0;
   out_2019667365437356819[7] = 0;
   out_2019667365437356819[8] = 0;
}
void h_26(double *state, double *unused, double *out_2227634517855805770) {
   out_2227634517855805770[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5631831737168172973) {
   out_5631831737168172973[0] = 0;
   out_5631831737168172973[1] = 0;
   out_5631831737168172973[2] = 0;
   out_5631831737168172973[3] = 0;
   out_5631831737168172973[4] = 0;
   out_5631831737168172973[5] = 0;
   out_5631831737168172973[6] = 0;
   out_5631831737168172973[7] = 1;
   out_5631831737168172973[8] = 0;
}
void h_27(double *state, double *unused, double *out_1431684425716165823) {
   out_1431684425716165823[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4194430677237781730) {
   out_4194430677237781730[0] = 0;
   out_4194430677237781730[1] = 0;
   out_4194430677237781730[2] = 0;
   out_4194430677237781730[3] = 1;
   out_4194430677237781730[4] = 0;
   out_4194430677237781730[5] = 0;
   out_4194430677237781730[6] = 0;
   out_4194430677237781730[7] = 0;
   out_4194430677237781730[8] = 0;
}
void h_29(double *state, double *unused, double *out_5339150800634185113) {
   out_5339150800634185113[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1509436021122964635) {
   out_1509436021122964635[0] = 0;
   out_1509436021122964635[1] = 1;
   out_1509436021122964635[2] = 0;
   out_1509436021122964635[3] = 0;
   out_1509436021122964635[4] = 0;
   out_1509436021122964635[5] = 0;
   out_1509436021122964635[6] = 0;
   out_1509436021122964635[7] = 0;
   out_1509436021122964635[8] = 0;
}
void h_28(double *state, double *unused, double *out_7994516235552909351) {
   out_7994516235552909351[0] = state[0];
}
void H_28(double *state, double *unused, double *out_7456551652532688279) {
   out_7456551652532688279[0] = 1;
   out_7456551652532688279[1] = 0;
   out_7456551652532688279[2] = 0;
   out_7456551652532688279[3] = 0;
   out_7456551652532688279[4] = 0;
   out_7456551652532688279[5] = 0;
   out_7456551652532688279[6] = 0;
   out_7456551652532688279[7] = 0;
   out_7456551652532688279[8] = 0;
}
void h_31(double *state, double *unused, double *out_4193276414875992619) {
   out_4193276414875992619[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8905711745052013146) {
   out_8905711745052013146[0] = 0;
   out_8905711745052013146[1] = 0;
   out_8905711745052013146[2] = 0;
   out_8905711745052013146[3] = 0;
   out_8905711745052013146[4] = 0;
   out_8905711745052013146[5] = 0;
   out_8905711745052013146[6] = 0;
   out_8905711745052013146[7] = 0;
   out_8905711745052013146[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6892551631316742231) {
  err_fun(nom_x, delta_x, out_6892551631316742231);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_632344986528208149) {
  inv_err_fun(nom_x, true_x, out_632344986528208149);
}
void car_H_mod_fun(double *state, double *out_6681531782435372252) {
  H_mod_fun(state, out_6681531782435372252);
}
void car_f_fun(double *state, double dt, double *out_5873591982632034770) {
  f_fun(state,  dt, out_5873591982632034770);
}
void car_F_fun(double *state, double dt, double *out_5222924237119091959) {
  F_fun(state,  dt, out_5222924237119091959);
}
void car_h_25(double *state, double *unused, double *out_7830725943130014969) {
  h_25(state, unused, out_7830725943130014969);
}
void car_H_25(double *state, double *unused, double *out_8936357706928973574) {
  H_25(state, unused, out_8936357706928973574);
}
void car_h_24(double *state, double *unused, double *out_7215733804563166273) {
  h_24(state, unused, out_7215733804563166273);
}
void car_H_24(double *state, double *unused, double *out_7337736767775078476) {
  H_24(state, unused, out_7337736767775078476);
}
void car_h_30(double *state, double *unused, double *out_4171406999558043376) {
  h_30(state, unused, out_4171406999558043376);
}
void car_H_30(double *state, double *unused, double *out_2019667365437356819) {
  H_30(state, unused, out_2019667365437356819);
}
void car_h_26(double *state, double *unused, double *out_2227634517855805770) {
  h_26(state, unused, out_2227634517855805770);
}
void car_H_26(double *state, double *unused, double *out_5631831737168172973) {
  H_26(state, unused, out_5631831737168172973);
}
void car_h_27(double *state, double *unused, double *out_1431684425716165823) {
  h_27(state, unused, out_1431684425716165823);
}
void car_H_27(double *state, double *unused, double *out_4194430677237781730) {
  H_27(state, unused, out_4194430677237781730);
}
void car_h_29(double *state, double *unused, double *out_5339150800634185113) {
  h_29(state, unused, out_5339150800634185113);
}
void car_H_29(double *state, double *unused, double *out_1509436021122964635) {
  H_29(state, unused, out_1509436021122964635);
}
void car_h_28(double *state, double *unused, double *out_7994516235552909351) {
  h_28(state, unused, out_7994516235552909351);
}
void car_H_28(double *state, double *unused, double *out_7456551652532688279) {
  H_28(state, unused, out_7456551652532688279);
}
void car_h_31(double *state, double *unused, double *out_4193276414875992619) {
  h_31(state, unused, out_4193276414875992619);
}
void car_H_31(double *state, double *unused, double *out_8905711745052013146) {
  H_31(state, unused, out_8905711745052013146);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
