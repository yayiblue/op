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
void err_fun(double *nom_x, double *delta_x, double *out_6013722011862595432) {
   out_6013722011862595432[0] = delta_x[0] + nom_x[0];
   out_6013722011862595432[1] = delta_x[1] + nom_x[1];
   out_6013722011862595432[2] = delta_x[2] + nom_x[2];
   out_6013722011862595432[3] = delta_x[3] + nom_x[3];
   out_6013722011862595432[4] = delta_x[4] + nom_x[4];
   out_6013722011862595432[5] = delta_x[5] + nom_x[5];
   out_6013722011862595432[6] = delta_x[6] + nom_x[6];
   out_6013722011862595432[7] = delta_x[7] + nom_x[7];
   out_6013722011862595432[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_825540609892149147) {
   out_825540609892149147[0] = -nom_x[0] + true_x[0];
   out_825540609892149147[1] = -nom_x[1] + true_x[1];
   out_825540609892149147[2] = -nom_x[2] + true_x[2];
   out_825540609892149147[3] = -nom_x[3] + true_x[3];
   out_825540609892149147[4] = -nom_x[4] + true_x[4];
   out_825540609892149147[5] = -nom_x[5] + true_x[5];
   out_825540609892149147[6] = -nom_x[6] + true_x[6];
   out_825540609892149147[7] = -nom_x[7] + true_x[7];
   out_825540609892149147[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3795912724315969018) {
   out_3795912724315969018[0] = 1.0;
   out_3795912724315969018[1] = 0;
   out_3795912724315969018[2] = 0;
   out_3795912724315969018[3] = 0;
   out_3795912724315969018[4] = 0;
   out_3795912724315969018[5] = 0;
   out_3795912724315969018[6] = 0;
   out_3795912724315969018[7] = 0;
   out_3795912724315969018[8] = 0;
   out_3795912724315969018[9] = 0;
   out_3795912724315969018[10] = 1.0;
   out_3795912724315969018[11] = 0;
   out_3795912724315969018[12] = 0;
   out_3795912724315969018[13] = 0;
   out_3795912724315969018[14] = 0;
   out_3795912724315969018[15] = 0;
   out_3795912724315969018[16] = 0;
   out_3795912724315969018[17] = 0;
   out_3795912724315969018[18] = 0;
   out_3795912724315969018[19] = 0;
   out_3795912724315969018[20] = 1.0;
   out_3795912724315969018[21] = 0;
   out_3795912724315969018[22] = 0;
   out_3795912724315969018[23] = 0;
   out_3795912724315969018[24] = 0;
   out_3795912724315969018[25] = 0;
   out_3795912724315969018[26] = 0;
   out_3795912724315969018[27] = 0;
   out_3795912724315969018[28] = 0;
   out_3795912724315969018[29] = 0;
   out_3795912724315969018[30] = 1.0;
   out_3795912724315969018[31] = 0;
   out_3795912724315969018[32] = 0;
   out_3795912724315969018[33] = 0;
   out_3795912724315969018[34] = 0;
   out_3795912724315969018[35] = 0;
   out_3795912724315969018[36] = 0;
   out_3795912724315969018[37] = 0;
   out_3795912724315969018[38] = 0;
   out_3795912724315969018[39] = 0;
   out_3795912724315969018[40] = 1.0;
   out_3795912724315969018[41] = 0;
   out_3795912724315969018[42] = 0;
   out_3795912724315969018[43] = 0;
   out_3795912724315969018[44] = 0;
   out_3795912724315969018[45] = 0;
   out_3795912724315969018[46] = 0;
   out_3795912724315969018[47] = 0;
   out_3795912724315969018[48] = 0;
   out_3795912724315969018[49] = 0;
   out_3795912724315969018[50] = 1.0;
   out_3795912724315969018[51] = 0;
   out_3795912724315969018[52] = 0;
   out_3795912724315969018[53] = 0;
   out_3795912724315969018[54] = 0;
   out_3795912724315969018[55] = 0;
   out_3795912724315969018[56] = 0;
   out_3795912724315969018[57] = 0;
   out_3795912724315969018[58] = 0;
   out_3795912724315969018[59] = 0;
   out_3795912724315969018[60] = 1.0;
   out_3795912724315969018[61] = 0;
   out_3795912724315969018[62] = 0;
   out_3795912724315969018[63] = 0;
   out_3795912724315969018[64] = 0;
   out_3795912724315969018[65] = 0;
   out_3795912724315969018[66] = 0;
   out_3795912724315969018[67] = 0;
   out_3795912724315969018[68] = 0;
   out_3795912724315969018[69] = 0;
   out_3795912724315969018[70] = 1.0;
   out_3795912724315969018[71] = 0;
   out_3795912724315969018[72] = 0;
   out_3795912724315969018[73] = 0;
   out_3795912724315969018[74] = 0;
   out_3795912724315969018[75] = 0;
   out_3795912724315969018[76] = 0;
   out_3795912724315969018[77] = 0;
   out_3795912724315969018[78] = 0;
   out_3795912724315969018[79] = 0;
   out_3795912724315969018[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6140302451994653393) {
   out_6140302451994653393[0] = state[0];
   out_6140302451994653393[1] = state[1];
   out_6140302451994653393[2] = state[2];
   out_6140302451994653393[3] = state[3];
   out_6140302451994653393[4] = state[4];
   out_6140302451994653393[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6140302451994653393[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6140302451994653393[7] = state[7];
   out_6140302451994653393[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8479765871903219056) {
   out_8479765871903219056[0] = 1;
   out_8479765871903219056[1] = 0;
   out_8479765871903219056[2] = 0;
   out_8479765871903219056[3] = 0;
   out_8479765871903219056[4] = 0;
   out_8479765871903219056[5] = 0;
   out_8479765871903219056[6] = 0;
   out_8479765871903219056[7] = 0;
   out_8479765871903219056[8] = 0;
   out_8479765871903219056[9] = 0;
   out_8479765871903219056[10] = 1;
   out_8479765871903219056[11] = 0;
   out_8479765871903219056[12] = 0;
   out_8479765871903219056[13] = 0;
   out_8479765871903219056[14] = 0;
   out_8479765871903219056[15] = 0;
   out_8479765871903219056[16] = 0;
   out_8479765871903219056[17] = 0;
   out_8479765871903219056[18] = 0;
   out_8479765871903219056[19] = 0;
   out_8479765871903219056[20] = 1;
   out_8479765871903219056[21] = 0;
   out_8479765871903219056[22] = 0;
   out_8479765871903219056[23] = 0;
   out_8479765871903219056[24] = 0;
   out_8479765871903219056[25] = 0;
   out_8479765871903219056[26] = 0;
   out_8479765871903219056[27] = 0;
   out_8479765871903219056[28] = 0;
   out_8479765871903219056[29] = 0;
   out_8479765871903219056[30] = 1;
   out_8479765871903219056[31] = 0;
   out_8479765871903219056[32] = 0;
   out_8479765871903219056[33] = 0;
   out_8479765871903219056[34] = 0;
   out_8479765871903219056[35] = 0;
   out_8479765871903219056[36] = 0;
   out_8479765871903219056[37] = 0;
   out_8479765871903219056[38] = 0;
   out_8479765871903219056[39] = 0;
   out_8479765871903219056[40] = 1;
   out_8479765871903219056[41] = 0;
   out_8479765871903219056[42] = 0;
   out_8479765871903219056[43] = 0;
   out_8479765871903219056[44] = 0;
   out_8479765871903219056[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8479765871903219056[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8479765871903219056[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8479765871903219056[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8479765871903219056[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8479765871903219056[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8479765871903219056[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8479765871903219056[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8479765871903219056[53] = -9.8000000000000007*dt;
   out_8479765871903219056[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8479765871903219056[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8479765871903219056[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8479765871903219056[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8479765871903219056[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8479765871903219056[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8479765871903219056[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8479765871903219056[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8479765871903219056[62] = 0;
   out_8479765871903219056[63] = 0;
   out_8479765871903219056[64] = 0;
   out_8479765871903219056[65] = 0;
   out_8479765871903219056[66] = 0;
   out_8479765871903219056[67] = 0;
   out_8479765871903219056[68] = 0;
   out_8479765871903219056[69] = 0;
   out_8479765871903219056[70] = 1;
   out_8479765871903219056[71] = 0;
   out_8479765871903219056[72] = 0;
   out_8479765871903219056[73] = 0;
   out_8479765871903219056[74] = 0;
   out_8479765871903219056[75] = 0;
   out_8479765871903219056[76] = 0;
   out_8479765871903219056[77] = 0;
   out_8479765871903219056[78] = 0;
   out_8479765871903219056[79] = 0;
   out_8479765871903219056[80] = 1;
}
void h_25(double *state, double *unused, double *out_6776265007677880447) {
   out_6776265007677880447[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7428847556636375798) {
   out_7428847556636375798[0] = 0;
   out_7428847556636375798[1] = 0;
   out_7428847556636375798[2] = 0;
   out_7428847556636375798[3] = 0;
   out_7428847556636375798[4] = 0;
   out_7428847556636375798[5] = 0;
   out_7428847556636375798[6] = 1;
   out_7428847556636375798[7] = 0;
   out_7428847556636375798[8] = 0;
}
void h_24(double *state, double *unused, double *out_4941408709385696404) {
   out_4941408709385696404[0] = state[4];
   out_4941408709385696404[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2608526051980387535) {
   out_2608526051980387535[0] = 0;
   out_2608526051980387535[1] = 0;
   out_2608526051980387535[2] = 0;
   out_2608526051980387535[3] = 0;
   out_2608526051980387535[4] = 1;
   out_2608526051980387535[5] = 0;
   out_2608526051980387535[6] = 0;
   out_2608526051980387535[7] = 0;
   out_2608526051980387535[8] = 0;
   out_2608526051980387535[9] = 0;
   out_2608526051980387535[10] = 0;
   out_2608526051980387535[11] = 0;
   out_2608526051980387535[12] = 0;
   out_2608526051980387535[13] = 0;
   out_2608526051980387535[14] = 1;
   out_2608526051980387535[15] = 0;
   out_2608526051980387535[16] = 0;
   out_2608526051980387535[17] = 0;
}
void h_30(double *state, double *unused, double *out_8033029537777648819) {
   out_8033029537777648819[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7299508609493135728) {
   out_7299508609493135728[0] = 0;
   out_7299508609493135728[1] = 0;
   out_7299508609493135728[2] = 0;
   out_7299508609493135728[3] = 0;
   out_7299508609493135728[4] = 1;
   out_7299508609493135728[5] = 0;
   out_7299508609493135728[6] = 0;
   out_7299508609493135728[7] = 0;
   out_7299508609493135728[8] = 0;
}
void h_26(double *state, double *unused, double *out_5347297017213560074) {
   out_5347297017213560074[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3687344237762319574) {
   out_3687344237762319574[0] = 0;
   out_3687344237762319574[1] = 0;
   out_3687344237762319574[2] = 0;
   out_3687344237762319574[3] = 0;
   out_3687344237762319574[4] = 0;
   out_3687344237762319574[5] = 0;
   out_3687344237762319574[6] = 0;
   out_3687344237762319574[7] = 1;
   out_3687344237762319574[8] = 0;
}
void h_27(double *state, double *unused, double *out_6668805798370825524) {
   out_6668805798370825524[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5124745297692710817) {
   out_5124745297692710817[0] = 0;
   out_5124745297692710817[1] = 0;
   out_5124745297692710817[2] = 0;
   out_5124745297692710817[3] = 1;
   out_5124745297692710817[4] = 0;
   out_5124745297692710817[5] = 0;
   out_5124745297692710817[6] = 0;
   out_5124745297692710817[7] = 0;
   out_5124745297692710817[8] = 0;
}
void h_29(double *state, double *unused, double *out_3031356270116803174) {
   out_3031356270116803174[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7809739953807527912) {
   out_7809739953807527912[0] = 0;
   out_7809739953807527912[1] = 1;
   out_7809739953807527912[2] = 0;
   out_7809739953807527912[3] = 0;
   out_7809739953807527912[4] = 0;
   out_7809739953807527912[5] = 0;
   out_7809739953807527912[6] = 0;
   out_7809739953807527912[7] = 0;
   out_7809739953807527912[8] = 0;
}
void h_28(double *state, double *unused, double *out_2049031121497866802) {
   out_2049031121497866802[0] = state[0];
}
void H_28(double *state, double *unused, double *out_2727340936737997338) {
   out_2727340936737997338[0] = 1;
   out_2727340936737997338[1] = 0;
   out_2727340936737997338[2] = 0;
   out_2727340936737997338[3] = 0;
   out_2727340936737997338[4] = 0;
   out_2727340936737997338[5] = 0;
   out_2727340936737997338[6] = 0;
   out_2727340936737997338[7] = 0;
   out_2727340936737997338[8] = 0;
}
void h_31(double *state, double *unused, double *out_6501070945393374558) {
   out_6501070945393374558[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7459493518513336226) {
   out_7459493518513336226[0] = 0;
   out_7459493518513336226[1] = 0;
   out_7459493518513336226[2] = 0;
   out_7459493518513336226[3] = 0;
   out_7459493518513336226[4] = 0;
   out_7459493518513336226[5] = 0;
   out_7459493518513336226[6] = 0;
   out_7459493518513336226[7] = 0;
   out_7459493518513336226[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_6013722011862595432) {
  err_fun(nom_x, delta_x, out_6013722011862595432);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_825540609892149147) {
  inv_err_fun(nom_x, true_x, out_825540609892149147);
}
void car_H_mod_fun(double *state, double *out_3795912724315969018) {
  H_mod_fun(state, out_3795912724315969018);
}
void car_f_fun(double *state, double dt, double *out_6140302451994653393) {
  f_fun(state,  dt, out_6140302451994653393);
}
void car_F_fun(double *state, double dt, double *out_8479765871903219056) {
  F_fun(state,  dt, out_8479765871903219056);
}
void car_h_25(double *state, double *unused, double *out_6776265007677880447) {
  h_25(state, unused, out_6776265007677880447);
}
void car_H_25(double *state, double *unused, double *out_7428847556636375798) {
  H_25(state, unused, out_7428847556636375798);
}
void car_h_24(double *state, double *unused, double *out_4941408709385696404) {
  h_24(state, unused, out_4941408709385696404);
}
void car_H_24(double *state, double *unused, double *out_2608526051980387535) {
  H_24(state, unused, out_2608526051980387535);
}
void car_h_30(double *state, double *unused, double *out_8033029537777648819) {
  h_30(state, unused, out_8033029537777648819);
}
void car_H_30(double *state, double *unused, double *out_7299508609493135728) {
  H_30(state, unused, out_7299508609493135728);
}
void car_h_26(double *state, double *unused, double *out_5347297017213560074) {
  h_26(state, unused, out_5347297017213560074);
}
void car_H_26(double *state, double *unused, double *out_3687344237762319574) {
  H_26(state, unused, out_3687344237762319574);
}
void car_h_27(double *state, double *unused, double *out_6668805798370825524) {
  h_27(state, unused, out_6668805798370825524);
}
void car_H_27(double *state, double *unused, double *out_5124745297692710817) {
  H_27(state, unused, out_5124745297692710817);
}
void car_h_29(double *state, double *unused, double *out_3031356270116803174) {
  h_29(state, unused, out_3031356270116803174);
}
void car_H_29(double *state, double *unused, double *out_7809739953807527912) {
  H_29(state, unused, out_7809739953807527912);
}
void car_h_28(double *state, double *unused, double *out_2049031121497866802) {
  h_28(state, unused, out_2049031121497866802);
}
void car_H_28(double *state, double *unused, double *out_2727340936737997338) {
  H_28(state, unused, out_2727340936737997338);
}
void car_h_31(double *state, double *unused, double *out_6501070945393374558) {
  h_31(state, unused, out_6501070945393374558);
}
void car_H_31(double *state, double *unused, double *out_7459493518513336226) {
  H_31(state, unused, out_7459493518513336226);
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
