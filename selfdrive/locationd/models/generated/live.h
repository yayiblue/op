#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_5893259616952359995);
void live_err_fun(double *nom_x, double *delta_x, double *out_4830021671970159035);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3922488873627355783);
void live_H_mod_fun(double *state, double *out_5660387373054586398);
void live_f_fun(double *state, double dt, double *out_8710621012237093742);
void live_F_fun(double *state, double dt, double *out_7370016201918838578);
void live_h_4(double *state, double *unused, double *out_1965350356141844935);
void live_H_4(double *state, double *unused, double *out_5757122479151917060);
void live_h_9(double *state, double *unused, double *out_8102558235427240065);
void live_H_9(double *state, double *unused, double *out_1530096456112530410);
void live_h_10(double *state, double *unused, double *out_5367289565694616813);
void live_H_10(double *state, double *unused, double *out_492949424303976863);
void live_h_12(double *state, double *unused, double *out_4164990670209354339);
void live_H_12(double *state, double *unused, double *out_6308363217514901560);
void live_h_31(double *state, double *unused, double *out_5706607973709820228);
void live_H_31(double *state, double *unused, double *out_2007896961205058444);
void live_h_32(double *state, double *unused, double *out_7762004935965751907);
void live_H_32(double *state, double *unused, double *out_6609159659451218754);
void live_h_13(double *state, double *unused, double *out_5405810379734963908);
void live_H_13(double *state, double *unused, double *out_6535263366373938463);
void live_h_14(double *state, double *unused, double *out_8102558235427240065);
void live_H_14(double *state, double *unused, double *out_1530096456112530410);
void live_h_33(double *state, double *unused, double *out_2591779336706256217);
void live_H_33(double *state, double *unused, double *out_760096582859547920);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}