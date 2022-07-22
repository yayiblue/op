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
void live_H(double *in_vec, double *out_1986321194722187952);
void live_err_fun(double *nom_x, double *delta_x, double *out_9060234064802883238);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6701773287529323107);
void live_H_mod_fun(double *state, double *out_771851589417143562);
void live_f_fun(double *state, double dt, double *out_8469395992817321048);
void live_F_fun(double *state, double dt, double *out_9007411877149068418);
void live_h_4(double *state, double *unused, double *out_4496058279147309907);
void live_H_4(double *state, double *unused, double *out_3144930195240663846);
void live_h_9(double *state, double *unused, double *out_2177212884378738909);
void live_H_9(double *state, double *unused, double *out_4142288740023783624);
void live_h_10(double *state, double *unused, double *out_1594179005802415580);
void live_H_10(double *state, double *unused, double *out_8017581178918911948);
void live_h_12(double *state, double *unused, double *out_1529497512963034638);
void live_H_12(double *state, double *unused, double *out_4522198118441786646);
void live_h_31(double *state, double *unused, double *out_6260980419806761527);
void live_H_31(double *state, double *unused, double *out_221731862131943530);
void live_h_32(double *state, double *unused, double *out_8562679850602801047);
void live_H_32(double *state, double *unused, double *out_7631138582955592532);
void live_h_13(double *state, double *unused, double *out_6633941345833725407);
void live_H_13(double *state, double *unused, double *out_2333890705352792507);
void live_h_14(double *state, double *unused, double *out_2177212884378738909);
void live_H_14(double *state, double *unused, double *out_4142288740023783624);
void live_h_33(double *state, double *unused, double *out_8680022473253276378);
void live_H_33(double *state, double *unused, double *out_8028425918303893657);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}