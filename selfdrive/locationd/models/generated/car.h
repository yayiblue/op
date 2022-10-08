#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_6892551631316742231);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_632344986528208149);
void car_H_mod_fun(double *state, double *out_6681531782435372252);
void car_f_fun(double *state, double dt, double *out_5873591982632034770);
void car_F_fun(double *state, double dt, double *out_5222924237119091959);
void car_h_25(double *state, double *unused, double *out_7830725943130014969);
void car_H_25(double *state, double *unused, double *out_8936357706928973574);
void car_h_24(double *state, double *unused, double *out_7215733804563166273);
void car_H_24(double *state, double *unused, double *out_7337736767775078476);
void car_h_30(double *state, double *unused, double *out_4171406999558043376);
void car_H_30(double *state, double *unused, double *out_2019667365437356819);
void car_h_26(double *state, double *unused, double *out_2227634517855805770);
void car_H_26(double *state, double *unused, double *out_5631831737168172973);
void car_h_27(double *state, double *unused, double *out_1431684425716165823);
void car_H_27(double *state, double *unused, double *out_4194430677237781730);
void car_h_29(double *state, double *unused, double *out_5339150800634185113);
void car_H_29(double *state, double *unused, double *out_1509436021122964635);
void car_h_28(double *state, double *unused, double *out_7994516235552909351);
void car_H_28(double *state, double *unused, double *out_7456551652532688279);
void car_h_31(double *state, double *unused, double *out_4193276414875992619);
void car_H_31(double *state, double *unused, double *out_8905711745052013146);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}