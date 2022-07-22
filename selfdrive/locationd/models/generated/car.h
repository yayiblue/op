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
void car_err_fun(double *nom_x, double *delta_x, double *out_6013722011862595432);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_825540609892149147);
void car_H_mod_fun(double *state, double *out_3795912724315969018);
void car_f_fun(double *state, double dt, double *out_6140302451994653393);
void car_F_fun(double *state, double dt, double *out_8479765871903219056);
void car_h_25(double *state, double *unused, double *out_6776265007677880447);
void car_H_25(double *state, double *unused, double *out_7428847556636375798);
void car_h_24(double *state, double *unused, double *out_4941408709385696404);
void car_H_24(double *state, double *unused, double *out_2608526051980387535);
void car_h_30(double *state, double *unused, double *out_8033029537777648819);
void car_H_30(double *state, double *unused, double *out_7299508609493135728);
void car_h_26(double *state, double *unused, double *out_5347297017213560074);
void car_H_26(double *state, double *unused, double *out_3687344237762319574);
void car_h_27(double *state, double *unused, double *out_6668805798370825524);
void car_H_27(double *state, double *unused, double *out_5124745297692710817);
void car_h_29(double *state, double *unused, double *out_3031356270116803174);
void car_H_29(double *state, double *unused, double *out_7809739953807527912);
void car_h_28(double *state, double *unused, double *out_2049031121497866802);
void car_H_28(double *state, double *unused, double *out_2727340936737997338);
void car_h_31(double *state, double *unused, double *out_6501070945393374558);
void car_H_31(double *state, double *unused, double *out_7459493518513336226);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}