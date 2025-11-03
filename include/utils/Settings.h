#ifndef SETTINGS_H
#define SETTINGS_H

#include <optional>

namespace mol_sim {

class SettingsParam {
   public:
    double delta_t;
    double start_time;
    double end_time;
    SettingsParam(double delta_t = 0.014, double start_time = 0., double end_time = 1000.0)
        : delta_t(delta_t), start_time(start_time), end_time(end_time) {}
};
}  // namespace mol_sim
#endif
