#include "../../include/common/common.h"
#include "../../include/generators/power_calc.h"

double computeElectricalPower(double Ed2, double id, double Eq2, double iq, double Xd2, double Xq2) {       
    double power_elec;
    double term1 = Ed2*id;
    double term2 = Eq2*iq;
    double term3 = (Xd2-Xq2)*(id)*(iq);
    
    power_elec = term1+term2+term3;
    return power_elec;
} 