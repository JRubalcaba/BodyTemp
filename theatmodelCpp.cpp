#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Q(double Tb,
                double A1,
                double a,
                double S,
                double EWL_estimate,
                double Ta,
                double hc,
                double Ag,
                double hg,
                double Tg){
  
  double Tb_t = Tb + 273;
  double Ta_t = Ta + 273;
  double Tg_t = Tg + 273;
  
  double qabs = A1 * abs * S;
  double qevap =  2257 * EWL; 
  double qrad = 5.670367e-8 * emis * (A1 * (Tb_t*Tb_t*Tb_t*Tb_t - Ta_t*Ta_t*Ta_t*Ta_t) + Ag * (Tb_t*Tb_t*Tb_t*Tb_t - Tg_t*Tg_t*Tg_t*Tg_t));  
  double qconv = A1 * hc * (Tb_t - Ta_t);                   
  double qcond = Ag * hg * (Tb_t - Tg_t);                   
  double qnet = qabs - qevap - qrad - qconv - qcond;    
  
  return qnet;
}


// [[Rcpp::export]]
NumericVector theatmodelCpp(double Tb,
                            double A,        
                            double M,        
                            double Ta,      
                            double Tg,      
                            double S,        
                            double v,       
                            double HR,      
                            double skin_humidity = 1,
                            double resistance = 60000,        
                            double posture = 1, 
                            double vent = 2/3, 
                            int delta = 60, 
                            double C = 3.6){
  NumericVector out(2);
  for(int i = 0; i < delta; i++){
    // Model parameters
    //// Animal geometry and color
    double Ag = vent * A;
    double A1 = (1-vent) * A * posture; 
    double l = pow(A, 1/2);
    double a = 0.8;   
    
    //// Estimation convection coefficient
    double nu = -1.1555e-14*pow(Ta+273, 3) + 9.5728e-11*pow(Ta+273, 2) + 3.7604e-08*(Ta+273) - 3.4484e-06;   
    double kf = 1.5207e-11*pow(Ta+273, 3) - 4.8574e-08*pow(Ta+273, 2) + 1.0184e-04*(Ta+273) - 3.9333e-04; 
   
    double Re = l * v / nu;   
    double c = 0.37;          
    double n = 0.7; 
    double Nu = 2 + c * pow(Re,n); 
      
    double hc = Nu * kf / l;

    //// Estimation of conduction coefficient
    double ksub = 0.027;       
    double ts = 0.025 * pow(0.001 * M / (3.1416 * 1000), 0.2);   
    double hg = ksub / ts; 
    
    //// Estimate evaporative water loss (g s-1)
    double HR_prop = HR * 0.01;
    double ps_a = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / pow(Ta+273, 8.2); 
    double rho_a = HR_prop * 2.2 * ps_a / (Ta+273);   
      
    double ps_b = exp(77.3450 + 0.0057 * (Tb+273) - 7235 / (Tb+273)) / pow(Tb+273, 8.2); 
    double rho_b = skin_humidity * 2.2 * ps_b / (Tb+273);
      
    double EWL_estimate = A1 * (rho_b - rho_a) / resistance; 
    
    // Model iteration
    Tb = Tb + 1 / (C * M) * Q(Tb, A1, a, S, EWL_estimate, Ta, hc, Ag, hg, Tg);
    
    // OUTPUT
    out[0] = Tb;
    out[1] = EWL_estimate;
  }
  return out;
}
