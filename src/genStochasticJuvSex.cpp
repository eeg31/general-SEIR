#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector stepModel(IntegerVector state, NumericVector r, IntegerVector spec, int nevents)
{
  IntegerVector var = clone(state);
  
  int* whichBeta; whichBeta = &spec[0];   
  int* whichRho; whichRho = &spec[1];     
  int* whichGamma; whichGamma = &spec[2];
  
  int* SF; SF = &var[0];
  int* EF; EF = &var[1];
  int* IF; IF = &var[2];
  int* RF; RF = &var[3];
  int* SJ0; SJ0 = &var[4];
  int* SJ; SJ = &var[5];
  int* MJ; MJ = &var[6];
  int* EJ; EJ = &var[7];
  int* IJ; IJ = &var[8];
  int* RJ; RJ = &var[9];
  int* EJ0; EJ0 = &var[10];
  int* IJ0; IJ0 = &var[11];
  int* RJ0; RJ0 = &var[12];
  int* SM; SM = &var[13];
  int* EM; EM = &var[14];
  int* IM; IM = &var[15];
  int* RM; RM = &var[16];
  
  double totalRate = sum(r);
  double x = R::runif(0,1);
  x *= totalRate;
  
  int mindex;
  bool isDone = false;
  for (int i = 0; !isDone; ++i){
    double sum = 0;
    for (int j = 0; j <= i; ++j){
      sum += r[j];
    }
    if (sum > x) {
      mindex = i;
      isDone = true;
    } 
  }
  
  NumericVector nEvents(nevents);
  for (int i = 0; i < nevents; i++) {
    if (i == mindex) {
      nEvents[i] = 1;
    } else {
      nEvents[i] = 0;                
    }
  }
  
  *SF   = *SF - nEvents[1] - nEvents[5] + nEvents[10] + nEvents[26];
  *EF   = *EF - nEvents[2] - nEvents[6] - nEvents[7] + nEvents[9] + nEvents[27];
  *IF   = *IF - nEvents[3] + nEvents[7] - nEvents[8] - nEvents[9] + nEvents[28];
  *RF   = *RF - nEvents[4] - nEvents[10] + nEvents[29];
  *SJ0 = *SJ0 + nEvents[0] - nEvents[12] - nEvents[21] - nEvents[30] + nEvents[40];
  *SJ  = *SJ - nEvents[13] + nEvents[21] + nEvents[22] - nEvents[26] - nEvents[31] + nEvents[41] - nEvents[52];
  *MJ  = *MJ + nEvents[11] - nEvents[14] - nEvents[22];
  *EJ0 = *EJ0 - nEvents[15] - nEvents[23] - nEvents[32] - nEvents[34] + nEvents[38];
  *EJ  = *EJ - nEvents[16] + nEvents[23] - nEvents[27] - nEvents[33] - nEvents[35] + nEvents[39] - nEvents[53];
  *IJ0 = *IJ0 - nEvents[17] - nEvents[24] + nEvents[34] - nEvents[36] - nEvents[38];
  *IJ  = *IJ - nEvents[18] + nEvents[24] - nEvents[28] + nEvents[35] - nEvents[37] - nEvents[39] - nEvents[54];
  *RJ0 = *RJ0 - nEvents[19] - nEvents[25] - nEvents[40];
  *RJ  = *RJ - nEvents[20] + nEvents[25] - nEvents[29] - nEvents[41] - nEvents[55];
  *SM   = *SM - nEvents[42] - nEvents[46] + nEvents[51] + nEvents[52];
  *EM   = *EM - nEvents[43] - nEvents[47] - nEvents[48] + nEvents[50] + nEvents[53];
  *IM   = *IM - nEvents[44] + nEvents[48] - nEvents[49] - nEvents[50] + nEvents[54];
  *RM   = *RM - nEvents[45] - nEvents[51] + nEvents[55];
  
  if (*whichBeta == 1) {
    *EF  = *EF + nEvents[5];
    *EM  = *EM + nEvents[46];
    *EJ0 = *EJ0 + nEvents[30];
    *EJ  = *EJ + nEvents[31];
  } else {
    *IF  = *IF + nEvents[5];
    *IM  = *IM + nEvents[46];
    *IJ0 = *IJ0 + nEvents[30];
    *IJ  = *IJ + nEvents[31];
  }
  
  if (*whichRho == 1) {
    *SF  = *SF + nEvents[6];
    *SM  = *SM + nEvents[47];
    *SJ0 = *SJ0 + nEvents[32];
    *SJ  = *SJ + nEvents[33];
  } else {
    *RF  = *RF + nEvents[6];
    *RM  = *RM + nEvents[47];
    *RJ0 = *RJ0 + nEvents[32];
    *RJ  = *RJ + nEvents[33];
  }
  
  if (*whichGamma == 1){
    *SF  = *SF + nEvents[8];
    *SM  = *SM + nEvents[49];
    *SJ0 = *SJ0 + nEvents[36]; 
    *SJ  = *SJ + nEvents[37];
  } else {
    *RF  = *RF + nEvents[8];
    *RM  = *RM + nEvents[49];
    *RJ0 = *RJ0 + nEvents[36];
    *RJ  = *RJ + nEvents[37];
  }
  
  return var;
}



//[[Rcpp::export]]
IntegerVector tauleap(IntegerVector state, NumericVector r, IntegerVector spec, int nevents, double tau)
{
  IntegerVector var = clone(state);
  
  int* whichBeta; whichBeta = &spec[0];   
  int* whichRho; whichRho = &spec[1];     
  int* whichGamma; whichGamma = &spec[2];
  
  int* SF; SF = &var[0];
  int* EF; EF = &var[1];
  int* IF; IF = &var[2];
  int* RF; RF = &var[3];
  int* SJ0; SJ0 = &var[4];
  int* SJ; SJ = &var[5];
  int* MJ; MJ = &var[6];
  int* EJ; EJ = &var[7];
  int* IJ; IJ = &var[8];
  int* RJ; RJ = &var[9];
  int* EJ0; EJ0 = &var[10];
  int* IJ0; IJ0 = &var[11];
  int* RJ0; RJ0 = &var[12];
  int* SM; SM = &var[13];
  int* EM; EM = &var[14];
  int* IM; IM = &var[15];
  int* RM; RM = &var[16];
  
  NumericVector nEvents(nevents);
  for (int i = 0; i < nevents; i++) {
    nEvents[i] = R::rpois(r[i] * tau);
  }
  
  *SF   = *SF - nEvents[1] - nEvents[5] + nEvents[10] + nEvents[26];
  *EF   = *EF - nEvents[2] - nEvents[6] - nEvents[7] + nEvents[9] + nEvents[27];
  *IF   = *IF - nEvents[3] + nEvents[7] - nEvents[8] - nEvents[9] + nEvents[28];
  *RF   = *RF - nEvents[4] - nEvents[10] + nEvents[29];
  *SJ0 = *SJ0 + nEvents[0] - nEvents[12] - nEvents[21] - nEvents[30] + nEvents[40];
  *SJ  = *SJ - nEvents[13] + nEvents[21] + nEvents[22] - nEvents[26] - nEvents[31] + nEvents[41] - nEvents[52];
  *MJ  = *MJ + nEvents[11] - nEvents[14] - nEvents[22];
  *EJ0 = *EJ0 - nEvents[15] - nEvents[23] - nEvents[32] - nEvents[34] + nEvents[38];
  *EJ  = *EJ - nEvents[16] + nEvents[23] - nEvents[27] - nEvents[33] - nEvents[35] + nEvents[39] - nEvents[53];
  *IJ0 = *IJ0 - nEvents[17] - nEvents[24] + nEvents[34] - nEvents[36] - nEvents[38];
  *IJ  = *IJ - nEvents[18] + nEvents[24] - nEvents[28] + nEvents[35] - nEvents[37] - nEvents[39] - nEvents[54];
  *RJ0 = *RJ0 - nEvents[19] - nEvents[25] - nEvents[40];
  *RJ  = *RJ - nEvents[20] + nEvents[25] - nEvents[29] - nEvents[41] - nEvents[55];
  *SM   = *SM - nEvents[42] - nEvents[46] + nEvents[51] + nEvents[52];
  *EM   = *EM - nEvents[43] - nEvents[47] - nEvents[48] + nEvents[50] + nEvents[53];
  *IM   = *IM - nEvents[44] + nEvents[48] - nEvents[49] - nEvents[50] + nEvents[54];
  *RM   = *RM - nEvents[45] - nEvents[51] + nEvents[55];
  
  if (*whichBeta == 1) {
    *EF  = *EF + nEvents[5];
    *EM  = *EM + nEvents[46];
    *EJ0 = *EJ0 + nEvents[30];
    *EJ  = *EJ + nEvents[31];
  } else {
    *IF  = *IF + nEvents[5];
    *IM  = *IM + nEvents[46];
    *IJ0 = *IJ0 + nEvents[30];
    *IJ  = *IJ + nEvents[31];
  }
  
  if (*whichRho == 1) {
    *SF  = *SF + nEvents[6];
    *SM  = *SM + nEvents[47];
    *SJ0 = *SJ0 + nEvents[32];
    *SJ  = *SJ + nEvents[33];
  } else {
    *RF  = *RF + nEvents[6];
    *RM  = *RM + nEvents[47];
    *RJ0 = *RJ0 + nEvents[32];
    *RJ  = *RJ + nEvents[33];
  }
  
  if (*whichGamma == 1){
    *SF  = *SF + nEvents[8];
    *SM  = *SM + nEvents[49];
    *SJ0 = *SJ0 + nEvents[36]; 
    *SJ  = *SJ + nEvents[37];
  } else {
    *RF  = *RF + nEvents[8];
    *RM  = *RM + nEvents[49];
    *RJ0 = *RJ0 + nEvents[36];
    *RJ  = *RJ + nEvents[37];
  }
  
  
  return var;
}

#include <Rcpp.h>
#include <math.h>
//[[Rcpp::export]]
NumericVector getRates(IntegerVector state, NumericVector par, int nevents, double t)
{
  IntegerVector var = clone(state);
  
  double* beta; beta = &par[0];
  double* gamma; gamma = &par[1];
  double* rho; rho = &par[2];
  double* p; p = &par[3];
  double* epsilon; epsilon = &par[4];
  double* omega; omega = &par[5];
  double* c; c = &par[6];
  double* m; m = &par[7];
  double* omegam; omegam = &par[8];
  double* mu; mu = &par[9];
  double* mj; mj = &par[10];
  double* phi; phi = &par[11];
  double* s; s = &par[12];
  double* k; k = &par[13];
  double* preg; preg = &par[14];
  
  double b = *c * exp((0.0 - *s) * pow(cos(M_PI * t - *phi), 2.0));
  double epsilon2 = *preg * exp((0.0 - *s) * pow(cos(M_PI * t - *phi), 2.0));
    
  NumericVector rates(nevents);
  
  int* SF; SF = &var[0];
  int* EF; EF = &var[1];
  int* IF; IF = &var[2];
  int* RF; RF = &var[3];
  int* SJ0; SJ0 = &var[4];
  int* SJ; SJ = &var[5];
  int* MJ; MJ = &var[6];
  int* EJ; EJ = &var[7];
  int* IJ; IJ = &var[8];
  int* RJ; RJ = &var[9];
  int* EJ0; EJ0 = &var[10];
  int* IJ0; IJ0 = &var[11];
  int* RJ0; RJ0 = &var[12];
  int* SM; SM = &var[13];
  int* EM; EM = &var[14];
  int* IM; IM = &var[15];
  int* RM; RM = &var[16];
  
  int S = *SF + *SM;
  int E = *EF + *EM;
  int I = *IF + *IM;
  int R = *RF + *RM;
  
  int N = S + E + I;
  int Nf = *SF + *EF + *IF;
  int Nall = N + R + *SJ0 + *SJ + *MJ + *EJ + *IJ + *RJ + *EJ0 + *IJ0 + *RJ0;
  
  rates[0] = 2.0 * b * Nf;
  rates[1] = *m * *SF * Nall/ *k;
  rates[2] = *m * *EF * Nall/ *k;
  rates[3] = *m * *IF * Nall/ *k;
  rates[4] = *m * *RF * Nall/ *k;
  rates[5] = *beta * *SF * (I + *IJ + *IJ0);
  rates[6] = *rho * *EF;
  rates[7] = (*epsilon + epsilon2) * *EF;
  rates[8] = *gamma * *IF;
  rates[9] = *p * *IF;
  rates[10] = *omega * *RF;
  rates[11] = 2.0 * b * *RF;
  rates[12] = *mj * *SJ0 * Nall/ *k;
  rates[13] = *mj * *SJ * Nall/ *k;
  rates[14] = *mj * *MJ * Nall/ *k;
  rates[15] = *mj * *EJ0 * Nall/ *k;
  rates[16] = *mj * *EJ * Nall/ *k;
  rates[17] = *mj * *IJ0 * Nall/ *k;
  rates[18] = *mj * *IJ * Nall/ *k;
  rates[19] = *mj * *RJ0 * Nall/ *k;
  rates[20] = *mj * *RJ * Nall/ *k;
  rates[21] = *omegam * *SJ0;
  rates[22] = *omegam * *MJ;
  rates[23] = *omegam * *EJ0;
  rates[24] = *omegam * *IJ0;
  rates[25] = *omegam * *RJ0;
  rates[26] = *mu * *SJ / 2.0;
  rates[27] = *mu * *EJ / 2.0;
  rates[28] = *mu * *IJ / 2.0;
  rates[29] = *mu * *RJ / 2.0;
  rates[30] = *beta * *SJ0 * (I + *IJ + *IJ0);
  rates[31] = *beta * *SJ * (I + *IJ + *IJ0);
  rates[32] = *rho * *EJ0;
  rates[33] = *rho * *EJ;
  rates[34] = *epsilon * *EJ0;
  rates[35] = *epsilon * *EJ;
  rates[36] = *gamma * *IJ0;
  rates[37] = *gamma * *IJ;
  rates[38] = *p * *IJ0;
  rates[39] = *p * *IJ;
  rates[40] = *omega * *RJ0;
  rates[41] = *omega * *RJ;
  rates[42] = *m * *SM * Nall/ *k;
  rates[43] = *m * *EM * Nall/ *k;
  rates[44] = *m * *IM * Nall/ *k;
  rates[45] = *m * *RM * Nall/ *k;
  rates[46] = *beta * *SM * (I + *IJ + *IJ0);
  rates[47] = *rho * *EM;
  rates[48] = *epsilon * *EM;
  rates[49] = *gamma * *IM;
  rates[50] = *p * *IM;
  rates[51] = *omega * *RM;
  rates[52] = *mu * *SJ / 2.0;
  rates[53] = *mu * *EJ / 2.0;
  rates[54] = *mu * *IJ / 2.0;
  rates[55] = *mu * *RJ / 2.0;

  return rates;
}

//[[Rcpp::export]]
NumericVector simRcpp(IntegerVector state, NumericVector par, IntegerVector spec, int nevents, int tmax, int inc, double maxtau)
{
  nevents = 56;
  
  NumericVector out(state.size() * tmax * inc);
  out[0] = state[0];
  out[tmax * inc] = state[1];
  out[2 * tmax * inc] = state[2];
  out[3 * tmax * inc] = state[3];
  out[4 * tmax * inc] = state[4];
  out[5 * tmax * inc] = state[5];
  out[6 * tmax * inc] = state[6];
  out[7 * tmax * inc] = state[7];
  out[8 * tmax * inc] = state[8];
  out[9 * tmax * inc] = state[9];
  out[10 * tmax * inc] = state[10];
  out[11 * tmax * inc] = state[11];
  out[12 * tmax * inc] = state[12];
  out[13 * tmax * inc] = state[13];
  out[14 * tmax * inc] = state[14];
  out[15 * tmax * inc] = state[15];
  out[16 * tmax * inc] = state[16];
  
  int minstate = 10;
  int nt = 1;
  double t = 0.0;
  double mintau = 1.0/365/1000;
  bool hasEI = true;
  bool allpos;
  bool thingsHappen = true;
  double tau;
  IntegerVector varnew(state.size());
  IntegerVector varold(state.size());
  varold = clone(state);
  NumericVector rates(nevents);
  
  //NumericVector rates = getRates(var,par,nevents,t);
  //Rcout << rates;
  while (nt < tmax * inc && hasEI && thingsHappen){
    NumericVector newrates = getRates(state,par,nevents,t);
    
    for (int i = 0; i < rates.size(); i++){
      if (newrates[i] < 0) {
        Rcout << " par = ";
        Rcout << par;
        Rcout << " penultimate state = ";
        Rcout << varold;
        Rcout << " penultimate rates = ";
        Rcout << rates;
        Rcout << "; end state = ";
        Rcout << state;
        Rcout << " end rates = ";
        Rcout << newrates;
        Rcout << " allpos = ";
        Rcout << allpos;
        stop("negative rates");
      }
    }
    
    rates = newrates;
    if (sum(rates) == 0){
      thingsHappen = false;
      Rcout << "nothing happens";
      break;
    }
    //Rcout << "rates = ";
    //Rcout << rates;
    
    //double sumR = sum(rates);
    //double tau = R::rexp(1/sumR);
    tau = maxtau;
    //double mintau = pow(10,-4);
    //bool smallEnough = FALSE;
    //while (!smallEnough){
    //    if (tau <= maxtau){
    //        smallEnough = TRUE;
    //    }
    //    rates = getRates(var,par,nevents,t + maxtau);
    //    t = t + maxtau;
    //    tau = R::rexp(1/sumR);
    //}
    
    //if(tau > maxtau){
    
    if (state[0] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[1],rates[5]));
    }
    if (state[1] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[2],rates[6],rates[7]));
    }
    if (state[2] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[3],rates[8],rates[9]));
    }
    if (state[3] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[4],rates[10]));
    }
    if (state[4] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[12],rates[21],rates[30]));
    }
    if (state[5] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[13],rates[26],rates[31],rates[52]));
    }
    if (state[6] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[14],rates[22]));
    }
    if (state[7] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[16],rates[27],rates[33],rates[35],rates[53]));
    }
    if (state[8] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[18],rates[28],rates[37],rates[39],rates[54]));
    }
    if (state[9] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[20],rates[29],rates[41],rates[55]));
    }
    if (state[10] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[15],rates[23],rates[32],rates[34]));
    }
    if (state[11] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[17],rates[24],rates[36],rates[38]));
    }
    if (state[12] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[19],rates[25],rates[40]));
    }
    if (state[13] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[42],rates[46]));
    }
    if (state[14] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[43],rates[47],rates[48]));
    }
    if (state[15] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[44],rates[49],rates[50]));
    }
    if (state[16] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[45],rates[51]));
    }
    
    
    
    bool stepped = false;
    while (!stepped && hasEI){
      varnew = tauleap(state,rates,spec,nevents,tau);
      //Rcout << " varnew = ";
      //Rcout << varnew;
      //Rcout << "; t=";
      //Rcout << t*inc;
      //Rcout << ", tau=";
      //Rcout << tau*inc;
      
      allpos = true;
      for (int i=0; i < 17; i++){
        if (varnew[i] < 0.0) {
          allpos = false;
          //Rcout << " varold = ";
          //Rcout << varold;
          //Rcout << " ! negative value in new state ! ";
        }
      }
      
      if (allpos) {
        for (int i=0; i < 17; i++){
          varold[i] = state[i];
          state[i] = varnew[i];
          //Rcout << " tau leaped ";
        }
        stepped = true;
      } else {
        double sumR = sum(rates);
        double taunew = R::rexp(1.0/sumR);
        if (taunew > maxtau | tau/2.0 > mintau) {
          tau = tau/2.0;
          //Rcout << " tau halved ";
        } else {
          varold = state;
          tau = taunew;
          //  Rcout << "; old state = ";
          //  Rcout << state;
          state = stepModel(state,rates,spec,nevents);
          //  Rcout << "; old state = ";
          //  Rcout << varold;
          //  Rcout << "; rates = ";
          //  Rcout << rates;
          //  Rcout << ", tau = ";
          //  Rcout << tau;//print diagnostic info
          //  Rcout << ", t = ";
          //  Rcout << t;
          stepped = true;
        }
        //Rcout << "some negative, ";
      }
      if (state[1]==0 && state[2]==0 && state[7]==0 && state[10]==0 && state[8]==0 && state[11]==0 
            && state[14]==0 && state[15]==0) {
        hasEI = false;
      }
    }
    //} else {
    //  var = stepModel(var,rates,spec,nevents); 
    //}
    
    t += tau;
    
    if(t * inc > nt){
      out[nt] = state[0]; 
      out[tmax * inc + nt] = state[1];
      out[2 * tmax * inc + nt] = state[2];
      out[3 * tmax * inc + nt] = state[3];
      out[4 * tmax * inc + nt] = state[4];
      out[5 * tmax * inc + nt] = state[5];
      out[6 * tmax * inc + nt] = state[6];
      out[7 * tmax * inc + nt] = state[7];
      out[8 * tmax * inc + nt] = state[8];
      out[9 * tmax * inc + nt] = state[9];
      out[10 * tmax * inc + nt] = state[10];
      out[11 * tmax * inc + nt] = state[11];
      out[12 * tmax * inc + nt] = state[12];
      out[13 * tmax * inc + nt] = state[13];
      out[14 * tmax * inc + nt] = state[14];
      out[15 * tmax * inc + nt] = state[15];
      out[16 * tmax * inc + nt] = state[16];
      ++nt;
      //Rcout << nt;
    }
    
  }
  return out;
}
