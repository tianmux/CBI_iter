#include <iostream>
#include <complex>
//#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <string.h>
#include <sstream>
#include <iterator>
#include <map>
#include <cstdlib>
//#include <aligned_new>

//#include <advisor-annotate.h>
//#include <pstl/algorithm>
//#include <pstl/numeric>
#include <parallel/algorithm>
#include <parallel/numeric>
//#include <gsl/gsl_histogram.h>

#include"inputPara.h"

const double c_light = 299792458;
const double c_inver = 1/c_light;
const double pi = M_PI;
const double me = 9.1093837015e-31;
const double mp = 1.67262192369e-27;
const double mAu = 196.9665687*1.660540e-27;
const double qe = 1.6021766208e-19;
const double qAu = 79*qe;
const double E0p = 938.2720813e6;
const double E0e = 0.510998950e6;
const double E0Au = 196.9665687*931.5e6;

double E0 = E0p;

double qovermp = qe/mp;
double qoverme = qe/me;
double qovermAu = qAu/mAu;
int qovermp_on = 1;
int OMP_NUM_THREADS = 1;//omp_get_max_threads()/1;
double max_step = 1e0;

using namespace std::complex_literals;
std::complex<double> sig1 = 1;
std::complex<double> sig2 = 1i;
std::complex<double> Z0(double R,double Q, double omega0,std::complex<double> omegas){
    return R*omegas/(omegas+1i*Q*(omega0-omegas*omegas/omega0));
}
std::complex<double> Z(double R,double Q, double omega0,std::complex<double> omegas){
    double ReOmegas = omegas.real();
    double ImOmegas = omegas.imag();
    std::complex<double> result;
    double abs = sqrt((ReOmegas*omega0+2*ReOmegas*ImOmegas*Q)*(ReOmegas*omega0+2*ReOmegas*ImOmegas*Q)+Q*Q*(ImOmegas*omega0+omega0*omega0-ReOmegas*ReOmegas+ImOmegas*ImOmegas));

    return R*omegas/(omegas+1i*Q*(omega0-omegas*omegas/omega0));
}
std::complex<double> Err(int ps,int nRF, int i,int nBunch,int mu,int m,
                        std::vector<std::vector<std::complex<double>>> &p_omegas,
                        std::vector<std::vector<std::vector<std::complex<double>>>> &Zs,
                        std::vector<std::complex<double>> &temp,
                        std::vector<double> &R, std::vector<double> &QL,std::vector<double> &omegac,
                        std::complex<double> OMEGA_Init,
                        double omega0,double omegas, double factor,double sig_phi){
    // for each point of OMEGA, calculate the p_omegas array
    // then calculate the Zs at that frequency
    // then sum them up and times the 'factor' to it.
    //std::cout<<"factor: "<<factor<<std::endl;
    temp[i]=0;    
    // dynamic part
    for(int iRF= 0; iRF<nRF; ++iRF){
#pragma omp parallel for simd schedule(simd: static)
#pragma vector always
#pragma ivdep
        for(int j = 0; j<ps*2; ++j){
            p_omegas[i][j] = ((-ps+j)*nBunch+mu)*omega0+OMEGA_Init;
            Zs[i][j][iRF] = Z(R[iRF],QL[iRF],omegac[iRF],p_omegas[i][j]);
        }
    }
    
    for(int iRF= 0; iRF<nRF; ++iRF){
        for(int j = 0; j<ps*2; ++j){
            temp[i]+=pow(p_omegas[i][j]/omega0,2*m-1)*Zs[i][j][iRF]*exp(-double((-ps+j)*nBunch+mu)*double((-ps+j)*nBunch+mu)*sig_phi*sig_phi);
        }
    }
    temp[i] *= factor*1i;
    //std::cout<<"RHS :"<<temp[i]<<std::endl;
    return (OMEGA_Init*OMEGA_Init-omegas*omegas*m*m-temp[i])/omegas;
}

std::complex<double> ApproxOMEGA(std::vector<std::complex<double>> &OMEGA_Init0,std::vector<std::complex<double>> &OMEGA_Init, double mu, double nRF, double factor,
                                    std::vector<double> &R, std::vector<double> &QL,std::vector<double> &omegac,
                                    double nPar,double r0,double eta, double Gamma0,double T0,double Qs, double f0, double omega0, double omegas,
                                    double epsilon, double sig_phi,int ps, int nBunch,int maxIter,int m ){
    std::vector<std::complex<double>> temp(OMEGA_Init.size(),std::complex<double>(0,0));
    std::vector<std::vector<std::vector<std::complex<double>>>> Zs(OMEGA_Init.size(),std::vector<std::vector<std::complex<double>>>(ps*2,std::vector<std::complex<double>>(nRF,std::complex<double>(0,0))));
    std::vector<std::vector<std::complex<double>>> p_omegas(OMEGA_Init.size(),std::vector<std::complex<double>>(ps*2,std::complex<double>(0,0))); // the frequency points where we sample the impedance 
    
    return -Err(ps,nRF,0,nBunch,mu,m,p_omegas,Zs,temp,R,QL,omegac,omegas,omega0,omegas,factor,sig_phi)/2.0/omegas;
}
std::complex<double> NSolveOMEGA(std::vector<std::complex<double>> &OMEGA_Init0,std::vector<std::complex<double>> &OMEGA_Init, 
                                    double mu, double nRF, double factor,
                                    std::vector<double> &R, std::vector<double> &QL,std::vector<double> &omegac,
                                    double nPar,double r0,double eta, double Gamma0,double T0,double Qs, double f0, double omega0, double omegas,
                                    double epsilon, double sig_phi, int ps, int nBunch,int maxIter,int m){
    std::vector<std::complex<double>> err(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> err1(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> err2(OMEGA_Init.size(),std::complex<double>(1e10,0));

    std::vector<std::complex<double>> errRe(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> errIm(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr11(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr12(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr21(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr22(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> det(OMEGA_Init.size(),std::complex<double>(1,0));
    std::vector<std::complex<double>> dErr11inv(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr12inv(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr21inv(OMEGA_Init.size(),std::complex<double>(1e10,0));
    std::vector<std::complex<double>> dErr22inv(OMEGA_Init.size(),std::complex<double>(1e10,0));

    std::vector<std::complex<double>> OMEGA_1(OMEGA_Init.size(),std::complex<double>(0,0));
    std::vector<std::complex<double>> OMEGA_2(OMEGA_Init.size(),std::complex<double>(0,0));// the solution 1 and 2 steps ago.
    
    // stepsize for convergency 
    //std::complex<double> sig1 = 0.1;
    //std::complex<double> sig2 = 0.1i;
    std::vector<std::complex<double>> temp(OMEGA_Init.size(),std::complex<double>(0,0));
    std::vector<std::vector<std::vector<std::complex<double>>>> Zs(OMEGA_Init.size(),std::vector<std::vector<std::complex<double>>>(ps*2,std::vector<std::complex<double>>(nRF,std::complex<double>(0,0))));
    std::vector<std::vector<std::complex<double>>> p_omegas(OMEGA_Init.size(),std::vector<std::complex<double>>(ps*2,std::complex<double>(0,0))); // the frequency points where we sample the impedance 
    std::complex<double> OMEGA_result(0,0);
    
    double maxImOMEGA = -1e300;
    std::vector<std::complex<double>> OMEGA_step(OMEGA_Init.size(),std::complex<double>(0,0));
    std::vector<double> OMEGA_step_abs(OMEGA_Init.size(),0);
    // find the OMEGA for each starting point, see if there is any instability (Im(OMEGA)>0).
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
    for(int i = 0 ; i<OMEGA_Init.size(); ++i){
        int nIter = 0;
        OMEGA_Init[i] = OMEGA_Init0[i];
        OMEGA_1[i] = OMEGA_Init[i];
        OMEGA_2[i] = OMEGA_Init[i];
        while (abs(errRe[i]+1i*errIm[i])>epsilon & nIter<maxIter){
            OMEGA_1[i] = OMEGA_Init[i];
            err[i] = Err(ps,nRF,i,nBunch,mu,m,p_omegas,Zs,temp,R,QL,omegac,OMEGA_Init[i],omega0,omegas,factor,sig_phi);
            errRe[i] = real(err[i]);
            errIm[i] = imag(err[i]);
            err1[i] = Err(ps,nRF,i,nBunch,mu,m,p_omegas,Zs,temp,R,QL,omegac,OMEGA_Init[i]+sig1,omega0,omegas,factor,sig_phi);
            dErr11[i] = (real(err1[i])-errRe[i])/abs(sig1);
            dErr21[i] = (imag(err1[i])-errIm[i])/abs(sig2);
            err2[i] = Err(ps,nRF,i,nBunch,mu,m,p_omegas,Zs,temp,R,QL,omegac,OMEGA_Init[i]+sig2,omega0,omegas,factor,sig_phi);
            dErr12[i] = (real(err2[i])-errRe[i])/abs(sig1);
            dErr22[i] = (imag(err2[i])-errIm[i])/abs(sig2);
            det[i] = dErr11[i]*dErr22[i]-dErr21[i]*dErr12[i];
            dErr11inv[i] = dErr22[i]/det[i];
            dErr12inv[i] = -dErr12[i]/det[i];
            dErr21inv[i] = -dErr21[i]/det[i];
            dErr22inv[i] = dErr11[i]/det[i];
            OMEGA_step[i] = -errRe[i]*dErr11inv[i]-errIm[i]*dErr12inv[i]-1i*errRe[i]*dErr21inv[i]-1i*errIm[i]*dErr22inv[i];

            OMEGA_step_abs[i] = abs(OMEGA_step[i])/max_step;
            //std::cout<<abs(OMEGA_step[i])<<std::endl;
            //std::cout<<abs(OMEGA_step_abs[i])<<std::endl;

            OMEGA_step[i] = OMEGA_step_abs[i]>1?OMEGA_step[i]/OMEGA_step_abs[i]:OMEGA_step[i];
            //std::cout<<"Step: "<<(OMEGA_step[i])<<std::endl;
              
            //std::cout<<"Err: "<<(errRe[i]+1i*errIm[i])<<std::endl;
            OMEGA_Init[i] = OMEGA_Init[i]+OMEGA_step[i];//real(OMEGA_Init[i])-errRe[i]*dErr11inv[i]-errIm[i]*dErr12inv[i]+1i*(imag(OMEGA_Init[i])-errRe[i]*dErr21inv[i]-errIm[i]*dErr22inv[i]);
            //std::cout<<"OMEGA: "<<(OMEGA_Init[i])<<std::endl;        
            //std::cout<<"OMEGA1: "<<(OMEGA_1[i])<<std::endl;        
            //std::cout<<"OMEGA2: "<<(OMEGA_2[i])<<std::endl;        
            //std::cout<<"OMEGA-OMEGA2: "<<(OMEGA_Init[i]-OMEGA_2[i])<<std::endl;
            nIter++;
            if(abs(OMEGA_Init[i]-OMEGA_2[i])<epsilon){
                OMEGA_Init[i] = (OMEGA_Init[i]+OMEGA_1[i])/2.0;
            }
            OMEGA_2[i] = OMEGA_1[i];
            temp[i] = 0;
        }
        
#if 0
        if(nIter == maxIter){
            std::cout<<"Final Error = "<<err[i]<<std::endl;
        }
        std::cout<<"# of iteration: "<<nIter<<std::endl;
        std::cout<<"Final Error = "<<err[i]<<std::endl;
        std::cout<<"Final OMEGA = "<<OMEGA_Init[i]<<std::endl;
#endif
    }
    //std::cout<<"Got one. "<<std::endl;
    // at least return the first result. 
    OMEGA_result = OMEGA_Init[0];
    // find the OMEGA that has the largest imaginary part.
    int idx = 0;
    for(int i = 0;i<OMEGA_Init.size();++i){
        if(abs(err[i])<epsilon){
            OMEGA_result = maxImOMEGA<OMEGA_Init[i].imag()?OMEGA_Init[i]:OMEGA_result;
            maxImOMEGA = maxImOMEGA<OMEGA_Init[i].imag()?OMEGA_Init[i].imag():maxImOMEGA;
            idx = i;
        }
    }
    //std::cout<<"Final Error = "<<err[idx]<<std::endl;
    //std::cout<<"Final OMEGA = "<<OMEGA_Init[idx]<<std::endl;
    return OMEGA_result;
}
double Fact(int N){
    if(N==0){
        return 1;
    }
    else{
        int temp=1;
        for(int i = 1;i<=N;++i){
            temp*=i;
        }
        return temp;
    }
}
int main(){
    std::cout.precision(17);
    omp_set_num_threads(OMP_NUM_THREADS);
    int nRF;
    std::vector<double> h;
    std::vector<double> RoQ;
    std::vector<double> RoQ_Acc;
    std::vector<double> NC;
    std::vector<double> Vs;
    std::vector<double> Vq;
    std::vector<double> Vc;
    std::vector<double> pBeam;
    std::vector<double> Phis;
    std::vector<double> dfMin;
    std::vector<double> dfMax;
    std::vector<double> n_df;// number of frequency samples for each RF
    double V0 = 0;//23.7e6/(NC[0]+NC[1]);
    
    int nBunch = 1260;
    
    double NperBunch = 17.2e10;    // number of real particles.
    double ringR = 610.1754;
    double GMTSQ = 961;
    double pRad0 = 10e6;
    double Gamma0 = 19600;
    double detuneMin = 0;
    double detuneMax = 0;//opt_Detune+1e4;
    double sig_z=2.3e-3; // rms bunch length in unit of seconds

    int maxIter = 30; // total number of iteration for convergency.
    int nIb = 10; // number samples for Ib
    int nDetune = 1; // number of samples for Detune
    int nMu = 2;// number of modes to invesitgate. 
    int m = 1;// the type of the coupled bunch mode, 1 is dipole, 2 is quad...
    int nOmegaGuess = 11;
    int ps = nBunch; //2*nBunch+1; // number of frequency samples, actually is 2*ps+1
    
    double epsilon = 1e-9; // error for convergency estimation
    
    //std::vector<double> h;
    // read the input
    std::ifstream inputfile;
    inputfile.open ("input.txt",std::ios::in);
    if (!inputfile) {
        std::cout<< "Unable to open file datafile.txt";
        exit(1);   // call system to stop
    }
    std::string x;
    std::string delimiter = "\t";
    // first must be 'nRF'
    while (getline (inputfile, x)) {
        std::istringstream ss(x);
        std::string substr;
        
        while (getline(ss, substr, '\t')){
            if (substr == "nRF"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                nRF = stod(substr);
                h.resize(nRF);
                RoQ.resize(nRF);
                RoQ_Acc.resize(nRF);
                NC.resize(nRF);
                dfMin.resize(nRF);
                dfMax.resize(nRF);
                n_df.resize(nRF);
                Vs.resize(nRF);
                Vq.resize(nRF);
                Vc.resize(nRF);
                pBeam.resize(nRF);
                Phis.resize(nRF);
                std::cout<<nRF<<std::endl;
                break;
            }
            if (substr == "h"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    h[i] = stod(substr);
                    std::cout<<h[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "RoQ"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    RoQ[i] = stod(substr);
                    RoQ_Acc[i] = RoQ[i]*2;
                    std::cout<<RoQ[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "NC"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    NC[i] = stod(substr);
                    std::cout<<NC[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "n_df"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    n_df[i] = stod(substr);
                    std::cout<<n_df[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "dfMin"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    dfMin[i] = stod(substr);
                    std::cout<<dfMin[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "dfMax"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    dfMax[i] = stod(substr);
                    std::cout<<dfMax[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "Vs"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    Vs[i] = stod(substr);
                    std::cout<<Vs[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "Vq"){
                std::cout<<substr<< ":\t";
                int i = 0;
                while(getline(ss, substr, '\t')){
                    Vq[i] = stod(substr);
                    std::cout<<Vq[i]<<'\t';
                    i++;
                }
                std::cout<<std::endl;
                break;
            }
            if (substr == "V00"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                V0 = stod(substr);
                std::cout<<V0<<std::endl;
                break;
            }
            
            if (substr == "nBunch"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                nBunch = stod(substr);
                std::cout<<nBunch<<std::endl;
                break;
            }
            if (substr == "nPB0"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                NperBunch = stod(substr);
                std::cout<<NperBunch<<std::endl;
                break;
            }
            if (substr == "R_ring"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                ringR = stod(substr);
                std::cout<<ringR<<std::endl;
                break;
            }
            if (substr == "GMTSQ"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                GMTSQ = stod(substr);
                std::cout<<GMTSQ<<std::endl;
                break;
            }
            if (substr == "pRad0"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                pRad0 = stod(substr);
                std::cout<<pRad0<<std::endl;
                break;
            }
            if (substr == "Gamma0"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                Gamma0 = stod(substr);
                std::cout<<Gamma0<<std::endl;
                break;
            }
            if (substr == "maxIter"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                maxIter = stod(substr);
                std::cout<<maxIter<<std::endl;
                break;
            }
            if (substr == "nIb"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                nIb = stod(substr);
                std::cout<<nIb<<std::endl;
                break;
            }
            
            if (substr == "nMu"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                nMu = stod(substr);
                std::cout<<nMu<<std::endl;
                break;
            }
            if (substr == "m"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                m = stod(substr);
                std::cout<<m<<std::endl;
                break;
            }
            if (substr == "nOmegaGuess"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                nOmegaGuess = stod(substr);
                std::cout<<nOmegaGuess<<std::endl;
                break;
            }
            if (substr == "ps"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                ps = stod(substr);
                std::cout<<ps<<std::endl;
                break;
            }
            if (substr == "epsilon"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                epsilon = stod(substr);
                std::cout<<epsilon<<std::endl;
                break;
            }
            if (substr == "max_step"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                max_step = stod(substr);
                std::cout<<max_step<<std::endl;
                break;
            }
            if (substr == "sig1"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                sig1 = stod(substr);
                std::cout<<sig1<<std::endl;
                break;
            }
            if (substr == "sig2"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                sig2 = stod(substr)*1i;
                std::cout<<sig2<<std::endl;
                break;
            }
            if (substr == "n_thread"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                OMP_NUM_THREADS = int(stod(substr));
                std::cout<<OMP_NUM_THREADS<<std::endl;
                break;
            }
            if (substr == "sig_z"){
                std::cout<<substr<< ":\t";
                getline(ss, substr, '\t');
                sig_z = stod(substr);
                std::cout<<sig_z<<std::endl;
                break;
            }
        }
    }

    
    pRad0 = pRad0/(NC[0]+NC[1]); // normalize to one cavity
    V0 = V0/(std::accumulate(NC.begin(),NC.end(),0));      // normalize to one cavity
    for(int iRF = 0; iRF<nRF; ++iRF){
        Vc[iRF] = sqrt(Vs[iRF]*Vs[iRF]+Vq[iRF]*Vq[iRF]);
        Phis[iRF] = acos(Vs[iRF]/Vc[iRF]);
    }
    double Ek = Gamma0*E0e;
    double f0 = (c_light*sqrt(1-1/Gamma0/Gamma0))/(2*M_PI*ringR);
    double omega0 = f0*2*pi;
    double T0 = 1/f0;
    double omegarf = f0*h[0]*2*pi;
    double sig_t = sig_z/c_light;
    double sig_phi = omega0*sig_t;
    double r0 = (qe*qe)/(me*c_light*c_light);
    double eta = 1/GMTSQ-1/(Gamma0*Gamma0);
    ps = int(5*sqrt(2)/sig_t/nBunch/omega0);
    double nPerBunch0 = NperBunch;

    double IbDC0 = nPerBunch0*f0*qe*nBunch;
    // Voltages per cavity needed
    
    double Vqtot = 0;

    double Qs = 0; 
    for(int i = 0;i<nRF;++i){
        Vqtot+=Vq[i]*NC[i];
    }
    std::cout<<"eta : "<<eta<<std::endl;
    std::cout<<"h : "<<h[0]<<std::endl;
    std::cout<<"Ek : "<<Ek<<std::endl;
    std::cout<<"sig_t : "<<sig_t<<std::endl;
    std::cout<<"omega0 : "<<omega0<<std::endl;

    std::cout<<"ps : "<<ps<<std::endl;

    Qs = sqrt(h[0]*eta*Vqtot/(2*pi*Ek));
    
    double omegas = Qs*f0*2*pi;
    
    
    std::cout<<"nRF : "<<nRF<<std::endl;
    std::cout<<"f0 : "<<f0<<std::endl;
    std::cout<<"Qs : "<<Qs<<std::endl;

    std::vector<double> nPerBunch(nIb,0);
    std::vector<double> IbDC(nIb,0);
    std::vector<double> factor(nIb,0); // the factor infront of the Sum(((pM+mu)omega0+omegas))Z((pM+mu)omega0+omegas))
    std::vector<double> pRad(nIb,0); // radiation power, use this to calculate optimum QL
    std::vector<double> mus(nMu,0);
    std::vector<std::vector<double>> R(nIb,std::vector<double>(nRF,0)); // Impedance of each cavity, for different cases of Ib
    std::vector<std::vector<double>> Q(nIb,std::vector<double>(nRF,0)); // Optimum QL for each cavity for different cases of Ib
    for(int i = 0;i<n_df.size();++i){
        nDetune *= n_df[i];
    }
    std::cout<<"Tot # detune : "<<nDetune<<std::endl;
    std::vector<std::vector<double>> omegac(nDetune,std::vector<double>(nRF,0)); // set of frequencies of different cavities. 
    std::vector<std::vector<std::vector<std::complex<double>>>> OMEGAS(nMu,std::vector<std::vector<std::complex<double>>>(nIb,std::vector<std::complex<double>>(nDetune,std::complex<double>(0,0)))); //store the solutions for OMEGA
    std::vector<std::vector<std::vector<std::complex<double>>>> approxOMEGAS(nMu,std::vector<std::vector<std::complex<double>>>(nIb,std::vector<std::complex<double>>(nDetune,std::complex<double>(0,0)))); //store the solutions for OMEGA
    
    std::vector<std::complex<double>> OMEGA_init0(nOmegaGuess,std::complex<double>(0,0)); // initial guess for OMEGA, to keep the record
    std::vector<std::complex<double>> OMEGA_init(nOmegaGuess,std::complex<double>(0,0)); // initial guess for OMEGA, for calculation


    // initialize the mus
    for (int i = 0; i<nMu;++i){
        mus[i] =(nBunch-i);
    }
    // initialize the parameters related to the beam current
    for (int i = 0; i < nIb; ++i){
        nPerBunch[i] = nPerBunch0/nIb*(i+1);
        IbDC[i] = nPerBunch[i]*qe*f0*nBunch;
        pRad[i] = pRad0/IbDC0*IbDC[i];
        factor[i] = nBunch*nPerBunch[i]*r0*eta/Gamma0/T0/T0/T0*4*pi*m/Fact(m-1)/pow(2,m)*pow(sig_phi,2*m-2);
        std::cout<<"Factor = "<<factor[i]<<std::endl;
        for (int j = 0; j<nRF; ++j){
            Q[i][j] = Vc[j]*Vc[j]/(RoQ_Acc[j]*pRad[i]);
            R[i][j] = RoQ[j]*Q[i][j]*NC[j];
            std::cout<<"Q an R : "<<Q[i][j]<<','<<R[i][j]<<std::endl;
        }
    }

    //Change the detune of each RF:
    if(nRF > 2){
        std::cout<<"Number of RF is not supported, coming soon. "<<std::endl;
        exit(1);
    }
    if(nRF==2){
        for(int i = 0;i<n_df[0];++i){
            for(int j = 0;j<n_df[1];++j){
                omegac[i*n_df[1]+j][0] = (dfMin[0]+(dfMax[0]-dfMin[0])/(n_df[0]-1)*i)*2*pi+omegarf;
                omegac[i*n_df[1]+j][1] = (dfMin[1]+(dfMax[1]-dfMin[1])/(n_df[1]-1)*j)*2*pi+omegarf;
                std::cout<<"dfC: "<<(omegac[i*n_df[1]+j][0]-omegarf)/2/pi<<','<<(omegac[i*n_df[1]+j][1]-omegarf)/2/pi<<std::endl;
            }
        }
    }
    else if(nRF==1){
        std::cout<<"Change the detune of  the focusing cavity."<<std::endl;
        for (int i = 0; i < nDetune; ++i){
            omegac[i][0] = (dfMin[0]+(dfMax[0]-dfMin[0])/n_df[0]*i)*2*pi+omegarf;
            std::cout<<"dfc : "<<(omegac[i][0]-omegarf)/2/pi<<std::endl;
        }
    }

    // initialize the initial guess for OMEGA
    for (int i = 0;i < nOmegaGuess; ++i){
        OMEGA_init0[i] = omegas*exp(1i/std::complex<double>(nOmegaGuess)*std::complex<double>((i-int(nOmegaGuess/2))*1.0*pi));//-1i*1.5*omegas;
        OMEGA_init[i] = OMEGA_init0[i];
        std::cout<<OMEGA_init[i]<<std::endl;
    }
    auto t_start = omp_get_wtime(); 
    
    for ( int k = 0; k<nMu; ++k){
        std::cout<<"Calculating mu = "<<k<<std::endl;
        for ( int i = 0; i < nIb; ++i){
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
            for(int  j = 0; j < nDetune; ++j){
                std::vector<std::complex<double>> OMEGA_init(nOmegaGuess,std::complex<double>(0,0)); // initial guess for OMEGA, for calculation
                //std::cout<<"Ib : "<<IbDC[i]<<','<<"Detune : "<<(omegac[j][0]-omegarf)/2/pi<<std::endl;
                OMEGAS[k][i][j] = NSolveOMEGA(OMEGA_init0,OMEGA_init,mus[k],nRF,factor[i],R[i],Q[i],omegac[j],nPerBunch[i],r0,eta,Gamma0,T0,Qs,f0,omega0,omegas,epsilon,sig_phi,ps,nBunch, maxIter,m);
            }
            std::cout<<i+1<<'/'<<nIb<<std::endl;
        }
    }
    auto t_end = omp_get_wtime(); 
    std::cout<<"Iterative solver takes : "<<(t_end-t_start) <<" [s]. "<<std::endl;

#if 0
    t_start = omp_get_wtime();
    for ( int k = 0; k<nMu; ++k){
        for ( int i = 0; i < nIb; ++i){
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
            for(int  j = 0; j < nDetune; ++j){
                //std::cout<<"Ib : "<<IbDC[i]<<','<<"Detune : "<<(omegac[j][0]-omegarf)/2/pi<<std::endl;
                approxOMEGAS[k][i][j] = ApproxOMEGA(OMEGA_init0,OMEGA_init,mus[k],nRF,factor[i],R[i],Q[i],omegac[j],nPerBunch[i],r0,eta,Gamma0,T0,Qs,f0,omega0,omegas,epsilon,ps,nBunch, maxIter,m);
            }
            std::cout<<i+1<<'/'<<nIb<<std::endl;
        }
    }
    t_end = omp_get_wtime(); 
    std::cout<<"Approximation solver takes : "<<(t_end-t_start) <<" [s]. "<<std::endl;
#endif


    // output
    // first write the header which is the detune.
    if(nRF==1){
        std::ofstream ImOmegafile;
        for ( int k = 0; k<nMu; ++k){
                ImOmegafile.open ("ImOmega"+std::to_string(int(mus[k]))+".txt",std::ios::out);

                for(int iRF = 0;iRF<nRF;++iRF){
                    //ImOmegafile<<'x'<<',';
                    for( int j = 0; j < nDetune; ++j){
                        ImOmegafile<<(omegac[j][iRF]-omegarf)/2/pi<<',';
                    }
                    ImOmegafile<<std::endl;
                }
                for ( int i = 0; i < nIb; ++i){
                    ImOmegafile<<IbDC[i]<<',';
                    for( int j = 0; j < nDetune; ++j){
                        ImOmegafile<<OMEGAS[k][i][j].imag()<<',';
                    }
                    ImOmegafile<<std::endl;
                }
                ImOmegafile.close();
            }
        std::ofstream ReOmegafile;
        // first write the header which is the detune.
        for ( int k = 0; k<nMu; ++k){
            ReOmegafile.open ("ReOmega"+std::to_string(int(mus[k]))+".txt",std::ios::out);
            for(int iRF = 0;iRF<nRF;++iRF){
                for( int j = 0; j < nDetune; ++j){
                    ReOmegafile<<(omegac[j][iRF]-omegarf)/2/pi<<',';
                }
                ReOmegafile<<std::endl;
            }
            
            for ( int i = 0; i < nIb; ++i){
                ReOmegafile<<IbDC[i]<<',';
                for( int j = 0; j < nDetune; ++j){
                    ReOmegafile<<OMEGAS[k][i][j].real()<<',';
                }
                ReOmegafile<<std::endl;
            }
            ReOmegafile.close();
    }
    }
    

    // write the data in the table format with colums equal to df of first RF and index is df of second RF
    if(nRF==2 && nIb==1){
        std::ofstream ImOmegafile;
        for ( int k = 0; k<nMu; ++k){
            ImOmegafile.open ("ImOmega"+std::to_string(int(mus[k]))+".txt",std::ios::out);
            ImOmegafile<<IbDC[0]<<',';
            for( int df0 = 0; df0 < n_df[0]; ++df0){
                ImOmegafile<<(omegac[df0*n_df[1]][0]-omegarf)/2/pi<<',';
            }
            ImOmegafile<<std::endl;
            for( int df1 = 0; df1 < n_df[1]; ++df1){
                ImOmegafile<<(omegac[df1][1]-omegarf)/2/pi<<',';
                for(int df0 = 0;df0<n_df[0];++df0){
                    ImOmegafile<<OMEGAS[k][0][df0*n_df[1]+df1].imag()<<',';
                }
                ImOmegafile<<std::endl;
            }
            ImOmegafile.close();
        }
        std::ofstream ReOmegafile;
        for ( int k = 0; k<nMu; ++k){
            ReOmegafile.open ("ReOmega"+std::to_string(int(mus[k]))+".txt",std::ios::out);
            ReOmegafile<<IbDC[0]<<',';
            for( int df0 = 0; df0 < n_df[0]; ++df0){
                ReOmegafile<<(omegac[df0*n_df[1]][0]-omegarf)/2/pi<<',';
            }
            ReOmegafile<<std::endl;
            for( int df1 = 0; df1 < n_df[1]; ++df1){
                ReOmegafile<<(omegac[df1][1]-omegarf)/2/pi<<',';
                for(int df0 = 0;df0<n_df[0];++df0){
                    ReOmegafile<<OMEGAS[k][0][df0*n_df[1]+df1].real()<<',';
                }
                ReOmegafile<<std::endl;
            }
            ReOmegafile.close();
        }
    }
    
    
#if 0
    // output
    std::ofstream ApproxImOmegafile;
    // first write the header which is the detune.
    for ( int k = 0; k<nMu; ++k){
        ApproxImOmegafile.open ("ApproxImOmega"+std::to_string(int(mus[k]))+".txt",std::ios::out);
        for( int j = 0; j < nDetune; ++j){
            ApproxImOmegafile<<(omegac[j][0]-omegarf)/2/pi<<',';
        }
        ApproxImOmegafile<<std::endl;
        for ( int i = 0; i < nIb; ++i){
            ApproxImOmegafile<<IbDC[i]<<',';
            for( int j = 0; j < nDetune; ++j){
                ApproxImOmegafile<<approxOMEGAS[k][i][j].imag()/2/pi<<',';
            }
            ApproxImOmegafile<<std::endl;
        }
        ApproxImOmegafile.close();
    }
    for ( int k = 0; k<nMu; ++k){
        ReOmegafile.open ("ApproxReOmega"+std::to_string(int(mus[k]))+".txt",std::ios::out);
        for( int j = 0; j < nDetune; ++j){
            ReOmegafile<<(omegac[j][0]-omegarf)/2/pi<<',';
        }
        ReOmegafile<<std::endl;
        for ( int i = 0; i < nIb; ++i){
            ReOmegafile<<IbDC[i]<<',';
            for( int j = 0; j < nDetune; ++j){
                ReOmegafile<<approxOMEGAS[k][i][j].real()/2/pi<<',';
            }
            ReOmegafile<<std::endl;
        }
        ReOmegafile.close();
    }
#endif
#if 0
    std::cout<<f0<<std::endl;
    std::cout<<Vs0<<std::endl;
    std::cout<<Vq0<<std::endl;
    std::cout<<"Qs : " <<Qs<<std::endl;
    std::cout<<sqrt(Vs0*Vs0+14.0/4*Vq0*14/4*Vq0)<<std::endl;
    std::cout<<Vnew0<<std::endl;
    std::cout<<Phis/pi*180<<std::endl;
    std::cout<<opt_Detune<<std::endl;
    std::cout<<NC[0]<<std::endl;
#endif
    return 0;
}
