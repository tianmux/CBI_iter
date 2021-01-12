#ifndef INPUTPARA_H
#define INPUTPARA_H


// Class of input parameters
class inputPara{
public:
    double c = 299792458;
    double E0 = 938.2720813e6; //rest energy of proton in unit of eV
    double q = 1.6021766208e-19; // charge of proton
    int cT = 0; // current turn

    // general parameters
    std::map<std::string,unsigned int> generalPara;
   
    // ring parameters
    std::map<std::string,double> ringPara;
    // RF parameters
    std::map<std::string,std::vector<double>> rfPara;
    std::map<std::string,std::vector<double>> apPara;
    // Modulation
    std::map<std::string,std::vector<double>> modPara;

    // bunch parameters

    std::map<std::string,double> bunchPara;
    
    inputPara(){
        c = 299792458;
        E0 = 938.2720813e6; //rest energy of proton
        q = 1.6021766208e-19; // charge of proton
        cT = 0; // current turn
    
        // general parameters
        generalPara={
            {"type",1}, // 1 for electron, 0 for proton
            {"dynamicOn",1},
            {"n_dynamicOn",9000},
            {"n_turns",1e8},
            {"turn_start",0},
            {"n_per_show",2},
            {"n_bunches",1},
            {"n_fill",10},// number of turn we start to fill the bunch into the bucket.
            {"n_q_ramp",2000},
            {"n_detune_start",40000}, // turn number when we start to ramp the frequency
            {"n_detune_ramp",80000}, // number of turns it takes to ramp the frequency.
            {"detune_slow_factor",10},
            {"n_detune_ramp_tot",100000}, // number of turns for second ramp.
            {"n_I_ramp_start",2000},// turn number when we start to ramp the Ig
            {"n_I_ramp_end",10000},// turn number when we end the ramp on Ig.
            {"mainRF",0}, // the main RF used for initial bunch distribution calculation.
            {"main_detune",0},
			{"Plot",0},
            {"init",0}, // whether this is a simulation from some previous result (0 means not from previous)
            {"write_break",1}, // whether we write the breaking point info to the disk.
            {"step_store",1},
            {"n_threads",18},
 			{"Slice",0}
        };
       
        // ring parameters
        ringPara={
			{"R",610.1754}, // radius of the ring
			{"GMTSQ",552.25}, // gamma_t square
			{"Gamma",293},
			{"Ek",0},//ringPara["Gamma"]*E0},
			{"f0",0},//c/(2*M_PI*ringPara["R"])},
			{"T0",0},// 1/ringPara["f0"]},
			{"omg0",0},//2*M_PI*ringPara["f0"]},
			{"eta",0},//1/ringPara["GMTSQ"]-1/(ringPara["Gamma"]*ringPara["Gamma"])}
			{"nRF",1},
            {"nRF1",1},
			{"nRF2",1},
			{"nRFc",1},
            {"Prad",1e7}, // synchrotron radiation power.
            {"t_rad_long",1},
            {"Ek_damp",3.0},
            {"nHOM",0}
        };
        // RF parameters
        rfPara={
            {"h",std::vector<double>(1)},
			{"RoQ",std::vector<double>(0)},
            {"Qs",std::vector<double>(0)},// This is the synchrotron tune, sqrt(h[0]*V[0]*eta*abs(cos(phis[0]))/(2*M_PI*Ek));
			{"QL",std::vector<double>(0)}, // this is the loaded Q of each mode. 
            {"Vref_I",std::vector<double>(0)},
			{"Vref_Q",std::vector<double>(0)} ,
            {"Iref_I",std::vector<double>(0)},
            {"Iref_Q", std::vector<double>(0)},
            {"I_I_ref_ini", std::vector<double>(0)},
            {"I_I_ref_final",std::vector<double>(0)},
            {"I_Q_ref_ini", std::vector<double>(0)},
            {"I_Q_ref_final",std::vector<double>(0)}
        };
        apPara={
			{"fa",std::vector<double>(1)},
			{"phia",std::vector<double>(1)},
            {"n_ramp",std::vector<double>(1)}, // number of turns it takes to ramp up the current.
			{"amplitude",std::vector<double>(1)},
            {"delay",std::vector<double>(1)}, // delay of the feedback.
			{"detune",std::vector<double>(1)},
            {"detune_ini",std::vector<double>(1)},
            {"detune_mid",std::vector<double>(1)}, // intermediate frequency where we switch tuning speed. 
            {"detune_final",std::vector<double>(1)},
            {"PA_cap",std::vector<double>(1)}, // the cap of the power amplifier output. in unit of Iref's
            {"nCav", std::vector<double>(1)},// number of cavities, for the radiation calculation. 
            {"n_fb_on",std::vector<double>(1)}
		};
        // Modulation
        modPara={
        {"Vepsilon",std::vector<double>(0)}, // modulation depth
        {"Mum",std::vector<double>(0)}, // modulation tune factor
        };
    
        // bunch parameters
    
        bunchPara={
            {"nBeam",1},
            {"beam_shift",0},
            {"Npar",1e6},
			{"NperBunch",2.7e9},
			{"N_bins",1e3},
            {"A",0.8},
            {"fill_step",1},
            {"siglong",1},
            {"epsln",bunchPara["A"]/6},//bunchPara["A"]/6},
            {"delta_hat",0},//sqrt(epsln)*sqrt(ringPara["omg0"]/(M_PI*ringPara["Ek"]))*std::pow((rfPara["h"][0]*rfPara["V"][0]*abs(cos(rfPara["phis"][0]))/(2*M_PI*ringPara["Ek"]*ringPara["eta"])),0.25)},
            {"t_hat",0}//bunchPara["delta_hat"]/rfPara["Qs"]*ringPara["eta"]/ringPara["omg0"]
            
        };
    };
	bool inRFPara(std::vector<std::string> V,std::vector<std::string> T, std::vector<std::string> TP,std::vector<std::string> Phi,std::string para){
		for (int i = 0;i<V.size();++i){
			if (V[i] == para|T[i] == para|TP[i] == para|Phi[i]==para){ 
				return true;
			}
		}
		return false;
	};
    int read(std::string path){
        // Object to write in file
        std::ifstream fin;
        fin.open(path, std::ios::in);
		std::vector<std::string> rfVLines;
		std::vector<std::string> rfTLines;
		std::vector<std::string> rfTPLines;
		std::vector<std::string> rfPhiLines;
        std::string temp;
        // Read the object's data in file
        while(std::getline(fin,temp)){
            if(temp.size()>1){
                std::stringstream ss(temp);
                std::istream_iterator<std::string> begin(ss);
                std::istream_iterator<std::string> end;
                std::vector<std::string> vstrings(begin, end);
                if(vstrings.size()>1){
                    if(vstrings[0]=="nRF"){
						int nRF = stoi(vstrings[1]);
						rfVLines.resize(nRF);
						rfTLines.resize(nRF);
						rfTPLines.resize(nRF);
						rfPhiLines.resize(nRF);

						for (int i = 0;i<nRF;++i){
							rfVLines[i] = "V"+std::to_string(i);
							rfTLines[i] = "TV"+std::to_string(i);
							rfTPLines[i] = "TP"+std::to_string(i);
							rfPhiLines[i] = "Phi"+std::to_string(i);
						}
                        rfPara["RoQ"].resize(nRF);
						rfPara["h"].resize(nRF);
						rfPara["QL"].resize(nRF);
                        rfPara["Vref_I"].resize(nRF);
                        rfPara["Vref_Q"].resize(nRF);
                        rfPara["Iref_I"].resize(nRF);
                        rfPara["Iref_Q"].resize(nRF);
                        rfPara["I_I_ref_ini"].resize(nRF);
                        rfPara["I_I_ref_final"].resize(nRF);
                        rfPara["I_Q_ref_ini"].resize(nRF);
                        rfPara["I_Q_ref_final"].resize(nRF);
                        apPara["fa"].resize(nRF);
                        apPara["phia"].resize(nRF);
                        apPara["n_ramp"].resize(nRF);
                        apPara["amplitude"].resize(nRF);
                        apPara["detune"].resize(nRF);
                        apPara["detune_ini"].resize(nRF);
                        apPara["detune_mid"].resize(nRF);
                        apPara["detune_final"].resize(nRF);
                        apPara["delay"].resize(nRF);
                        apPara["n_fb_on"].resize(nRF);
                        apPara["PA_cap"].resize(nRF);
                        apPara["nCav"].resize(nRF);
                        apPara["gII"].resize(nRF);
                        apPara["gQQ"].resize(nRF);
                        apPara["gIQ"].resize(nRF);
                        apPara["gQI"].resize(nRF);
                        apPara["gIIi"].resize(nRF);
                        apPara["gQQi"].resize(nRF);
                        apPara["gIQi"].resize(nRF);
                        apPara["gQIi"].resize(nRF);
                    }
                    std::cout<<vstrings[0]<<std::endl;
					if(inRFPara(rfVLines,rfTLines,rfTPLines,rfPhiLines,vstrings[0])){
						rfPara[vstrings[0]].resize(vstrings.size()-1);
						for(unsigned int i=0;i<vstrings.size()-1;++i){
                            rfPara[vstrings[0]][i]=stod(vstrings[i+1]); // store the ramping of the voltage.
                        }
					}

                    if(vstrings[0]=="I_I_ref_ini"|vstrings[0]=="I_I_ref_final"|vstrings[0]=="I_Q_ref_ini"|vstrings[0]=="I_Q_ref_final"|vstrings[0]=="Iref_I"|vstrings[0]=="Iref_Q"|vstrings[0]=="RoQ"|vstrings[0]=="h"|vstrings[0]=="QL"|vstrings[0]=="Vref_I"|vstrings[0]=="Vref_Q"){
                        for(unsigned int i=0;i<rfPara[vstrings[0]].size();++i){
                            rfPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }

                    if(vstrings[0]=="n_dynamicOn"|vstrings[0]=="dynamicOn"|vstrings[0]=="n_I_ramp_start"|vstrings[0]=="n_I_ramp_end"|vstrings[0]=="n_threads"|vstrings[0]=="main_detune"|vstrings[0]=="detune_slow_factor"|vstrings[0]=="n_detune_ramp_tot"|vstrings[0]=="mainRF"|vstrings[0]=="type"|vstrings[0]=="n_detune_start"|vstrings[0]=="n_detune_ramp"|vstrings[0]=="n_q_start"|vstrings[0]=="n_q_ramp"|vstrings[0]=="turn_start"|vstrings[0]=="init"|vstrings[0]=="step_store"|vstrings[0]=="Slice"|vstrings[0]=="Plot"|vstrings[0]=="n_turns"|vstrings[0]=="n_per_show"|vstrings[0]=="n_bunches"|vstrings[0]=="n_fill"|vstrings[0]=="turn_btw_fill"){
                        generalPara[vstrings[0]]=int(stod(vstrings[1]));
                    }

                    if(vstrings[0]=="Ek_damp"|vstrings[0]=="nRF1"|vstrings[0]=="nRFc"|vstrings[0]=="nRF2"|vstrings[0]=="t_rad_long"|vstrings[0]=="Prad"|vstrings[0]=="nHOM"|vstrings[0]=="R"|vstrings[0]=="GMTSQ"|vstrings[0]=="Gamma"|vstrings[0]=="nRF"){
                        ringPara[vstrings[0]]=stod(vstrings[1]);
                    }

                    if(vstrings[0]=="Vepsilon"|vstrings[0]=="Mum"||vstrings[0]=="Tesp"||vstrings[0]=="Tmu"){
                        modPara[vstrings[0]].resize(vstrings.size()-1);
						for(unsigned int i=0;i<vstrings.size()-1;++i){
							modPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }

                    if(vstrings[0]=="beam_shift"|vstrings[0]=="nBeam"|vstrings[0]=="NperBunch"|vstrings[0]=="siglong"|vstrings[0]=="fill_step"|vstrings[0]=="Charge"|vstrings[0]=="Npar"|vstrings[0]=="A"|vstrings[0]=="N_bins"|vstrings[0]=="N_bnches_p_train"|vstrings[0]=="xmin"|vstrings[0]=="xmax"|vstrings[0]=="ymin"|vstrings[0]=="ymax"){
                        bunchPara[vstrings[0]]=stod(vstrings[1]);
                    }

					if(vstrings[0]=="nCav"|vstrings[0]=="PA_cap"|vstrings[0]=="detune_mid"|vstrings[0]=="detune_ini"|vstrings[0]=="detune_final"|vstrings[0]=="n_fb_on"|vstrings[0]=="gIIi"|vstrings[0]=="gQQi"|vstrings[0]=="gIQi"|vstrings[0]=="gQIi"|vstrings[0]=="gII"|vstrings[0]=="gQQ"|vstrings[0]=="gIQ"|vstrings[0]=="gQI"|vstrings[0]=="delay"|vstrings[0]=="fa"|vstrings[0]=="phia"|vstrings[0]=="amplitude"|vstrings[0]=="detune"|vstrings[0]=="n_ramp"){
                        for(unsigned int i=0;i<apPara[vstrings[0]].size();++i){
                            apPara[vstrings[0]][i]=stod(vstrings[i+1]);
                        }
                    }

                }
            }
        }
        std::cout<<"Read parameters in."<<std::endl;
        ringPara["Ek"]=ringPara["Gamma"]*E0;
        ringPara["f0"]=(c*sqrt(1-1/ringPara["Gamma"]/ringPara["Gamma"]))/(2*M_PI*ringPara["R"]);
        ringPara["T0"]=1/ringPara["f0"];
        ringPara["omg0"]=2*M_PI*ringPara["f0"];
        ringPara["eta"]=1/ringPara["GMTSQ"]-1/(ringPara["Gamma"]*ringPara["Gamma"]);

        rfPara["Qs"].push_back(sqrt(rfPara["h"][0]*rfPara["Vref_I"][0]*ringPara["eta"]*fabs(cos((rfPara["Phi0"][0]-90)/180*M_PI))/(2*M_PI*ringPara["Ek"])));
		bunchPara["epsln"]=bunchPara["A"]/6.0;
        bunchPara["delta_hat"]=sqrt(bunchPara["epsln"])*sqrt(ringPara["omg0"]/(M_PI*ringPara["Ek"]))*std::pow((rfPara["h"][0]*rfPara["Vref_I"][0]*fabs(cos((rfPara["Phi0"][0]-90)/180*M_PI))/(2*M_PI*ringPara["Ek"]*ringPara["eta"])),0.25);
        bunchPara["t_hat"]=bunchPara["delta_hat"]/rfPara["Qs"][0]*ringPara["eta"]/ringPara["omg0"];
		//std::cout<<"Read parameters in."<<std::endl;
		return 0;
    };
    int printout(){
        std::cout<<"General Parameters:"<<std::endl;
        for(auto& x:generalPara){
            std::cout<<x.first<<"="<<x.second<<std::endl;
        }
        std::cout<<"Ring Parameters:"<<std::endl;
        for(auto& x:ringPara){
            std::cout<<x.first<<"="<<x.second<<std::endl;
        }
        std::cout<<"Rf Parameters:"<<std::endl;
        for(auto& x:rfPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
		for(auto& x:apPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
        std::cout<<"Modulation Parameters:"<<std::endl;
        for(auto& x:modPara){
            std::cout<<x.first;
            for(auto&data:x.second){
                std::cout<<":"<<data;
            }
            std::cout<<std::endl;
        }
        std::cout<<"Bunch Parameters:"<<std::endl;
        for(auto& x:bunchPara){
            std::cout<<x.first<<"="<<x.second<<std::endl;
        }
        return 0;
    };
};

#endif