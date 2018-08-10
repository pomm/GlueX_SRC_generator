#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom2.h"

void Gen_gn_pimp_carbon(int N = 100000, double Ndays = 10, bool printOutput = true){

  TString out;
  out.Form("gn_pimp_N%d_carbon_test_th40_140.root",N);

  // masses:                                                                                                                                                    
  const Double_t md = 1.875612;
  const Double_t mp = 0.938272;
  const Double_t mn=  0.939565;
  const Double_t mpi = 0.13957018;
  const Double_t mC12 = 11.1750;
  const Double_t mB10 = 9.3245;
  int A = 12;

  // Half of the step in theta_cm                                                                                                                            
  double tsh = 1./2.;// degrees                                                                                                                             
  double theta_min = 40.;
  double theta_max = 140.;
  int Nthpo = (theta_max - theta_min)/(2.*tsh);

  //parameters for output:                                                                 
  int Runno = 9001;
  int NofFSparticles = 2;
  int Nevent = 0;

  // output tree
  TFile* f = new TFile(out,"recreate");
  TTree* T = new TTree("T","T");

  double theta_recoil, phi_recoil, Precoil;
  double theta_miss, phi_miss, Pmiss, Pmiss_rec;
  double theta_P3, phi_P3, absP3;
  double theta_P4, phi_P4, absP4;
  TVector3 P3, P4;
  double t,u,s, s_init, s_check, k_i, k_f;
  double theta_cm, weight_kk;
  double Original_E, Effective_E;
  double gamma_cm, Nucleon_cm;
  double cross_section;
  double E_beam;

  T->Branch("weight_kk",&weight_kk,"weight_kk/D");
  T->Branch("theta_cm",&theta_cm,"theta_cm/D");
  T->Branch("Effective_E",&Effective_E,"Effective_E/D");   
  T->Branch("Original_E",&Original_E,"Original_E/D");      
  T->Branch("gamma_cm",&gamma_cm,"gamma_cm/D");            
  T->Branch("Nucleon_cm",&Nucleon_cm,"Nucleon_cm/D");      

  T->Branch("k_i",&k_i,"k_i/D");
  T->Branch("k_f",&k_f,"k_f/D");
  T->Branch("s_init",&s_init,"s_init/D");

  T->Branch("cross_section",&cross_section,"cross_section/D");
  T->Branch("theta_recoil",&theta_recoil,"theta_recoil/D");
  T->Branch("phi_recoil",&phi_recoil,"phi_recoil/D");
  T->Branch("Precoil",&Precoil,"Precoil/D");         
  T->Branch("theta_miss",&theta_miss,"theta_miss/D");
  T->Branch("phi_miss",&phi_miss,"phi_miss/D");
  T->Branch("P_miss",&Pmiss,"Pmiss/D");        
  T->Branch("P_miss_rec",&Pmiss_rec,"Pmiss_rec/D");
  T->Branch("theta_P3",&theta_P3,"theta_P3/D");
  T->Branch("phi_P3",&phi_P3,"phi_P3/D");      
  T->Branch("absP3",&absP3,"absP3/D");         
  T->Branch("theta_P4",&theta_P4,"theta_P4/D");
  T->Branch("phi_P4",&phi_P4,"phi_P4/D");   
  T->Branch("absP4",&absP4,"absP4/D");      
  T->Branch("P3","TVector3",&P3);
  T->Branch("P4","TVector3",&P4);
  T->Branch("t",&t,"t/D");
  T->Branch("u",&u,"u/D");
  T->Branch("s",&s,"s/D");

  // Parametrization of the momentum distribution for nucleons in a nucleus (including SRC nucleons) 
  TF1* SRCtail = new TF1("SRCtail","1/(x**4)",0.25,0.7);
  TF1* n_k     = new TF1("","((x<0.25)?0.8*3/((0.25)**3):0.5*0.2*2.533333/(1./(0.25)-1./5)/(x**4))       ",0,1);
  TF1* n_k_k2  = new TF1("","((x<0.25)?0.8*3/((0.25)**3):0.5*0.2*2.533333/(1./(0.25)-1./5)/(x**4))*(x**2)",0,1);
  TF1* R_2 = new TF1("","x**2",0,0.25);

  TVector3 Prec3, Pcm3, Pmiss3;
  double con, Sr, Beam, cs_theta_cm;
  double PLCrecoil, PLCcm;

  // variables for the reaction particles in different c.s.:                                                                         
  TVector3 m1, m; // boost vectors                                                                                         
  Double_t rot_phi, rot_theta, theta_ll, phi_ll;
  TLorentzVector Vbeam_nrf, Vmiss_nrf;
  TLorentzVector Vbeam_cm, Vmiss_cm;
  TVector3 Vbeam_cm_oz3;
  TVector3 Vmiss_cm_oz3;
  TVector3 Vpart3_cm_oz3, Vpart4_cm_oz3;
  Double_t Epart3_cm, Epart4_cm;
  TVector3 Vpart3_cm3, Vpart4_cm3;
  Double_t PmissX, PmissY, PmissZ;
  Double_t PrecX, PrecY, PrecZ;
  Double_t Ecm;

  for (int i=0; i<N; i++){

    //       if(i%(N/100)==0) cout << i/(N/100)<<"%" << endl; 
    E_beam = 8.5;
    Original_E = E_beam;

    //
    // struck neutron
    //
    Pmiss = n_k_k2->GetRandom();
    gRandom->Sphere(PmissX, PmissY, PmissZ, Pmiss);
    theta_miss = TMath::ACos(PmissZ/Pmiss)*TMath::RadToDeg();  
    phi_miss = TMath::ATan2(PmissY,PmissX)*TMath::RadToDeg();
    Pmiss3.SetXYZ(PmissX, PmissY, PmissZ);

    //
    // recoil proton
    //
    if(Pmiss > 0.25){
      Pcm3.SetXYZ(gRandom->Gaus(0,0.14),gRandom->Gaus(0,0.14),gRandom->Gaus(0,0.14));
      Prec3 = Pcm3 - Pmiss3;
      Precoil = Prec3.Mag();
    }else{
      Precoil = R_2->GetRandom();
      gRandom->Sphere(PrecX, PrecY, PrecZ, Precoil);
      Prec3.SetXYZ(PrecX, PrecY, PrecZ);
      Pcm3 = Prec3 + Pmiss3;
    }
    Ecm = mC12 - TMath::Sqrt( pow(Pcm3.Mag(),2) + pow(mB10,2) );
    theta_recoil = Prec3.Theta()*TMath::RadToDeg();
    phi_recoil = Prec3.Phi()*TMath::RadToDeg();

    //    cout<<"Pmiss = "<<Pmiss<<", Prec3.Theta - "<<Prec3.Theta()<<", theta_recoil"<<theta_recoil<<", PrecX = "<<PrecX<<", PrecY = "<<PrecY<<", PrecZ = "<<PrecZ<<endl;

    // light cone momentum fractions:                                                                      
    PLCrecoil = (sqrt(mn*mn + pow(Precoil,2)) - Prec3.Z())/mn;
    PLCcm = (2*mn - Pcm3.Z())/mn;

    // Factor:
    s_init = 0.;
    if(Pmiss > 0.75 || (PLCcm - PLCrecoil) < 0.) continue;
    if(Pmiss < 0.25){
      s_init = pow( E_beam + mn - 0.015 - Pmiss*Pmiss/(2.*mn*(A-1)) ,2) - pow(PmissX,2) - pow(PmissY,2) - pow(PmissZ+E_beam,2);
    }else if(Pmiss >= 0.25){
      s_init = pow( E_beam + (mn+mp) - sqrt(pow(Precoil,2) + mp*mp) ,2) - pow(PmissX,2) - pow(PmissY,2) - pow(PmissZ+E_beam,2);
    }
    if(s_init < pow(mn+mpi,2)) continue;

    // k_i, k_f relative c.m. momentum in the initial and final states:                            
    k_i = (s_init - 0.940*0.940)/2./sqrt(s_init);
    k_f = sqrt((s_init - pow((0.938-0.140),2))*(s_init - pow((0.938+0.140),2))/4./s_init);

    // 4-momentum of the beam and nucleon in the lab:                                               
    TLorentzVector Vbeam_lab(TVector3(0.,    0.,    E_beam), E_beam                               );
    TLorentzVector Vmiss_lab(Pmiss3                        , Ecm - sqrt(pow(mp,2)+pow(Precoil,2)) );
    s_check = (Vbeam_lab + Vmiss_lab)*(Vbeam_lab + Vmiss_lab);

    // boost to the c.m. system of the nucleon + beam photon:                                
    m = (Vmiss_lab + Vbeam_lab).BoostVector();
    Vbeam_cm = Vbeam_lab;
    Vmiss_cm = Vmiss_lab;
    Vbeam_cm.Boost(-m);
    Vmiss_cm.Boost(-m);
    gamma_cm = Vbeam_cm.Vect().Mag();
    Nucleon_cm = Vmiss_cm.Vect().Mag();

    for(int k=0; k<=Nthpo; k++){

      theta_cm = theta_min+k*tsh*2.;
      // rotate vectors so that the collision is along z axis:                                                                   
      rot_phi = Vbeam_cm.Vect().Phi();
      rot_theta = Vbeam_cm.Vect().Theta();
      Vbeam_cm_oz3 = Vbeam_cm.Vect();
      Vmiss_cm_oz3 = Vmiss_cm.Vect();
      Vbeam_cm_oz3.RotateZ(-rot_phi);
      Vbeam_cm_oz3.RotateY(-rot_theta);
      Vmiss_cm_oz3.RotateZ(-rot_phi);
      Vmiss_cm_oz3.RotateY(-rot_theta);

      // scattering angle - theta_cm                                                           
      theta_ll = theta_cm*TMath::DegToRad();
      phi_ll = gRandom->Uniform(0,2.*TMath::Pi());

      // outgoing particles:               
      Vpart3_cm_oz3.SetXYZ( Vbeam_cm_oz3.Z()*cos(phi_ll)*sin(theta_ll),  Vbeam_cm_oz3.Z()*sin(phi_ll)*sin(theta_ll),  Vbeam_cm_oz3.Z()*cos(theta_ll));
      Vpart4_cm_oz3.SetXYZ(-Vbeam_cm_oz3.Z()*cos(phi_ll)*sin(theta_ll), -Vbeam_cm_oz3.Z()*sin(phi_ll)*sin(theta_ll), -Vbeam_cm_oz3.Z()*cos(theta_ll));
      // make sure the outgoing particle are on mass shell
      Epart3_cm = (s_check - 0.938*0.938 + 0.140*0.140)/(2.*sqrt(s_check));
      Epart4_cm = (s_check + 0.938*0.938 - 0.140*0.140)/(2.*sqrt(s_check));

      // rotate back:                                           
      Vpart3_cm3 = Vpart3_cm_oz3;
      Vpart4_cm3 = Vpart4_cm_oz3;
      Vpart3_cm3.RotateY(rot_theta);
      Vpart3_cm3.RotateZ(rot_phi);
      Vpart4_cm3.RotateY(rot_theta);
      Vpart4_cm3.RotateZ(rot_phi);

      // boost outgoing particles to the lab system:                     
      TLorentzVector Vpart3_lab(Vpart3_cm3, Epart3_cm);//pion
      TLorentzVector Vpart4_lab(Vpart4_cm3, Epart4_cm);//proton
      Vpart3_lab.Boost(m);
      Vpart4_lab.Boost(m);

      // save stuff to the tree:      
      Pmiss_rec = (Vpart3_lab.Vect() + Vpart4_lab.Vect() - Vbeam_lab.Vect()).Mag();                                                                                
      P3 = Vpart3_lab.Vect();
      P4 = Vpart4_lab.Vect();
      absP3 = Vpart3_lab.Vect().Mag();
      absP4 = Vpart4_lab.Vect().Mag();
      theta_P3 = P3.Theta()*TMath::RadToDeg();
      theta_P4 = P4.Theta()*TMath::RadToDeg();
      phi_P3 = P3.Phi()*TMath::RadToDeg();
      phi_P4 = P4.Phi()*TMath::RadToDeg();

      s = (Vpart3_lab + Vpart4_lab)*(Vpart3_lab + Vpart4_lab);
      t = (Vbeam_lab - Vpart3_lab)*(Vbeam_lab - Vpart3_lab);
      u = (Vbeam_lab - Vpart4_lab)*(Vbeam_lab - Vpart4_lab);

      // define components for the weight:                                                    
      //    Target       N days       nb->b   b->cm2                         
      con = 1.45e23 * 3600*24*Ndays * 1e-9 * 1e-24; // target carbon 1.45e23 atoms/cm2, assume target thickness 7% X0=1.3 cm 
      Sr = (2*3.14*(cos((theta_cm-tsh)*3.14/180)-cos((theta_cm+tsh)*3.14/180))) ; // theta_cm +- half of the step                  
      Beam = 2.e7;
      cs_theta_cm = pow((1-cos(theta_cm*3.14/180)),-5) * pow((1+cos(theta_cm*3.14/180)),-4);

      weight_kk = con * (1.25*1.e7*pow(s_init,-7.)*(k_i*k_f)/3.14) * Sr * Beam * cs_theta_cm * (A/2.) * pow(A, -1./3.) / N;

      if(abs(t)>1. && abs(u)>1.){ // only consider 'hard' scattering events.                                  
	T->Fill();
      }

      if(printOutput && abs(t) > 2 && abs(u) > 2){
        if(gRandom->Uniform(0,1) <= weight_kk){
          Nevent++;
          if(Precoil == 0){
            printf("%d %d %d\n",Runno, Nevent, NofFSparticles);
            printf("%d %d %lf\n", 1, 9, 0.140);
            printf("   %d %lf %lf %lf %lf\n",-1, P3.x(), P3.y(), P3.z(), sqrt(pow(absP3,2) + pow(0.140,2)));
            printf("%d %d %lf\n", 2, 14, 0.938);
            printf("   %d %lf %lf %lf %lf\n", 1, P4.x(), P4.y(), P4.z(), sqrt(pow(absP4,2) + pow(0.938,2)));
          }else{
            printf("%d %d %d\n",Runno, Nevent, NofFSparticles+1);
            printf("%d %d %lf\n", 1, 9, 0.140);
            printf("   %d %lf %lf %lf %lf\n",-1, P3.x(), P3.y(), P3.z(), sqrt(pow(absP3,2) + pow(0.140,2)));
            printf("%d %d %lf\n", 2, 14, 0.938);
            printf("   %d %lf %lf %lf %lf\n", 1, P4.x(), P4.y(), P4.z(), sqrt(pow(absP4,2) + pow(0.938,2)));
            printf("%d %d %lf\n", 3, 14, 0.938);
            printf("   %d %lf %lf %lf %lf\n",0, Prec3.X(), Prec3.Y(), Prec3.Z(), sqrt(pow(Precoil,2) + pow(0.938,2)));
          }
        }
      }
    } // for theta_cm
  } // for N
  T->Write();
  f->Write();
  f->Close();
}
