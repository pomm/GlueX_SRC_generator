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

void Gen_gn_pimp_deu(int N = 10000, double Ndays = 10, bool printOutput = true){

  // Half of the step in theta_cm
  double tsh = 1./2.;// degrees                                                                                                      
  double theta_min = 40.;//40.; it used to be 30 - 150 - in the proposal
  double theta_max = 140.;//130.;
  int Nthpo = (theta_max - theta_min)/(2.*tsh);

  //deuteron
  int A = 2;

  // masses:
  const Double_t md = 1.875612; // deuteron
  const Double_t mp = 0.938272;
  const Double_t mn=  0.939565;
  const Double_t mpi = 0.13957018;
  const Double_t mC12 = 11.1750;
  const Double_t mB10 = 9.3245;

  //parameters for output:                                           
  int Runno = 10000;
  int NofFSparticles = 2;
  int Nevent = 0;

  TString out;
  out.Form("gn_pimp_N%d_deu_test3_th40_140.root",N);

  double E_beam = 9;

  TFile* f = new TFile(out,"recreate");
  TTree* T = new TTree("T","T");

  double theta_recoil, phi_recoil, Precoil;
  double theta_miss, phi_miss, Pmiss, Pmiss_rec;
  double theta_P3, phi_P3, absP3;
  double theta_P4, phi_P4, absP4;
  TVector3 P3, P4;
  double t,u,s, s_init, s_check;
  double theta_cm;
  double weight_kk = 0.;
  double Effective_E;
  double Original_E;
  double gamma_cm;
  double Nucleon_cm;
  double cross_section;

  double k_i;
  double k_f;

  T->Branch("weight_kk",&weight_kk,"weight_kk/D");
  T->Branch("theta_cm",&theta_cm,"theta_cm/D");
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

  // Final 3-momentum vector for the scattered protons in the lab                                                                                                                    
  T->Branch("P3","TVector3",&P3);
  T->Branch("P4","TVector3",&P4);

  T->Branch("t",&t,"t/D");
  T->Branch("u",&u,"u/D");
  T->Branch("s",&s,"s/D");

  // Parametrization of the momentum distribution for nucleons in a nucleus (including SRC nucleons)
    TF1* n_k_k2  = new TF1("","((x<0.25)?0.8*3/((0.25)**3):0.5*0.2*2.533333/(1./(0.25)-1./5)/(x**4))*(x**2)",0,1);
  //TF1* n_k_k2 = new TF1("","( 157.4*exp(-1.24*((x/0.1973)**2))/((1+18.3*((x/0.1973)**2))**2)+0.234*exp(-1.27*((x/0.1973)**2))+0.00623*exp(-0.22*((x/0.1973)**2)) ) * (x**2)",0,1);

  double PmissX, PmissY, PmissZ;
  TVector3 Prec3, Pcm3;
  double con, Sr, Beam, cs_theta_cm;
  double PLCrecoil, PLCcm;

  // variables for the reaction particles in different c.s.:
  TVector3 m; // boost vectors
  Double_t rot_phi, rot_theta, theta_ll, phi_ll;
  TLorentzVector Vbeam_nrf;
  TLorentzVector Vmiss_nrf;
  TLorentzVector Vbeam_cm;
  TLorentzVector Vmiss_cm;
  TLorentzVector m1;
  TVector3 Vbeam_cm_oz3;
  TVector3 Vmiss_cm_oz3;
  TVector3 Vpart3_cm_oz3;
  TVector3 Vpart4_cm_oz3;
  Double_t Epart3_cm, Epart3_cm1;
  Double_t Epart4_cm, Epart4_cm1;
  TVector3 Vpart3_cm3;
  TVector3 Vpart4_cm3;

  // Event generation:
  for (int i=0; i<N; i++){

    //    if(i%(N/100)==0) cout << i/(N/100)<<"%" << endl; 

    Original_E = E_beam;

    Pmiss = n_k_k2->GetRandom();
    gRandom->Sphere(PmissX, PmissY, PmissZ, Pmiss);
    theta_miss = TMath::ACos(PmissZ/Pmiss)*180/3.14;      // Calculate theta [degrees]                                                                                      
    phi_miss = TMath::ATan2(PmissY,PmissX)*180/3.14;      // Calculate phi   [degrees]

    Prec3.SetXYZ(-PmissX, -PmissY, -PmissZ);
    Precoil = Prec3.Mag();
    theta_recoil = Prec3.Theta()/TMath::Pi()*180.;
    phi_recoil = Prec3.Phi()/TMath::Pi()*180.;
    // light cone momentum fraction
    PLCrecoil = (sqrt(0.940*0.940 + pow(Precoil,2)) - Prec3.Z())/0.940;
    PLCcm = (2*0.94 - Pcm3.Z())/0.940;

    // Factor from Misak:
    s_init = 0.;
    // scattering off proton inside nucleus:
    // calculate light cone momentum fractions:
    if(Pmiss > 0.75 || (PLCcm - PLCrecoil) < 0.) continue;
    if(Pmiss < 0.25){
      s_init = pow( E_beam + 0.940 - 0.015 - Pmiss*Pmiss/(2.*0.940*(A-1)) ,2) - pow(PmissX,2) - pow(PmissY,2) - pow(PmissZ+E_beam,2);
    }else if(Pmiss >= 0.25){
      s_init = pow( E_beam + 2*0.940 - sqrt(pow(Precoil,2) + 0.940*0.940) ,2) - pow(PmissX,2) - pow(PmissY,2) - pow(PmissZ+E_beam,2);
    }
    if(s_init < pow(0.94+0.14,2)) continue;

    // k_i, k_f relative c.m. momentum in the initial and final states:
    k_i = (s_init - 0.940*0.940)/2./sqrt(s_init);//sqrt(0.5*E_beam*0.940);
    k_f = sqrt((s_init - pow((0.938-0.140),2))*(s_init - pow((0.938+0.140),2))/4./s_init);

    // 4-momentum of the beam and nucleon in the lab:
    TLorentzVector Vbeam_lab(TVector3(0.,0.,E_beam),         E_beam                            );
    TLorentzVector Vmiss_lab(TVector3(PmissX,PmissY,PmissZ), md - sqrt(pow(Precoil,2) + 0.938*0.938));
    s_check = (Vbeam_lab + Vmiss_lab)*(Vbeam_lab + Vmiss_lab);

    // boost to the c.m. system of the nucleon + beam photon:
    m = (Vmiss_lab + Vbeam_lab).BoostVector();
    Vbeam_cm = Vbeam_lab;
    Vmiss_cm = Vmiss_lab;
    Vbeam_cm.Boost(-m);
    Vmiss_cm.Boost(-m);
    gamma_cm = Vbeam_cm.Vect().Mag();
    Nucleon_cm = Vmiss_cm.Vect().Mag();
    
    // loop over theta_cm scattering angles:
    for(int k=0; k<=Nthpo; k++){

      theta_cm = theta_min+k*tsh*2.;
      //        cout<<"theta = "<<theta_cm<<endl;

      // rotate vectors so that the collision is along z axis:
      rot_phi = Vbeam_cm.Vect().Phi();
      rot_theta = Vbeam_cm.Vect().Theta();
      Vbeam_cm_oz3 = Vbeam_cm.Vect();
      Vmiss_cm_oz3 = Vmiss_cm.Vect();
      Vbeam_cm_oz3.RotateZ(-rot_phi);
      Vbeam_cm_oz3.RotateY(-rot_theta);
      Vmiss_cm_oz3.RotateZ(-rot_phi);
      Vmiss_cm_oz3.RotateY(-rot_theta);

      //cout<<"rot_phi = "<<rot_phi<<", rot_theta = "<<rot_theta<<endl;
      //      cout<<"beam_cm: "<<Vbeam_cm_oz3.z()<<", mag = "<<Vbeam_cm_oz3.Mag()<<endl;
      //      cout<<"miss_cm: "<<Vmiss_cm_oz3.z()<<", mag = "<<Vmiss_cm_oz3.Mag()<<", angle with beam = "<<Vbeam_cm_oz3.Angle(Vmiss_cm_oz3)/3.1415*180.<<endl;
      
      // scattering angle - theta_cm
      theta_ll = theta_cm*TMath::Pi()/180.;
      phi_ll = gRandom->Uniform(0,2*TMath::Pi());
      
      // outgoing particles:
      Vpart3_cm_oz3.SetXYZ( Vbeam_cm_oz3.Z()*cos(phi_ll)*sin(theta_ll),  Vbeam_cm_oz3.Z()*sin(phi_ll)*sin(theta_ll),  Vbeam_cm_oz3.Z()*cos(theta_ll));
      Vpart4_cm_oz3.SetXYZ(-Vbeam_cm_oz3.Z()*cos(phi_ll)*sin(theta_ll), -Vbeam_cm_oz3.Z()*sin(phi_ll)*sin(theta_ll), -Vbeam_cm_oz3.Z()*cos(theta_ll));
      // make sure the outgoing particle are on mass shell
      //Epart3_cm = (pow(Vbeam_cm.E() + Vmiss_cm.E(),2) - 0.938*0.938 + 0.140*0.140)/2./(Vbeam_cm.E() + Vmiss_cm.E());
      //Epart4_cm = (pow(Vbeam_cm.E() + Vmiss_cm.E(),2) + 0.938*0.938 - 0.140*0.140)/2./(Vbeam_cm.E() + Vmiss_cm.E());//Vbeam_cm.E() + Vmiss_cm.E() - Epart3_cm;
      Epart3_cm = (s_check - 0.938*0.938 + 0.140*0.140)/(2.*sqrt(s_check));
      Epart4_cm = (s_check + 0.938*0.938 - 0.140*0.140)/(2.*sqrt(s_check));//Vbeam_cm.E() + Vmiss_cm.E() - Epart3_cm;

      // rotate back:
      Vpart3_cm3 = Vpart3_cm_oz3;
      Vpart4_cm3 = Vpart4_cm_oz3;
      Vpart3_cm3.RotateY(rot_theta);
      Vpart3_cm3.RotateZ(rot_phi);
      Vpart4_cm3.RotateY(rot_theta);
      Vpart4_cm3.RotateZ(rot_phi);
      
      // boost outgoing particles to the lab system:
      TLorentzVector Vpart3_lab(Vpart3_cm3, Epart3_cm);//sqrt( pow(Vpart3_cm3.Mag(),2) + pow(0.140,2) ) );
      TLorentzVector Vpart4_lab(Vpart4_cm3, Epart4_cm);//sqrt( pow(Vpart4_cm3.Mag(),2) + pow(0.938,2) ) );
      Vpart3_lab.Boost(m);
      Vpart4_lab.Boost(m);

      //cout<<"tmp3: "<<Vpart3_lab.X()<<", "<<Vpart3_lab.Y()<<" "<<Vpart3_lab.Z()<<endl;
      //cout<<"tmp4: "<<Vpart4_lab.X()<<", "<<Vpart4_lab.Y()<<" "<<Vpart4_lab.Z()<<endl;

      //      cout<<"E conservation: "<<E_beam + md - Vpart3_lab.E() - Vpart4_lab.E() - sqrt(pow(Precoil,2) + 0.938*0.938)<<endl;
      //cout<<"P conservation: "<<E_beam - (Vpart3_lab.Vect() + Vpart4_lab.Vect() + Prec3).Mag()<<endl;
      //cout<<"Pmiss = "<<Pmiss<<", Pmiss_rec = "<<(Vpart3_lab.Vect() + Vpart4_lab.Vect() - Vbeam_lab.Vect()).Mag()<<", in cm = "<<(Vpart3_cm3 + Vpart4_cm3 - Vbeam_cm.Vect()).Mag()<<endl;

      // save stuff to the tree:
      Pmiss_rec = (Vpart3_lab.Vect() + Vpart4_lab.Vect() - Vbeam_lab.Vect()).Mag();
      P3 = Vpart3_lab.Vect();
      P4 = Vpart4_lab.Vect();
      absP3 = Vpart3_lab.Vect().Mag();
      absP4 = Vpart4_lab.Vect().Mag();
      theta_P3 = P3.Theta()/TMath::Pi()*180.;
      theta_P4 = P4.Theta()/TMath::Pi()*180.;
      phi_P3 = P3.Phi()/TMath::Pi()*180.;
      phi_P4 = P4.Phi()/TMath::Pi()*180.;
  
      s = (Vpart3_lab + Vpart4_lab)*(Vpart3_lab + Vpart4_lab);
      t = (Vbeam_lab - Vpart3_lab)*(Vbeam_lab - Vpart3_lab);
      u = (Vbeam_lab - Vpart4_lab)*(Vbeam_lab - Vpart4_lab);

      // define components for the weight:
      //    Target       N days       nb->b   b->cm2
      con = 1.516e24  * 3600.*24. * Ndays * 1e-9 * 1e-24; // target deuterium N days!!! assume target thickness 30 cm = 4.1% X0
      Sr = (2*3.14*(cos((theta_cm-tsh)*3.14/180)-cos((theta_cm+tsh)*3.14/180))) ; // theta_cm +- half of the step
      Beam = 2e7;
      cs_theta_cm = pow((1-cos(theta_cm*3.14/180)),-5) * pow((1+cos(theta_cm*3.14/180)),-4);

      weight_kk = con * (1.25*1e7*pow(s_init,-7.)*(k_i*k_f)/3.14) * Sr * Beam * cs_theta_cm * (A/2.) * pow(A, -1./3.) / N;
      //      cout<<"w = "<<weight<<", w_kk = "<<weight_kk<<", weight_kks = "<<weight_kks<<endl;

      if(abs(t)>1. && abs(u)>1.){ // only consider 'hard' scattering events.
	T->Fill();
      }
      if(printOutput && abs(t) > 2. && abs(u) > 2.){
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

    } // loop over theta_cm
    //hga->Fill(gamma_cm*gamma_cm);

  }// loop over N

  T->Write();
  f->Write();
  f->Close();
}
