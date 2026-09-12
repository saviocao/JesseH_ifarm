
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Hist3D.h"

Hist3D::Hist3D( const vector< string >& args ) :
UserAmplitude< Hist3D >( args )
{
	assert( args.size() == 4 );
	fileName = args[0].c_str();
	histName = args[1].c_str();
	histType = args[2].c_str();
	particleList = args[3].c_str();

	cout<<"Opening ROOT file "<<fileName.data()<<endl;
	cout<<"Model provided in histogram named "<<histName.data()<<endl;
	cout<<"Histogram type for generator "<<histType.data()<<endl;
	cout<<"Summing particle indices "<<particleList.data()<<" for invariant mass"<<endl;
	for(uint i=0; i<particleList.length(); i++) {
		string num; num += particleList[i];
		int index = atoi(num.c_str());
		cout<<index<<endl;
	}

	TFile *finput = TFile::Open(fileName.data());
	if(!finput->IsOpen()) {
		cout<<"Can't find file "<<fileName.data()<<endl;
		exit(1);
	}

	hist3D = (TH3*)finput->Get(histName.data());
	if(!hist3D) {
		cout<<"Can't find histogram "<<histName.data()<<" in file "<<fileName.data()<<endl;
		exit(1);
	}		

	if( histType != "ResMassVsCosThetaVsBeamE" && histType != "ResMassVsCosThetaVsT") {
		cout<<"Type of 3D histogram is not currently supported, please add necessary kinematics and options to Hist3D amplitude"<<endl;
		exit(1);
	}

	// keep in memory after file is closed
        hist3D->SetDirectory(gROOT);
	finput->Close();
}


complex< GDouble >
Hist3D::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector target  ( 0., 0., 0., 0.938);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1 ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2 ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 

	// compute particle P4 sum for invariant mass
	TLorentzVector res;
	for(uint i=0; i<particleList.length(); i++) {
		string num; num += particleList[i];
		int index = atoi(num.c_str());
		TLorentzVector particleP4 ( pKin[index][1], pKin[index][2], pKin[index][3], pKin[index][0] ); 
		res += particleP4;
	}
	
	double beamE = beam.E();
	double t = fabs((beam - recoil).M2()); 

	//COM frame
	TLorentzVector com = beam + target;
	TLorentzRotation comRestBoost ( -com.BoostVector() );

	TLorentzVector beam_com   = comRestBoost * beam;
	TLorentzVector recoil_com = comRestBoost * recoil;
	TLorentzVector res_com = comRestBoost * res;
	TLorentzVector p1_com = comRestBoost * p1;
	TLorentzVector p2_com = comRestBoost * p2;
	
	//Resonance Rest Frame from COM
	TLorentzRotation resRestBoost( -res_com.BoostVector() );
	TLorentzVector beam_res   = resRestBoost * beam_com;
	TLorentzVector recoil_res = resRestBoost * recoil_com;
	TLorentzVector p1_res = resRestBoost * p1_com;
	TLorentzVector p2_res = resRestBoost * p2_com;
	
	// helicity frame: z-axis is propagation of resonance X => opposite recoil proton in X rest frame
	TVector3 z = -1 * res_com.Vect().Unit();

	// y axis perpendicular to production plane
	TVector3 y = beam_com.Vect().Cross(z).Unit();

	TVector3 x = y.Cross(z);
  
	TVector3 angles( (p2_res.Vect()).Dot(x),
			 (p2_res.Vect()).Dot(y),
			 (p2_res.Vect()).Dot(z) );
  
	double cosTheta = angles.CosTheta();
	//double phi = angles.Phi();

	double userVarX = 0;
	double userVarY = 0;
	double userVarZ = 0;

	if(histType == "ResMassVsCosThetaVsBeamE") {
	        userVarZ = beamE;
	        userVarY = cosTheta;
		userVarX = res.M();
	}	


	if(histType == "ResMassVsCosThetaVsT") {
	        userVarZ = t;
	        userVarY = cosTheta;
		userVarX = res.M();
	}	

	// weighted model of intensity from histogram 
	GDouble W = 0.; // initialized to zero

	int bin = hist3D->FindBin(userVarX, userVarY, userVarZ); // generic bin index from 2D histogram (negative value if values outside defined range)
	if(bin > 0) W = hist3D->GetBinContent(bin); 

	return complex< GDouble > ( sqrt(W) );
}
