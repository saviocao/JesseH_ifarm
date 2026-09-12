#if !defined(HIST3D)
#define HIST3D

#include "TROOT.h"
#include "TFile.h"
#include "TH3.h"

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class Hist3D : public UserAmplitude< Hist3D >
{
    
public:
	
	Hist3D() : UserAmplitude< Hist3D >() { };
	Hist3D( const vector< string >& args );
	
	string name() const { return "Hist3D"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
private:
	
        string fileName, histName, histType, particleList;
	TH3 *hist3D;
};

#endif
