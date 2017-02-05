#include <iostream>
#include<libint2.hpp>
#include<string.h>
#include<vector>
#include<fstream>
#include<Eigen/Dense>

using namespace std;
using namespace libint2;
using Eigen::MatrixXd;

int overlap (BasisSet obs);
int kinetic (BasisSet obs);
int main(int argc, char *argv[]) {
    
    	libint2::initialize();  // safe to use libint now	
	int i,j,k;
//path to xyz stored in string
	string h2 = "/Users/kg/HF/hf/h2.xyz";
//reads file at path and stores in variable
	ifstream ip (h2);
//cartesian coords read from file and stored in vector
	vector<Atom> atoms = read_dotxyz(ip);
//init basis set obs and assoc with sto-3g basis predet in libint
	BasisSet obs("STO-3G",atoms);
//print basis set
//	std::copy(begin(obs),end(obs),std::ostream_iterator<Shell>(std::cout,"\n"));
	cout << "Overlap" << endl;
	overlap (obs);
	cout << "Kinetic" << endl;
	kinetic (obs);
	libint2::finalize();  // do not use libint after this 
	return 0;
}
/*---------------------*/

int overlap (BasisSet obs) {
//Overlap Integral Engine::enum class Operator that spec set of one of more operators

		Engine s_engine(Operator::overlap,  // will compute overlap ints
                	obs.max_nprim(),    // max # of primitives in shells this engine will accept
                	obs.max_l()         // max angular momentum of shells this engine will accept
              	 );
		auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
                                // shell2bf[0] = index of the first basis function in shell 0
                                // shell2bf[1] = index of the first basis function in shell 1
                                // ...
		const auto& buf_vec = s_engine.results(); // will point to computed shell sets
                                          // const auto& is very important!

		for(auto s1=0; s1!=obs.size(); ++s1) {
  			for(auto s2=0; s2!=obs.size(); ++s2) {
    				cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    				s_engine.compute(obs[s1], obs[s2]);
    				cout << "done" << endl;
    				auto ints_shellset = buf_vec[0];  // location of the computed integrals
    				if (ints_shellset == nullptr)
      					continue;  // nullptr returned if the entire shell-set was screened out

    				auto bf1 = shell2bf[s1];  // first basis function in first shell
    				auto n1 = obs[s1].size(); // number of basis functions in first shell
    				auto bf2 = shell2bf[s2];  // first basis function in second shell
    				auto n2 = obs[s2].size(); // number of basis functions in second shell
    				// integrals are packed into ints_shellset in row-major (C) form
    				// this iterates over integrals in this order
				MatrixXd matS(2,2);
    				for(auto f1=0; f1!=n1; ++f1)
      					for(auto f2=0; f2!=n2; ++f2)
        					cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
					//	cout << "" << matS(f1,f2) << endl;
			}
		}

	return 0;
}
/*-------------*/
int kinetic (BasisSet obs) {
//KINETIC ENERGY OPERATOR
		Engine k_engine(Operator::kinetic,
                	obs.max_nprim(),
                	obs.max_l()
              	 );
		auto shell2bf = obs.shell2bf();
		const auto& buf_vec = k_engine.results(); // will point to computed shell sets
                                          // const auto& is very important!

		for(auto s1=0; s1!=obs.size(); ++s1) {
  			for(auto s2=0; s2!=obs.size(); ++s2) {
    				cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    				k_engine.compute(obs[s1], obs[s2]);
    				cout << "done" << endl;
    				auto ints_shellset = buf_vec[0];  // location of the computed integrals
    				if (ints_shellset == nullptr)
      					continue;  // nullptr returned if the entire shell-set was screened out

    				auto bf1 = shell2bf[s1];  // first basis function in first shell
    				auto n1 = obs[s1].size(); // number of basis functions in first shell
    				auto bf2 = shell2bf[s2];  // first basis function in second shell
    				auto n2 = obs[s2].size(); // number of basis functions in second shell
    				// integrals are packed into ints_shellset in row-major (C) form
    				// this iterates over integrals in this order
				MatrixXd matS(2,2);
    				for(auto f1=0; f1!=n1; ++f1)
      					for(auto f2=0; f2!=n2; ++f2)
        					cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
		
  			}			
		}	
	return 0;
}
