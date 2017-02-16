#include <iostream>
#include<libint2.hpp>
#include<string.h>
#include<vector>
#include<fstream>
#include<Eigen/Dense>

using namespace std;
using namespace libint2;
using Eigen::MatrixXd;

MatrixXd overlap (BasisSet obs);
MatrixXd kinetic (BasisSet obs);
MatrixXd nuclearAtt (BasisSet obs, vector<Atom> atoms);
MatrixXd coulomb (BasisSet obs, MatrixXd density);
int main() {
    
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
	cout << "Basis Set" << endl;
	std::copy(begin(obs),end(obs),std::ostream_iterator<Shell>(std::cout,"\n"));
	cout << "Overlap" << endl;
	MatrixXd S;
	MatrixXd density = MatrixXd::Zero(obs.nbf(),obs.nbf());
	S = overlap (obs);
	cout << S << endl; 

	cout << "Kinetic" << endl;
	MatrixXd Kin;
	Kin = kinetic (obs);
	cout << Kin << endl;
 
	cout << "Nuclear Attraction" << endl;
	MatrixXd Nuc;
	Nuc = nuclearAtt (obs, atoms);
	cout << Nuc << endl;

	cout << "Coulomb" << endl;
	MatrixXd Coul;
	Coul = coulomb(obs,density);
	cout << Coul << endl; 
 
	libint2::finalize();  // do not use libint after this 
	return 0;
}
/*---------------------*/

MatrixXd overlap (BasisSet obs) {
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
//define eigen matrix
	MatrixXd S = MatrixXd::Zero(obs.nbf(),obs.nbf());
		for(auto s1=0; s1!=obs.size(); ++s1) {
  			for(auto s2=0; s2!=obs.size(); ++s2) {
  //  				cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    				s_engine.compute(obs[s1], obs[s2]);
  //  				cout << "done" << endl;
    				auto ints_shellset = buf_vec[0];  // location of the computed integrals
    				if (ints_shellset == nullptr)
      					continue;  // nullptr returned if the entire shell-set was screened out

    				auto bf1 = shell2bf[s1];  // first basis function in first shell
    				auto n1 = obs[s1].size(); // number of basis functions in first shell
    				auto bf2 = shell2bf[s2];  // first basis function in second shell
    				auto n2 = obs[s2].size(); // number of basis functions in second shell
    				// integrals are packed into ints_shellset in row-major (C) form
    				// this iterates over integrals in this order
    				for(auto f1=0; f1!=n1; ++f1) {
      					for(auto f2=0; f2!=n2; ++f2) {
    //    					cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
						S(bf1+f1,bf2+f2) = ints_shellset[f1*n2+f2];	
					}	
				}	
			}
		}

	return S;
}
/*-------------*/
MatrixXd kinetic (BasisSet obs) {
//KINETIC ENERGY OPERATOR
		Engine k_engine(Operator::kinetic,
                	obs.max_nprim(),
                	obs.max_l()
              	 );
		auto shell2bf = obs.shell2bf();
		const auto& buf_vec = k_engine.results(); // will point to computed shell sets
                                          // const auto& is very important!

	MatrixXd kin = MatrixXd::Zero(obs.nbf(),obs.nbf());
		for(auto s1=0; s1!=obs.size(); ++s1) {
  			for(auto s2=0; s2!=obs.size(); ++s2) {
    //				cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    				k_engine.compute(obs[s1], obs[s2]);
    //				cout << "done" << endl;
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
      					for(auto f2=0; f2!=n2; ++f2) {
      //  					cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
						kin(bf1+f1,bf2+f2) = ints_shellset[f1*n2+f2];	
					}	
  			}			
		}	
	return kin;
}

MatrixXd nuclearAtt (BasisSet obs, vector<Atom> atoms) {
//Engine::set_params();
// this engine will compute nuclear attraction ints
// by default there are no charges
Engine v_engine(Operator::nuclear,
                obs.max_nprim(), obs.max_l());
v_engine.set_params(make_point_charges(atoms));  // convert `atoms` to point charges
                                                 // classical charges in QM/MM need extra work
		auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
                                // shell2bf[0] = index of the first basis function in shell 0
                                // shell2bf[1] = index of the first basis function in shell 1
                                // ...
		const auto& buf_vec = v_engine.results(); // will point to computed shell sets
                                          // const auto& is very important!
//define eigen matrix
	MatrixXd S = MatrixXd::Zero(obs.nbf(),obs.nbf());
		for(auto s1=0; s1!=obs.size(); ++s1) {
  			for(auto s2=0; s2!=obs.size(); ++s2) {
  //  				cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    				v_engine.compute(obs[s1], obs[s2]);
  //  				cout << "done" << endl;
    				auto ints_shellset = buf_vec[0];  // location of the computed integrals
    				if (ints_shellset == nullptr)
      					continue;  // nullptr returned if the entire shell-set was screened out

    				auto bf1 = shell2bf[s1];  // first basis function in first shell
    				auto n1 = obs[s1].size(); // number of basis functions in first shell
    				auto bf2 = shell2bf[s2];  // first basis function in second shell
    				auto n2 = obs[s2].size(); // number of basis functions in second shell
    				// integrals are packed into ints_shellset in row-major (C) form
    				// this iterates over integrals in this order
    				for(auto f1=0; f1!=n1; ++f1) {
      					for(auto f2=0; f2!=n2; ++f2) {
    //    					cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
						S(bf1+f1,bf2+f2) = ints_shellset[f1*n2+f2];	
					}	
				}	
			}
		}

	return S;
}

/*--------------------------*/
MatrixXd coulomb (BasisSet obs, MatrixXd density) {
//Coulomb integrals
		Engine engine(Operator::coulomb,
                	obs.max_nprim(),
                	obs.max_l()
              	 );
		auto shell2bf = obs.shell2bf();
	const auto& buf_vec = engine.results(); // will point to computed shell sets
                                           // const auto& is very important!
	MatrixXd S = MatrixXd::Zero(obs.nbf(),obs.nbf());
	int cnt=0;
for(auto s1=0; s1!=obs.size(); ++s1) {
  for(auto s2=0; s2!=obs.size(); ++s2) {
    for(auto s3=0; s3!=obs.size(); ++s3) {
      for(auto s4=0; s4!=obs.size(); ++s4) {
//    cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
//    cout << "done" << endl;

    auto bf1 = shell2bf[s1];  // first basis function in first shell
    auto n1 = obs[s1].size(); // number of basis functions in first shell
    auto bf2 = shell2bf[s2];  // first basis function in second shell
    auto n2 = obs[s2].size(); // number of basis functions in second shell
    auto bf3 = shell2bf[s3];  // first basis function in first shell
    auto n3 = obs[s3].size(); // number of basis functions in first shell
    auto bf4 = shell2bf[s4];  // first basis function in first shell
    auto n4 = obs[s4].size(); // number of basis functions in first shell
 
    auto shellset = buf_vec[0];
    for(auto f1=0; f1!=n1; ++f1)
	for(auto f2=0; f2!=n2; ++f2)
	     for(auto f3=0; f3!=n3; ++f3)
      		for(auto f4=0; f4!=n4; ++f4) {
//      cout << "  " << shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4]  << endl;
	S(bf3+f3,bf4+f4) = shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4] * density(bf3+f3,bf4+f4);
	}	
  }
}
}
}	
return S;
}


