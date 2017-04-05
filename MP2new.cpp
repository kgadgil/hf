#include <iostream>
#include<libint2.hpp>
#include<string.h>
#include<vector>
#include<fstream>
#include<Eigen/Dense>
#include<unsupported/Eigen/CXX11/Tensor>

using namespace std;
using namespace libint2;
using Eigen::MatrixXd;
using Eigen::Tensor;

void printVec (vector <double> vec);
void printInt (vector <int> vec);
void printMat(vector<vector<double>> mat);
vector<double> tensor (BasisSet obs);
MatrixXd overlap (BasisSet obs);
MatrixXd kinetic (BasisSet obs);
MatrixXd nuclearAtt (BasisSet obs, vector<Atom> atoms);
MatrixXd coulomb (BasisSet obs, MatrixXd density);
MatrixXd exchange (BasisSet obs, MatrixXd density);
double hf(BasisSet obs, MatrixXd overlap, MatrixXd Core, int itr, vector<Atom> atoms, int occ, MatrixXd &evec, MatrixXd &eval);
int noddy (BasisSet obs, MatrixXd evec, MatrixXd eval, int occ, double enuc);
int mp2 (BasisSet obs, MatrixXd evec, MatrixXd eval, int occ, double enuc);
int main() {
    
    	libint2::initialize();  // safe to use libint now	
	int i,j,k;
//path to xyz stored in string
	string h2 = "/Users/kg/HF/hf/h2.xyz";
//reads file at path and stores in variable
	ifstream ip (h2);
//cartesian coords read from file and stored in vector
	vector<Atom> atoms = read_dotxyz(ip);
//no viable conversion. how do i figure out type? auto doesnt work either
//	pair<double, double> aa = make_point_charges(atoms);
//init basis set obs and assoc with sto-3g basis predet in libint
	BasisSet obs("STO-3G",atoms);
//print basis set
	cout << "number basis functions" <<obs.nbf() << endl;
	cout << "Basis Set" << endl;
	std::copy(begin(obs),end(obs),std::ostream_iterator<Shell>(std::cout,"\n"));
	
	MatrixXd density = MatrixXd::Zero(obs.nbf(),obs.nbf());
	
//	cout << "Overlap" << endl;
	MatrixXd S;
	S = overlap (obs);
//	cout << S << endl; 

//	cout << "Kinetic" << endl;
	MatrixXd Kin;
	Kin = kinetic (obs);
//	cout << Kin << endl;
 
//	cout << "Nuclear Attraction" << endl;
	MatrixXd Nuc;
	Nuc = nuclearAtt (obs, atoms);
//	cout << Nuc << endl;
	
	MatrixXd coreH;
	coreH= Kin + Nuc;
//	cout << "Core Hamiltonian" << endl;
//	cout << coreH << endl;
	
	int occ = 1; 
	cout << "HF" << endl;
	MatrixXd evec;
	MatrixXd eval;
	double enuc = hf(obs, S, coreH, 10, atoms, occ, evec, eval);
//	cout << "Eigenvectors \n" << evec << endl;
//	cout << "Eigenvals \n" << eval << endl;
	cout << "MP2 Starts Here" << endl;
	noddy(obs,evec,eval,occ, enuc);
	cout << "The Real MP2 please stand up" << endl;
	mp2(obs,evec,eval,occ, enuc);
	libint2::finalize();  // do not use libint after this 
	return 0;
}

/*
	SMART ALGORITHM FOR MP2
*/
int mp2 (BasisSet obs, MatrixXd evec, MatrixXd eval, int occ, double enuc) {
/*
________ 	________	________
_cip_0_|    *	_p_|_q_|  = 	
_0_|_0_|	_r_|_s_|

*/
	vector<double> intg;
	intg = tensor(obs);
	vector<double> iqrs, iqjs, iajs, iajb;
//	printVec(intg);
//	cout << "nbf" << obs.nbf() << endl;
	int cnt=0;
	double sum = 0;	
	for (int i = 0, cnt=0; i < occ; i++) {
		sum = 0;
		for (int p = 0; p < obs.nbf(); p++) {
			for (int q = 0; q < obs.nbf(); q++) { 
				for (int r = 0; r < obs.nbf(); r++) {
					for (int s = 0; s < obs.nbf(); s++) {
						sum += evec(p,i) * intg[p+q+r+s];
					}
				}
			}	
		} iqrs.push_back(sum); 
	cnt++;
	}
	printVec(iqrs);
/*	sum = 0;
//	cout << "cnt " << cnt << endl;
	for (int j = 0; j <= cnt; j++) {
		sum = 0;
		for (int p = 0; p < obs.nbf(); p++) {
			for (int q = 0; q < obs.nbf(); q++) { 
				for (int r = 0; r < obs.nbf(); r++) {
					for (int s = 0; s < obs.nbf(); s++) {
						sum += evec(q,j) * iqrs[p+q+r+s];
					}
				}
			}	
		} iqjs.push_back(sum); 
	}
	sum = 0;
	for (int k = 0; k <= cnt; k++) {
		sum = 0;
		for (int p = 0; p < obs.nbf(); p++) {
			for (int q = 0; q < obs.nbf(); q++) { 
				for (int r = 0; r < obs.nbf(); r++) {
					for (int s = 0; s < obs.nbf(); s++) {
						sum += evec(r,k) * iqjs[p+q+r+s];
					}
				}
			}	
		} iajs.push_back(sum);  
	}
	sum = 0;
	for (int l = 0; l <= cnt; l++) {
		sum = 0;
		for (int p = 0; p < obs.nbf(); p++) {
			for (int q = 0; q < obs.nbf(); q++) { 
				for (int r = 0; r < obs.nbf(); r++) {
					for (int s = 0; s < obs.nbf(); s++) {
						sum += evec(s,l) * iajs[p+q+r+s];
					}
				}
			}	
		} iajb.push_back(sum); 
	}

*/

//N = # electrons
//range of ijk 0 - N/2
//range abc N/2 - N
	double emp2 = 0;
	for (int i = 0; i < occ; i++){
		for (int a = occ; a < obs.nbf(); a++) {
			for (int j = 0; j < occ; j++) { 
				for (int b = occ; b < obs.nbf(); b++){
					emp2 += (iajb[i+a+j+b]*(2*iajb[i+a+j+b] - iajb[i+b+j+a]))/(eval(i) + eval(j) - eval(a) - eval(b));
				}
			}
		}
	}
//	printVec(pqrs);
	cout << "\t Smart EMP2 " << emp2 << endl;
//	cout << "\t Total E = " << emp2+enuc <<endl;
	return 0;
}

/*
	NODDY ALGORITHM FOR MP2
*/

int noddy (BasisSet obs, MatrixXd evec, MatrixXd eval, int occ, double enuc) {
	vector<double> intg;
	intg = tensor(obs);
	vector<double> pqrs;
//	printVec(intg);
//	cout << "nbf" << obs.nbf() << endl;	
	for (int i = 0; i < obs.nbf(); i++) {
		for (int j = 0; j <= i; j++) {
			for (int k = 0; k <= i; k++) { 
				for (int l = 0; l <= (i==k ? j : k); l++) {
					double sum = 0;
					for (int p = 0; p < obs.nbf(); p++) {
						for (int q = 0; q < obs.nbf(); q++) { 
							for (int r = 0; r < obs.nbf(); r++) {
								for (int s = 0; s < obs.nbf(); s++) {
									sum += evec(p,i) * evec(q,j) * evec(r,k) * evec(s,l) * intg[p+q+r+s];
								}
							}
						}
					} pqrs.push_back(sum);
				}
			}
		}
	}
//N = # electrons
//range of ijk 0 - N/2
//range abc N/2 - N
	double emp2 = 0;
	for (int i = 0; i < occ; i++){
		for (int a = occ; a < obs.nbf(); a++) {
			for (int j = 0; j < occ; j++) { 
				for (int b = occ; b < obs.nbf(); b++){
					emp2 += pqrs[i+a+j+b]*(2*pqrs[i+a+j+b] - pqrs[i+b+j+a])/(eval(i) + eval(j) - eval(a) - eval(b));
//					cout << "Inside loop" << pqrs[i+a+j+b] << endl;
				}
			}
		}
	}
//	cout << "size of " << pqrs.size() << endl;
//	printVec(pqrs);
//	cout << pqrs[0+1+0+1] << endl;
	cout << "\t Noddy EMP2 " << emp2 << endl;
//	cout << "\t Total E = " << emp2+enuc <<endl;
	return 0;

}


/*---------------------*/
/*
	OVERLAP FUNCTION FOR HARTREE-FOCK
*/
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
/*
	KINETIC ENERGY FUNCTION FOR HARTREE-FOCK
*/
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
/*
	NUCLEAR ATTRACTION ENERGY FUNCTION FOR HARTREE-FOCK
*/

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
/*
	FUNCTION TO CALCULATE COULOMB INTEGRALS
*/

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
    for(auto f1=0, cnt=0; f1!=n1; ++f1)
	for(auto f2=0; f2!=n2; ++f2)
	     for(auto f3=0; f3!=n3; ++f3)
      		for(auto f4=0; f4!=n4; ++f4, ++cnt) {
//	cout << cnt << endl;
//      cout << "  " << shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4]  << endl;
	S(bf1+f1,bf2+f2) += shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4] * density(bf3+f3,bf4+f4);
	}	
  }
}
}
}	
return S;
}

/*--------------------------*/
/*
	FUNCTION TO CALCULATE EXCHANGE INTEGRALS
*/
MatrixXd exchange (BasisSet obs, MatrixXd density) {
//K integrals
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
    auto bf2 = shell2bf[s3];  // first basis function in second shell
    auto n2 = obs[s2].size(); // number of basis functions in second shell
    auto bf3 = shell2bf[s2];  // first basis function in first shell
    auto n3 = obs[s3].size(); // number of basis functions in first shell
    auto bf4 = shell2bf[s4];  // first basis function in first shell
    auto n4 = obs[s4].size(); // number of basis functions in first shell
 
    auto shellset = buf_vec[0];
    for(auto f1=0; f1!=n1; ++f1)
	for(auto f3=0; f3!=n3; ++f3)
	     for(auto f2=0; f2!=n2; ++f2)
      		for(auto f4=0; f4!=n4; ++f4) {
//      cout << "  " << shellset[f1*n4*n3*n2 + f4*n4*n3 + f3*n4 + f2]  << endl;
	S(bf1+f1,bf2+f2) += shellset[f1*n3*n2*n4 + f3*n2*n4 + f2*n4 + f4] * density(bf3+f3,bf4+f4);
	}	
  }
}
}
}	
return S;
}

/*-----------------------------------*/
/*
	HARTREE FOCK FUNCTION
*/
double hf(BasisSet obs, MatrixXd overlap, MatrixXd Core, int itr, vector<Atom> atoms, int occ, MatrixXd &evec, MatrixXd &eval) {
//initial guess
	Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> es(Core,overlap);
//eigenvalues and corresponding eigenvectors already sorted
//get eigenvectors of lowest eigenvalue in matrix
	MatrixXd c = es.eigenvectors().col(0);
//	cout << c << endl;
	MatrixXd ct = c.transpose();
//	cout << ct << endl;	
//form density matrix newP
	MatrixXd newP = c*ct;
//	cout << "Initial Density Matrix" << endl;
//	cout << newP << endl;
//initial energy value
	double e1 = 0;
	for(int i = 0; i < obs.nbf(); i++){
		for(int j = 0; j < obs.nbf(); j++) {
			e1 += Core(i,j) * newP(i,j);	
		}
	}
	cout << "First Energy = " << e1 << endl;
	int cnt = 0;
	MatrixXd density = MatrixXd::Zero(obs.nbf(),obs.nbf());
	MatrixXd Fock;
	double enuc, elec; 
	while((newP-density).norm() >= (10^-6) && cnt < 1) {
		density = newP;
		MatrixXd coul = MatrixXd::Zero(obs.nbf(),obs.nbf());
		MatrixXd exch = MatrixXd::Zero(obs.nbf(),obs.nbf());
		coul = coulomb(obs, density);
		exch = exchange(obs, density);
//		cout << "Coulomb \n" << coul << endl;
//		cout << "Exchange \n" << exch << endl;
		Fock = Core + 2*coul - exch;
//		cout << "Fock \n" << Fock << endl;
		Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> es1(Fock,overlap);
		MatrixXd c1 = es1.eigenvectors().col(0);
//		cout << c1 << endl;
		MatrixXd ct1 = c1.transpose();
//		cout << ct1 << endl;	
//form density matrix newP
		MatrixXd newP = c1*ct1;
//		cout << "Density Matrix" << endl;
//		cout << newP << endl;
//iterated energy value
		double e1 = 0;
		for(int i = 0; i < obs.nbf(); i++){
			for(int j = 0; j < obs.nbf(); j++) {
				e1 += (Core(i,j) + Fock(i,j)) * newP(i,j);	
			}
		}	
		elec = e1;
		cout << "Electronic Energy = " << e1 << endl;
//Nuclear Energy
   // compute the nuclear repulsion energy
    enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij * xij + yij * yij + zij * zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << enuc
         << endl;		
		cnt++;
		
	}
	cout << "Total Energy = " << elec+enuc << endl;
//	vector<double> atoms = make_point_charges(atoms);
//	cout << atoms << endl
//for mp2, get all eigenvectors in matrix
//	MatrixXd evec;
	evec = es.eigenvectors();
//	cout << "evec \n" << evec << endl;
	eval = es.eigenvalues();
//	Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> es1(Fock,overlap);
	return enuc;
}

/*-----------------------------------*/
/*
	FUNCTION TO CALCULATE ALL TWO-ELECTRON INTEGRALS FOR MP2
*/
vector<double> tensor (BasisSet obs) {
	vector<double> integral;
		Engine engine(Operator::coulomb,
                	obs.max_nprim(),
                	obs.max_l()
              	 );
		auto shell2bf = obs.shell2bf();
	const auto& buf_vec = engine.results(); // will point to computed shell sets
                                           // const auto& is very important!
//	MatrixXd S = MatrixXd::Zero(obs.nbf(),obs.nbf());
/*	Tensor<double,4> meh(4,4,4,4);
	meh.setZero();
	cout << meh << endl;
*/
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
    for(auto f1=0, cnt=0; f1!=n1; ++f1)
	for(auto f2=0; f2!=n2; ++f2)
	     for(auto f3=0; f3!=n3; ++f3)
      		for(auto f4=0; f4!=n4; ++f4, ++cnt) {
//	cout << cnt << endl;
//      cout << "  " << shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4]  << endl;
//	meh(bf1+f1,bf2+f2,bf3+f3,bf4+f4) = shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4];
	integral.push_back(shellset[f1*n4*n3*n2 + f2*n4*n3 + f3*n4 + f4]);
	}	
  }
}
}
}	
//printVec(integral);
//cout << "1st " << integral[0] << endl;
return integral;
}

/*-----------------------------------*/
/*
	PRINT FUNCTIONS
*/
void printVec (vector <double> vec) {
	for (int i = 0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}

void printInt (vector <int> vec) {
	for (int i = 0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}

}

void printMat(vector<vector<double>> mat) {
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			cout << mat[i][j] << endl;
		}
	} 
}

