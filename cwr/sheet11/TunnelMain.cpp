#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>

using namespace std;

// declaration of the triangular_solver
// look after main() for definition

template <class T>
void triangular_solve(const std::vector<T>& diag,
		      const std::vector<T>& lower,
		      const std::vector<T>& upper,
		      const std::vector<T>& rhs,
		      std::vector<T>& solution);

// declaration of the probability
// look after main() for definition
template <class V>
double probability(const V& psi, int nl, int nr, double dx);

// this is specialication function for the simplified output
// of std::valarray's without the need of writing all the
// elements oneself:
template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) 
{
  int len = v.size();
  
  for(int i = 0; i < len - 1; ++i)
    o << v[i] << " ";
  
  if(len > 0)
    o << v[len-1];
  
  return o;
}

int main(int argc,char **argv) 
{
  double hbar = 1.;
  double mass = 1/2.;
  std::complex<double> complex_I = std::complex<double> (0,1);

  // spatial extent
  double xl = -256., xr = 0.; // width = 256
  int ndx = 1024;              //number of intervalls
  int nx = ndx + 1;           // number of points
  
  double dx = (xr-xl)/ndx;
  std::vector<double> x(nx);
  
  for(int i = 0; i < nx; ++i)
    x[i] = xl + i * dx;
  
  //  tunnel barriere: R_nucleus = 1, V_0 = 10
  double R_n = 50.0, V0 = -1.0;
  std::cerr << "R_nucleus, V_0  = " << R_n << ", " << V0 << std::endl;
  int ibl = int(((xr-xl)-R_n)/dx+0.5);
  int ibr = 0;
  
  double alpha;
  std::cerr << "coefficient alpha = " << std::flush;
  std::cin >> alpha;
  
  // inintial wave package
  double x0 = xl + (xr -xl) * 0.25;
  double sigma0 = (xr-xl) * 0.05;
  int nk = 32; // Wavenumber
  double zk = 2 * M_PI * nk/(xr-xl);
	
  double E0 = hbar * hbar * zk * zk/(2. * mass);
  std::cerr << "E0 = " << E0 << std::endl;
  
  // timestep and ending time
  double dt = 1.0;
  double tfinal = 250.;
  
   // wave functions
  std::vector<std::complex<double> > psi_now(nx), psi_before(nx);
  
  // the diagonal and the next-to-diagonal elements of the Hamiltonian matrix:
  std::vector<std::complex<double> > dH(nx), aH(ndx), cH(ndx);
  
  // matrices in the equations
  std::vector<std::complex<double> > dA(nx), dB(nx), aA(ndx), aB(ndx), cA(ndx), cB(ndx);
  
  double acoeff = hbar/(4. * mass * dx * dx);
  std::complex<double> cacoeff = dt * complex_I * acoeff;
  
  // initialize of the matrices
  for(int i = 0; i < nx; ++i) 
    {      
  
      double V = 0;
      if(i < ibl)
	V = - alpha/x[i];
      else
	V = V0;

      double bcoeff = V/(hbar*2.0);
      std::complex<double> cbcoeff = dt * complex_I * bcoeff;
      
      dH[i] = (2 * acoeff + bcoeff)*2*hbar;
      dA[i] = 1. + 2. * cacoeff + cbcoeff;
      dB[i] = 1. - (2. * cacoeff + cbcoeff);
    }
	

  for(int i = 0; i < ndx; ++i) 
    {
      aH[i] = cH[i] = -acoeff*2*hbar;
      aA[i] = cA[i] = -cacoeff;
      aB[i] = cB[i] = cacoeff;
    }
  

  //inintial wave function;
  for(int i = 0; i < nx; ++i) 
    {
      psi_now[i] = 0.;
      psi_before[i] = std::exp(complex_I * zk * x[i]) * std::exp( -(x[i] - x0) *(x[i] - x0)/(2 * sigma0 * sigma0) );
    }
  
  // normalize the wave function
  double pre_norm = sqrt(probability(psi_before,0,nx-1,dx));
  
  for(int i = 0; i < nx; ++i)
    psi_before[i] /= pre_norm;
  
  std::cerr << "pre_norm = " << pre_norm << std::endl;
  std::cerr << "Norm = " << probability(psi_before,0,nx-1,dx) << std::endl;
  
  // calculate the avergage exact energy
  {
    std::vector<std::complex<double> > psi_tmp(psi_before.size(), 0.);
    
    for(int i = 0; i < nx; ++i)
      psi_tmp[i] += dH[i] * psi_before[i];
    
    for(int i = 0; i < ndx; ++i) 
      {
	psi_tmp[i] += cH[i] * psi_before[i+1];
	psi_tmp[i+1] += aH[i] * psi_before[i];
      }
        
    // energy = quantum mechanical expectation value of the Hamiltonian:
    std::complex<double> erg = 0;
    for(int i = 0; i < nx; ++i)
      erg += std::conj(psi_before[i]) * psi_tmp[i] * dx;
    
    std::cerr << "<H> = " << std::real(erg) << std::endl;
  }
  
  std::cout << 0. << " "
            << probability(psi_before,0,ibl-1,dx) << " " // probability direct before the nucleus
            << probability(psi_before,ibl+1,nx-1,dx) << std::endl; // probability in the nucleus
  
	std::ofstream ofs("out/out.dat");



    //std::ostringstream oss;
    
    //oss << "out/psi_alpha_" << alpha << "_t_" << 0 << ".dat";
    
    //std::ofstream ofs(oss.str().c_str());
    
    for(int i = 0; i < nx; ++i)
      ofs << x[i] << " " << std::abs(psi_before[i]) * std::abs(psi_before[i]) << std::endl;
    
    //ofs.close();
 
  


  
  for(double time = 0; time < tfinal - dt/2.; time += dt) 
    {
      
      // Crank-Nicolson step Eqn
      // implement the Crank-Nicolson algorithm
      std::vector<std::complex<double> > rhs, solution;

	rhs.push_back(dB.at(0)*psi_before.at(0) + cB.at(0)*psi_before.at(1));
	for (size_t i=1; i<dH.size()-1; i++)
	{
	   rhs.push_back(aB.at(i)*psi_before.at(i-1) + dB.at(i)*psi_before.at(i) + cB.at(i)*psi_before.at(i+1));
	}
	rhs.push_back(dB.at(dH.size()-1)*psi_before.at(dH.size()-1) + cB.at(dH.size()-2)*psi_before.at(dH.size()-2));

	triangular_solve(dA, aA, cA, rhs, solution);

	for (size_t i=0; i<solution.size(); i++)
	{
	   psi_before.at(i) = psi_now.at(i);
	   psi_now.at(i) = solution.at(i);	
	}
      			

      // output of the resulting probabilities at this time step:
      std::cout << time+dt << " "
                << probability(psi_now,0,ibl-1,dx) << " " // probability direct before the nucleus
                << probability(psi_now,ibl+1,nx-1,dx) << std::endl; // probability in the nucleus
      
      if( std::abs(int(time+dt+0.1) - (time+dt)) < 1e-5 ) 
	{
	  
	  //std::ostringstream oss;
	  //oss << "out/psi_alpha_" << alpha << "_t_" << time + dt << ".dat";
	  
	  //std::ofstream ofs(oss.str().c_str());
	  
	  for(int i = 0; i < nx; ++i)
	    ofs << x[i] << " " << std::abs(psi_now[i]) * std::abs(psi_now[i]) << std::endl;
	  
	  //ofs.close();
	}
      
      psi_before = psi_now;
    } // end of time evolution loop

} // end of main(...)

//
// calculate the probability to find a particle in
// an interval
//
template <class V>
double probability(const V& psi, int nl, int nr, double dx) 
{
  double retval = 0;
  
  for(int i=nl;i<=nr;++i)
    retval += std::real( psi[i] * std::conj(psi[i]) );
  
  return retval * dx;
}

//
// definition of triangular_solve
//
template <class T>
void triangular_solve(const std::vector<T>& diag,
		      const std::vector<T>& lower,
		      const std::vector<T>& upper,
		      const std::vector<T>& rhs,
		      std::vector<T>& solution) 
{
  
  std::vector<T> new_diag = diag;
  std::vector<T> new_rhs = rhs;
  
  // forward elimination
  for(int i = 1; i < diag.size(); ++i) 
    {
      T pivot = lower[i-1]/new_diag[i-1];
      new_diag[i] -= pivot * upper[i-1];
      new_rhs[i] -= pivot * new_rhs[i-1];
    }
  
  solution.resize(diag.size());
  
  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];
  
  // back substitution
  for(int i = diag.size() - 2; i >= 0; --i) 
    {
      solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
    }
}
