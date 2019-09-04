/*
//@HEADER
// ************************************************************************
//
// codes.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the 
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions 
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
// notice, this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its 
// contributors may be used to endorse or promote products derived 
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

//FRizzi: @Chi, please fix all this 


/* assume that the bases is arranged as follows:
[basis_1__timeStep_1, basis_2__timeStep_1, ...,basis_nst__timeStep_1,   basis_1__timeStep_2, basis_2__timeStep_2, ...,basis_nst__timeStep_2, ... ... basis_1__timeStep_Nt, basis_2__timeStep_Nt, ...,basis_nst__timeStep_Nt ]
*/

#ifndef MULTI_SPACE_TIME_BASES_HPP
#define MULTI_SPACE_TIME_BASES_HPP

#include "../rom_ConfigDefs.hpp"

namespace pressio{ namespace rom{ namespace experimental{

template <typename T>
class MultiSpaceTimeBases : public pressio::containers::MultiVector<T>
{

private:
	using this_t = MultiSpaceTimeBases<T>;
	using wrap_t = typename pressio::containers::details::traits<this_t>::wrapped_t;

public:

	MultiSpaceTimeBases(const this_t & other)
		: data_(*other.data()){}

	explicit MultiSpaceTimeBases(const wrap_t & vecobj) : data_(vecobj){}	

	~MultiSpaceTimeBases() = delete;
	
	void GetDimensions(int Ns, int Nt, int nst) const {
		Ns_ = Ns;
		Nt_ = Nt;
		nst_ = nst;
	}


	// calculate uTilde = Phi_st * uhat (strictly for Phi_st only)
	void MultiplyRomVector(const romVec_t & romVec, MultiVector_t & Ust){
		// check size of romVec?
		if (romVec.size()!=nst_){
			std::cout << "CH error: size of romVec is not equal to number of space time bases!" << std::endl;
		}

		// multiply: multiply each column block of multiVector with romVec. Note that all indices here (i,j,ir,ic) are 0-based.
		MultiVector_t tmpMat(Ns_,nst_);		// (#rows, #cols)
		MultiVector_t tmpVec(Ns_,1);
		for (int i = 0; i < Nt_; ++i){
			
			// extract each column block
			for (int ir = 0; ir < Ns_; ++ir){
				int colCount=0;
				for (int ic = i*nst_; ic < (i+1)*nst_; ++i){
					tmpMat[ir][colCount] = data_[ir][ic];
					colCount+=1;
				}
			}

			// multiply: use dot operator already defined in library
			tmpVec = tmpMat * (*romVec.data());

			// assign to column "i" of Ust
			for (int ir = 0; ir < Ns_; ++ir){
				*Ust.data()[ir][i] = tmpVec[ir];
			}
		}
	}	// end MultiplyRomVector


	// calculate JP = Jac*Phi_st (strictly for Phi_st only)
	this_t & JP CalculateJacobianTimesPhiSpaceTime(const MultiVector_t & u_st, const probType & probObj){
		// NOTE: all indices here are 0-based as opposed with 1-based indices in the report

		MultiVector_t u_s(Ns_,1);
		MultiVector_t phiS_ij(Ns_,1);
		MultiVector_t phiS_im1j(Ns_,1);
		MultiVector_t tmpVec(Ns_,1);
		MultiVector_t tmpVec_im1(Ns_,1);
		int numNonZeros = Ns_*9;
		sparseMatrix_t J_ii(Ns_, Ns_, numNonZeros);	// pseudo-code, fix later
		sparseMatrix_t J_iim1 = identityMatrix(Ns_); // pseudo-code, fix later

		for (int i = 0; i < Nt_; ++i){

			// extract u at time step i
			for (int ir = 0; ir < Ns_; ++ir){
				u_s[ir] = *u_st.data()[ir][i];
			}

			// compute spatial jacobian
			J_ii = computeJacobian(u_s, probObj);	// pseudo-code, fix later

			for (int j = 0; j < nst_; ++j){
				int colInd = i*Nt_+j;

				// extract phiS_ij
				for (int ir = 0; ir < Ns_; ++ir){
					phiS_ij[ir] = data_[ir][colInd];
				}

				// multiply: use dot operator already defined in library
				tmpVec = J_ii * phiS_ij;

				if (i==0){
					//------ only need J_i,i and phiS_i,j ------//
				}
				else {
					//------ need J_i,i J_i,i-1 and phiS_i,j phiS_i-1,j  ------//
					// need jacobian J_i,i-1 -- only for Backward Euler

					// need phiS_i-1,j
					int colInd_im1j = (i-1)*Nt_+j;
					for (int ir = 0; ir < Ns_; ++ir){
						phiS_im1j[ir] = data_[ir][colInd_im1j];
					}

					// multiply and sum
					tmpVec_im1 = J_iim1 * phiS_im1j;
					tmpVec += tmpVec_im1;
				}

				// assign to final result
				for (int ir = 0; ir < Ns_; ++ir){
					JP[ir][colInd] = tmpVec[ir];
				}

			} // end for loop j
		} // end for loop i

		return *JP;
	} // end CalculateJacobianTimesPhiSpaceTime


	// calculate RHS = JP^T * r_st (strictly for JP only)
	vec_type & vec MultiplySpaceTimeResidual(const MultiVector_t & r_st){

		MultiVector_t tmpVecL(Ns_,1);
		MultiVector_t tmpVecR(Ns_,1);

		// loop over all nst_ (rows)
		for (int i = 0; i < nst_; ++i){
			
			double sum=0.0;
			for (int j = 0; j < Nt_; ++j){
				int colInd_JP = i + j*nst_;
				int colInd_r  = j;

				// product of two vectors
				double prod=0.0;
				for (int ir = 0; ir < Ns_; ++ir){
					tmpVecL[ir] = data_[ir][colInd_JP];
					tmpVecR[ir] = *r_st.data()[ir][colInd_r];
					prod += tmpVecL[ir] * tmpVecR[ir];
				}

				// sum
				sum += prod;
			} // end for loop j

			// assign to result (vec)
			*vec.data()[i] = sum;
		} // end for loop i

		return *vec;
	} // end MultiplySpaceTimeResidual


	// calculate LHS = JP^T * JP (strictly for JP only)
	mvec_type & mat MultiplyItself(){

		MultiVector_t tmpVecL(Ns_,1);
		MultiVector_t tmpVecR(Ns_,1);		

		// loop over all nst_ (rows)
		for (int i = 0; i < nst_; ++i){
			// loop over all nst_ (cols)
			for (int j = 0; j < nst_; ++j){
				
				double sum=0.0;
				for (int it = 0; it < Nt_; ++it){
					int colInd_left = i + it*nst_;
					int colInd_right= j + it*nst_;

					// product of two vectors
					double prod=0.0;
					for (int ir = 0; ir < Ns_; ++ir){
						tmpVecL[ir] = data_[ir][colInd_left];
						tmpVecR[ir] = data_[ir][colInd_right];
						prod += tmpVecL[ir] * tmpVecR[ir];
					}
					sum += prod;
				} // end for loop Nt

				// assign to result (mat)
				*mat.data()[i][j] = sum;
			}
		}

		return *mat
	} // end MultiplyItself


private:
	wrap_t data_ = {};
	int Ns_;
	int Nt_;
	int nst_;
};

}}}//end namespace pressio::rom::experimental
#endif
