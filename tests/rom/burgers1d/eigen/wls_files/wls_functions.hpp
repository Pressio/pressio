/// Outdated
template <typename Mat, typename MatOrVec, typename MatMat>
void local_mat_update(Mat & A, MatOrVec & B, MatMat & C, int scol, int srow, int colSize, int rowSize){
  (*C.data()).block(scol,srow,colSize,rowSize) += (*A.data()).transpose() * (*B.data());
}

template <typename Mat, typename Vec>
void local_vec_update(Mat & A, Vec & b, Vec  & c, int scol,  int colSize){
  auto gradientView = ::pressio::containers::span(c,scol,colSize);
  auto tmp = pressio::containers::ops::dot(A, b);
  for (int k =0; k< colSize; k++){
    gradientView[k] += tmp[k];}
}


template <typename MatMat>
void block_sym_update(MatMat & C,int scol, int srow, int colSize, int rowSize){
  if (scol != srow){
      (*C.data()).block(srow,scol,rowSize,colSize) =  (*C.data()).block(scol,srow,colSize,rowSize).transpose();
  }
}
