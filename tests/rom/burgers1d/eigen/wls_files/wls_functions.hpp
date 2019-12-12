template <typename Mat, typename MatOrVec, typename MatMat>
void local_mat_update(Mat & A, MatOrVec & B, MatMat & C, int scol, int srow, int colSize, int rowSize){
  (*C.data()).block(scol,srow,colSize,rowSize) += (*A.data()).transpose() * (*B.data());
}

template <typename MatMat>
void block_sym_update(MatMat & C,int scol, int srow, int colSize, int rowSize){
  if (scol != srow){
      (*C.data()).block(srow,scol,rowSize,colSize) =  (*C.data()).block(scol,srow,colSize,rowSize).transpose();
  }
}
