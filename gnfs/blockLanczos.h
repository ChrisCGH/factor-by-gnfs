#ifndef BLOCK_LANCZOS_H
#define BLOCK_LANCZOS_H
#include "SparseMatrix.h"

typedef BitMatrix64 BITMATRIX;
typedef BitOperations64 BITOPERATIONS;
//typedef BitMatrix BITMATRIX;
//typedef BitOperations BITOPERATIONS;
typedef SparseMatrix3 SPARSEMATRIX;
// class to encapsulate the block Lanczos algorithm,
// and add extra features like checkpointing

class BlockLanczos
{
   public:
      BlockLanczos(const std::string& matrix_file, const std::string& checkpoint_file = "", int checkpoint_interval = 0, bool split = false);
      ~BlockLanczos();

      void kernel(BITMATRIX& kerL, BITMATRIX& kerR);

   private:
      bool check_A_invertible(const BITMATRIX& Si, const BITMATRIX& VAVi);
      void checkpoint();
      void readMatrix(const std::string& matrix_file);
      void readCheckpoint(const std::string& checkpoint_file);
      int checkpoint_interval_;
// Static data that remains unchanged in loop :
      bool split_;
      size_t n_;
      int N_;
      SPARSEMATRIX* B_;
      BITMATRIX* V0_;
      BITMATRIX* Y_;
// Data that is updated inside loop :
      BITMATRIX* X_;
// Data that is updated at end of loop :
      BITMATRIX* Vim2_;
      BITMATRIX* Vim1_;
      BITMATRIX* Vi_;
      BITMATRIX* Sim1_;
      BITMATRIX* Winvim2_;
      BITMATRIX* Winvim1_;
      BITMATRIX* VAVim1_;
      BITMATRIX* VA2Vim1_;
      int iteration_;
};
#endif
