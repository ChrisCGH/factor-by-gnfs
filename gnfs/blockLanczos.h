#ifndef BLOCK_LANCZOS_H
#define BLOCK_LANCZOS_H
#include <memory>
#include "SparseMatrix.h"

typedef BitMatrix64 BITMATRIX;
typedef BitOperations64 BITOPERATIONS;
typedef SparseMatrix3 SPARSEMATRIX;
// class to encapsulate the block Lanczos algorithm,
// and add extra features like checkpointing

class BlockLanczos
{
public:
    BlockLanczos(const std::string& matrix_file, const std::string& checkpoint_file = "", int checkpoint_interval = 0, int validation_interval = 0, bool split = false);
    ~BlockLanczos();

    void kernel(BITMATRIX& kerL, BITMATRIX& kerR);

private:
    bool check_A_invertible(const BITMATRIX& Si, const BITMATRIX& VAVi) const;
    bool check_A_orthogonal(const BITMATRIX& Si, const BITMATRIX& Sim1, const BITMATRIX& Vi, const BITMATRIX& Vim1) const;
    void checkpoint();
    void readMatrix(const std::string& matrix_file);
    void readCheckpoint(const std::string& checkpoint_file);
    int checkpoint_interval_;
    int validation_interval_;
// Static data that remains unchanged in loop :
    bool split_;
    size_t n_;
    int N_;
    std::unique_ptr<SPARSEMATRIX> B_;
    std::unique_ptr<BITMATRIX> V0_;
    std::unique_ptr<BITMATRIX> Y_;
// Data that is updated inside loop :
    std::unique_ptr<BITMATRIX> X_;
// Data that is updated at end of loop :
    std::unique_ptr<BITMATRIX> Vim2_;
    std::unique_ptr<BITMATRIX> Vim1_;
    std::unique_ptr<BITMATRIX> Vi_;
    std::unique_ptr<BITMATRIX> Sim1_;
    std::unique_ptr<BITMATRIX> Winvim2_;
    std::unique_ptr<BITMATRIX> Winvim1_;
    std::unique_ptr<BITMATRIX> VAVim1_;
    std::unique_ptr<BITMATRIX> VA2Vim1_;
    int iteration_;
};
#endif
