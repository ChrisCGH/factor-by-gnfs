#include "blockLanczos.h"
#include <iostream>
#include "MemoryMappedFile.h"
#include <time.h>
namespace
{
char* time_str()
{
    time_t the_time = time(0);
    static char the_time_str[132];
    strftime(the_time_str, 132, "%Y%m%d%H%M%S", localtime(&the_time));
    return the_time_str;
}

void die(const std::string& msg)
{
    std::cerr << msg << std::endl;
    std::exit(1);
}
}

BlockLanczos::BlockLanczos(const std::string& matrix_file, const std::string& checkpoint_file, int checkpoint_interval, int validation_interval, bool split) : split_(split)
{
    std::cerr << "Starting ..." << std::endl;
    readMatrix(matrix_file);

    if (checkpoint_file == "")
    {
        V0_ = new BITMATRIX(n_, N_);
        X_ = new BITMATRIX(n_, N_);
        Vim2_ = new BITMATRIX(n_, N_);
        Vim1_ = new BITMATRIX(n_, N_);
        Vi_ = new BITMATRIX(n_, N_);
        Sim1_ = new BITMATRIX(N_);
        Winvim2_ = new BITMATRIX(N_, N_);
        Winvim1_ = new BITMATRIX(N_, N_);
        VAVim1_ = new BITMATRIX(N_, N_);
        VA2Vim1_ = new BITMATRIX(N_, N_);

        // choose a random n x N matrix
        Y_ = new BITMATRIX(n_, N_);
        randomise(*Y_);
        sym_multiply(*B_, *Y_, *Vi_);
        *V0_ = *Vi_;
        if (*V0_ != *Vi_) std::cerr << "Problem 0" << std::endl;
        iteration_ = 0;
    }
    else
    {
        readCheckpoint(checkpoint_file);
    }
    checkpoint_interval_ = checkpoint_interval;
    const long int default_checkpoint_interval = 1000L;
    if (!checkpoint_interval_) checkpoint_interval_ = default_checkpoint_interval;
    validation_interval_ = validation_interval;
    if (!validation_interval_) validation_interval_ = 100L;
}

BlockLanczos::~BlockLanczos()
{
    delete B_;
    delete V0_;
    delete X_;
    delete Vim2_;
    delete Vim1_;
    delete Vi_;
    delete Sim1_;
    delete Winvim2_;
    delete Winvim1_;
    delete VAVim1_;
    delete VA2Vim1_;
    delete Y_;
}

void BlockLanczos::readMatrix(const std::string& matrix_file)
{
    std::cerr << "reading sparse matrix ..." << std::endl;
    B_ = new SPARSEMATRIX(matrix_file, split_);
    std::cerr << "rows = " << static_cast<unsigned int>(B_->rows()) << std::endl;
    std::cerr << "cols = " << static_cast<unsigned int>(B_->cols()) << std::endl;
    n_ = B_->cols();
    N_ = BITOPERATIONS::BITS_IN_WORD;
}

void BlockLanczos::kernel(BITMATRIX& kerL, BITMATRIX& kerR)
{
    std::cerr << "blockLanczos: started" << std::endl;
    size_t n = n_;
    int N = N_;
    SPARSEMATRIX& B = *B_;

    BITMATRIX& V0 = *V0_;
    BITMATRIX& Vi = *Vi_;
    BITMATRIX& Vim1 = *Vim1_;
    BITMATRIX& Vim2 = *Vim2_;

    BITMATRIX Si(N);
    BITMATRIX& Sim1 = *Sim1_;

    BITMATRIX Winvi(N, N);
    BITMATRIX& Winvim1 = *Winvim1_;
    BITMATRIX& Winvim2 = *Winvim2_;

    BITMATRIX VAVi(N);
    BITMATRIX& VAVim1 = *VAVim1_;
    BITMATRIX VA2Vi(N);
    BITMATRIX& VA2Vim1 = *VA2Vim1_;

    // reference identify matrix
    const BITMATRIX I_N(N);

    BITMATRIX& Y = *Y_;
    BITMATRIX& X = *X_;

    int done = 0;
    while (!done)
    {
        std::cerr << "blockLanczos: iteration <" << iteration_ << ">" << std::endl;
        BITMATRIX AV;
        sym_multiply(B, Vi, AV); // B^t B Vi
        innerProduct(Vi, AV, VAVi);   // Vi^t B^t B Vi
        if (VAVi.isZero()) done = 1;
        else
        {
            BITMATRIX tmp;
            BITMATRIX tmp1;

            chooseS(VAVi, Sim1, Si, Winvi);

            // calculate next term in X, before doing next Block-Lanczos iteration:
            //             -1 t
            // X <- X + V W  V V
            //           i i  i 0
            innerProduct(Vi, V0, tmp);
            multiply(Winvi, tmp, tmp1);
            multiply(Vi, tmp1, tmp);
            X += tmp;

            // Calculate D, E and F
            BITMATRIX SS_t;
            sym_multiply(Si, SS_t);

            //          -1  t 2     t   t
            // D = I - W  (V A V S S + V AV )
            //          i   i   i i i   i  i
            BITMATRIX D(N);
            innerProduct(AV, AV, VA2Vi);
            multiply(VA2Vi, SS_t, tmp);
            tmp += VAVi;
            multiply(Winvi, tmp, tmp1);
            D -= tmp1;

            //        t -1  t
            // E = S S W   V AV
            //      i i i-1 i  i
            BITMATRIX E;
            multiply(Winvim1, VAVi, tmp);
            multiply(tmp, SS_t, E);

            //              t
            // newV = AV S S  + V D + V   E
            //          i i i    i     i-1
            BITMATRIX newV(n, N);
            multiply(AV, SS_t, newV);
            multiply(Vi, D, tmp);
            newV += tmp;
            multiply(Vim1, E, tmp);
            newV += tmp;

            // If Sim1 is identity matrix, then we don't need to
            // calculate F
            if (Sim1 != I_N)
            {
                //      -1       t         -1    t   2         t      t           t
                // F = W   (I - V   A V   W   )(V   A V   S   S    + V   AV   )S S
                //      i-2      i-1   i-1 i-1   i-1   i-1 i-1 i-1    i-1  i-1  i i
                BITMATRIX F;
                multiply(VAVim1, Winvim1, tmp);
                BITMATRIX tmp2(N);
                tmp2 -= tmp;
                multiply(Winvim2, tmp2, tmp);

                BITMATRIX SSm1_t;
                sym_multiply(Sim1, SSm1_t);
                multiply(VA2Vim1, SSm1_t, tmp1);
                tmp1 += VAVim1;
                multiply(tmp, tmp1, tmp2);
                multiply(tmp2, SS_t, F);
                multiply(Vim2, F, tmp);
                //
                // newV <- newV + V   F
                //                 i-2
                newV += tmp;
            }

            // We could check here that:
            //
            //       t                          t t
            // (a)  W A W  is invertible, i.e. S V A V S is invertible
            //       i   i                      i i   i i
            //
            //       t
            // (b)  W A W = 0 for i != j
            //       i   j
            //
            //
            if (!check_A_invertible(Si, VAVi))
            {
                std::cerr << "          t"                       << std::endl;
                std::cerr << "Problem: W A W  is not invertible" << std::endl;
                std::cerr << "          i   i"                   << std::endl;
            }

            if (!check_A_orthogonal(Si, Sim1, AV, Vim1))
            {
                std::cerr << "          t"                         << std::endl;
                std::cerr << "Problem: W   A W  is not zero" << std::endl;
                std::cerr << "          i-1   i"                   << std::endl;
            }

            Vim2 = Vim1;
            Vim1 = Vi;
            Vi = newV;
            Sim1 = Si;
            Winvim2 = Winvim1;
            Winvim1 = Winvi;
            VAVim1 = VAVi;
            VA2Vim1 = VA2Vi;
        }
        iteration_++;
        if (!done) checkpoint();
    }

    // build an n x 2N matrix Z from X - Y and V
    //                                          m
    // and calculate BZ.
    if (Vi.isZero())
    {
        std::cerr << "blockLancos: V  is zero" << std::endl;
        std::cerr << "              m" << std::endl;
        std::cerr << "blockLanczos: building BZ ... ";
        BITMATRIX Z(X);
        Z += Y;
        // check that AZ = 0;
        BITMATRIX AZ;
        sym_multiply(B, Z, AZ);
        if (!AZ.isZero())
        {
            std::cerr << ": Problem: AZ != 0 : ";
            std::cerr << std::endl << "AZ = " << AZ << std::endl;
        }
        BITMATRIX AX;
        sym_multiply(B, X, AX);
        BITMATRIX AY;
        sym_multiply(B, Y, AY);

        if (AX != AY)
        {
            std::cerr << ": Problem: AX != AY : ";
        }
        BITMATRIX BZ;
        multiply(B, Z, BZ);
        std::cerr << "done" << std::endl;
        if (BZ.isZero())
        {
            std::cerr << "blockLanczos: Z is already the kernel" << std::endl;
        }
        else
        {
            std::cerr << "blockLanczos: BZ = " << std::endl;
        }

        // Find U whose columns span the null space of BZ using Gaussian reduction

        BITMATRIX U;
        std::cerr << "blockLanczos: finding kernel U of BZ ... ";
        ::kernel(BZ, U);
        std::cerr << "done" << std::endl;
        // check
        BITMATRIX ZZ;
        multiply(BZ, U, ZZ);
        if (!ZZ.isZero())
        {
            std::cerr << "Problem in kernel" << std::endl;
        }

        // Now calculate ZU
        BITMATRIX ZU;
        std::cerr << "blockLanczos: calculating ZU ... ";
        multiply(Z, U, ZU);
        std::cerr << "done" << std::endl;

        kerL = ZU;

        if (split_)
        {
            // Now process dense rows
            // Calculate product of dense part of B by kernel
            BITMATRIX B2X;
            B.multiply_dense_part_by_bit_matrix(ZU, B2X);

            // Find kernel of B2X
            BITMATRIX Y;
            ::kernel(B2X, Y);

            // Calculate ZU Y
            multiply(ZU, Y, kerL);
        }
    }
    else
    {
        std::cerr << "blockLanczos: building BZ ... ";
        BITMATRIX ZL(X);
        ZL += Y;

        BITMATRIX ZR(Vi);
        BITMATRIX BZL;
        multiply(B, ZL, BZL);
        BITMATRIX BZR;
        multiply(B, ZR, BZR);
        std::cerr << "done" << std::endl;

        // Find U whose columns span the null space of BZ using Gaussian reduction

        BITMATRIX UL;
        BITMATRIX UR;
        std::cerr << "blockLanczos: finding kernel U of BZ ... ";
        ::kernel(BZL, BZR, UL, UR);
        std::cerr << "done" << std::endl;

        // Now calculate ZU
        BITMATRIX ZUL;
        BITMATRIX ZUR;
        std::cerr << "blockLanczos: calculating ZU ... ";
        multiply(ZL, ZR, UL, UR, ZUL, ZUR);
        std::cerr << "done" << std::endl;

        //std::cerr << "blockLanczos: ZUR = " << ZUR << std::endl;
        // Check
        BITMATRIX BZUL;
        BITMATRIX BZUR;
        multiply(B, ZUL, BZUL);
        multiply(B, ZUR, BZUR);
        if (!BZUL.isZero())
        {
            std::cerr << "Problem in kernel" << std::endl;
        }
        if (!BZUR.isZero())
        {
            std::cerr << "Problem in kernel" << std::endl;
        }
        kerL = ZUL;
        kerR = ZUR;

        if (split_)
        {
            // Now process dense rows
            // Calculate product of dense part of B by kernel
            BITMATRIX B2XL;
            BITMATRIX B2XR;
            B.multiply_dense_part_by_bit_matrix(ZUL, ZUR, B2XL, B2XR);

            // Find kernel of B2X
            BITMATRIX Y;
            BITMATRIX Z;
            ::kernel(B2XL, B2XR, Y, Z);

            multiply(ZUL, ZUR, Y, Z, kerL, kerR);
        }
    }
}

void BlockLanczos::checkpoint()
{
    if (iteration_ % checkpoint_interval_ != 0L) return;
    std::string checkpoint_file("blcp_");
    char* t = time_str();
    checkpoint_file += t;
    checkpoint_file += ".dat";
    std::cerr << "writing checkpoint " << checkpoint_file << " ..." << std::endl;
    /*
    // Static data that remains unchanged in loop :
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
    */
    std::fstream cp(checkpoint_file.c_str(), std::ios::out);
    cp << t << " : Block Lanczos checkpoint file" << std::endl;
    cp << "iteration " << iteration_ << std::endl;
    cp << "V0" << std::endl;
    V0_->write(cp);
    cp << "Y" << std::endl;
    Y_->write(cp);
    cp << "X" << std::endl;
    X_->write(cp);
    cp << "Vim2" << std::endl;
    Vim2_->write(cp);
    cp << "Vim1" << std::endl;
    Vim1_->write(cp);
    cp << "Vi" << std::endl;
    Vi_->write(cp);
    cp << "Sim1" << std::endl;
    Sim1_->write(cp);
    cp << "Winvim2" << std::endl;
    Winvim2_->write(cp);
    cp << "Winvim1" << std::endl;
    Winvim1_->write(cp);
    cp << "VAVim1" << std::endl;
    VAVim1_->write(cp);
    cp << "VA2Vim1" << std::endl;
    VA2Vim1_->write(cp);
    cp << "End of Block Lanczos checkpoint file" << std::endl;
}

void BlockLanczos::readCheckpoint(const std::string& checkpoint_file)
{
    std::cerr << "reading checkpoint " << checkpoint_file << " ..." << std::endl;
    //std::fstream cp(checkpoint_file.c_str(), std::ios::in);
    MemoryMappedFile cp(checkpoint_file.c_str());
    std::string str;
    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    // should be
    // 0123456789012345
    // YYYYMMDDHHMMSS : Block Lanczos checkpoint file"
    if (strcmp(str.c_str() + 14, " : Block Lanczos checkpoint file") != 0) die("corrupted checkpoint file header");

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    // should be
    // 01234567890
    // iteration nnnnnn
    if (strncmp(str.c_str(), "iteration ", 10) != 0) die("corrupted checkpoint file - iteration");
    if (!isdigit(*(str.c_str() + 10))) die("corrupted checkpoint file - iteration count");
    iteration_ = std::atol(str.c_str() + 10);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "V0") != 0) die("corrupted checkpoint file - V0");
    V0_ = new BITMATRIX;
    V0_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Y") != 0) die("corrupted checkpoint file - Y");
    Y_ = new BITMATRIX;
    Y_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "X") != 0) die("corrupted checkpoint file - X");
    X_ = new BITMATRIX;
    X_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Vim2") != 0) die("corrupted checkpoint file - Vim2");
    Vim2_ = new BITMATRIX;
    Vim2_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Vim1") != 0) die("corrupted checkpoint file - Vim1");
    Vim1_ = new BITMATRIX;
    Vim1_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Vi") != 0) die("corrupted checkpoint file - Vi");
    Vi_ = new BITMATRIX;
    Vi_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Sim1") != 0) die("corrupted checkpoint file - Sim1");
    Sim1_ = new BITMATRIX;
    Sim1_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Winvim2") != 0) die("corrupted checkpoint file - Winvim2");
    Winvim2_ = new BITMATRIX;
    Winvim2_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "Winvim1") != 0) die("corrupted checkpoint file - Winvim1");
    Winvim1_ = new BITMATRIX;
    Winvim1_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "VAVim1") != 0) die("corrupted checkpoint file - VAVim1");
    VAVim1_ = new BITMATRIX;
    VAVim1_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "VA2Vim1") != 0) die("corrupted checkpoint file - VA2Vim1");
    VA2Vim1_ = new BITMATRIX;
    VA2Vim1_->read(cp);

    if (!getline(cp, str)) die("unexpected end of checkpoint file");
    if (strcmp(str.c_str(), "End of Block Lanczos checkpoint file") != 0) die("corrupted checkpoint file trailer");
}

bool BlockLanczos::check_A_invertible(const BITMATRIX& Si, const BITMATRIX& VAVi) const
{
    if (iteration_ % validation_interval_ != 0L) return true;
    std::cerr << "               t"                   << std::endl;
    std::cerr << "Checking that W A W  is invertible" << std::endl;
    std::cerr << "               i   i"               << std::endl;
    BITMATRIX VAViSi;
    multiply(VAVi, Si, VAViSi);     // Vi^t A Vi Si
    BITMATRIX WAWi;
    innerProduct(Si, VAViSi, WAWi); // Si^t Vi^t A Vi Si = Wi^t A Wi
    BITMATRIX kerWAWi;
    ::kernel(WAWi, kerWAWi);
    return kerWAWi.isZero();
}

bool BlockLanczos::check_A_orthogonal(const BITMATRIX& Si, const BITMATRIX& Sim1, const BITMATRIX& AVi, const BITMATRIX& Vim1) const
{
    if (iteration_ % validation_interval_ != 0L) return true;
    std::cerr << "               t"               << std::endl;
    std::cerr << "Checking that W   A W  is zero" << std::endl;
    std::cerr << "               i-1   i"         << std::endl;
    BITMATRIX Vim1AVi;
    innerProduct(Vim1, AVi, Vim1AVi); // Vi-1^t A Vi
    BITMATRIX VAViSi;
    multiply(Vim1AVi, Si, VAViSi);     // Vi^t A Vi Si
    BITMATRIX Wim1AWi;
    innerProduct(Sim1, VAViSi, Wim1AWi); // Si^t Vi^t A Vi Si = Wi^t A Wi
    return Wim1AWi.isZero();
}
