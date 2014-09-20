#include "Matrix.h"
#include "Quotient.h"
#include "VeryLong.h"
#include "VeryLongModular.h"
#include "NumberField.h"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"

template <> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::W_mult_(1,1);
template <> Matrix<VeryLongModular> AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::M_(1,1);
template <> VeryLong AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::p_(0L);
template <> VeryLongModular AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::w01_(0L);
template <> VeryLongModular AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::w11_(0L);
template <> bool AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::optimisation_ok_(false);
template <> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::W_mult_(1,1);
template <> Matrix<LongModular> AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::M_(1,1);
template <> long int AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::p_(0L);
template <> LongModular AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::w01_(0L);
template <> LongModular AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::w11_(0L);
template <> bool AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::optimisation_ok_(false);
