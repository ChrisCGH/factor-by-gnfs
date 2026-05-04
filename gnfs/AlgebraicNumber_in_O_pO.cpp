#include "Matrix.h"
#include "Quotient.h"
#include "VeryLong.h"
#include "VeryLongModular.h"
#include "NumberField.h"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"

template <> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::W_mult_(1,1);
template <> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<long int, VeryLong, LongModular>::W_mult_(1,1);
