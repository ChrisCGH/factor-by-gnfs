OUTPUTDIR = gbin
GMPLIBDIR = /usr/local/lib
#LIBS = -L/usr/lib -L/usr/local/lib -lgmp -Wl,--stack,200000000
#LIBS = -L/usr/lib -L/usr/local/lib -lgmp -Wl,--verbose
#LIBS = -L/usr/lib -L/usr/local/lib -lgmp -lstdc++ -ldmalloc
COVERAGE_LIBS = -lgcov
ifeq ($(RUN_GCOV),yes)
   LIBS = -L/usr/lib64 -L/usr/local/lib $(COVERAGE_LIBS) -luuid -lgmp -lstdc++
   TEST_LIBS = -L/usr/lib64 -L/usr/local/lib $(COVERAGE_LIBS) -luuid -lgmp -lcppunit -ldl -lstdc++
else
   LIBS = -L/usr/lib64 -L/usr/local/lib -luuid -lgmp -lstdc++
   TEST_LIBS = -L/usr/lib64 -L/usr/local/lib -luuid -lgmp -lcppunit -ldl -lstdc++
endif
#LIBS = -L/usr/lib -L/usr/local/lib -lgmp -lstdc++ -Wl,--stack,209715200,--heap,1073741824
INCLUDES = -I/usr/local/include
#OPT = -O3 -mfpmath=sse
#OPT = -O3
#PROFILE = -pg
COVERAGE = -fprofile-arcs -ftest-coverage
COVERAGE_OPT = 
#PROFILE = -fprofile-arcs -fbranch-probabilities
#PROFILE = 
ARCH = -march=athlon64
USING_CPPUNIT = -DUSING_CPPUNIT 
#ARCH = -march=pentium4 -mfpmath=sse -momit-leaf-frame-pointer 
#ARCH = -march=pentium4 
CCFLAGS = -Wall -O3 -g -fbranch-probabilities 
CCFLAGS = -Wall -O3 -g -fprofile-arcs
CCFLAGS = -Wall -g -fprofile-arcs -ftest-coverage
ifeq ($(RUN_GCOV),yes)
  CCFLAGS = -Wall -g $(COVERAGE_OPT) $(COVERAGE) $(PROFILE) $(ARCH) -fno-operator-names -Wno-non-template-friend -Wno-uninitialized -DUSING_GCC $(USING_CPPUNIT)
  CFLAGS = -Wall -g $(COVERAGE_OPT) $(COVERAGE) $(PROFILE) $(ARCH)
else
  CCFLAGS = -Wall -g $(OPT) $(PROFILE) $(ARCH) -fno-operator-names -Wno-non-template-friend -Wno-uninitialized -DUSING_GCC $(USING_CPPUNIT)
  CFLAGS = -Wall -g $(OPT) $(PROFILE) $(ARCH)
endif
CC = g++

.SUFFIXES = .o .cpp .c
.cpp.o: 
	g++ -c $(CCFLAGS) -Wa,-ahls=$*.lst $<

.c.o:
	gcc -c $(CFLAGS) -Wa,-ahls=$*.lst $<

%.d: %.cpp
	set -e; $(CC) -MM $(CPPFLAGS) $< \
       | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; [ -s $@ ] || rm -f $@

%.d: %.c
	set -e; gcc -MM $(CPPFLAGS) $< \
       | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; [ -s $@ ] || rm -f $@

ALL_TESTS = AlgebraicNumberUnitTest BitOperationsUnitTest LongModularUnitTest VeryLongUnitTest VeryLongModularUnitTest MatrixUnitTest PolynomialUnitTest CombinationsUnitTest HashTableUnitTest QuotientUnitTest PolynomialOptimizerUnitTest LatticeSieverUnitTest MemoryMappedFileUnitTest SparseMatrixUnitTest

all : $(ALL_SRCS) $(OUTPUTDIR) skewed lsieve bl calcroot sieve filter buildMatrix processDependenciespl postProcessing DenseRowTest stest testan gen testhnf inv testmmf convert testfactor processDependencies $(ALL_TESTS)

test : $(ALL_TESTS)
	$(OUTPUTDIR)/AlgebraicNumberUnitTest.exe 
	$(OUTPUTDIR)/BitOperationsUnitTest.exe
	$(OUTPUTDIR)/QuotientUnitTest.exe
	$(OUTPUTDIR)/LongModularUnitTest.exe
	$(OUTPUTDIR)/VeryLongUnitTest.exe
	$(OUTPUTDIR)/VeryLongModularUnitTest.exe
	$(OUTPUTDIR)/MatrixUnitTest.exe
	$(OUTPUTDIR)/SparseMatrixUnitTest.exe
	$(OUTPUTDIR)/HashTableUnitTest.exe
	$(OUTPUTDIR)/CombinationsUnitTest.exe
	$(OUTPUTDIR)/PolynomialUnitTest.exe
	$(OUTPUTDIR)/PolynomialOptimizerUnitTest.exe
	$(OUTPUTDIR)/LatticeSieverUnitTest.exe
	$(OUTPUTDIR)/MemoryMappedFileUnitTest.exe

$(OUTPUTDIR) :
	mkdir $(OUTPUTDIR)

ALL_HEADERS = Combinations.h Config.h ContinuedFraction.h ExceptionalPrimes.h FactorBase.h Ideal.h LatticeSiever.h LongModular.h MPFloat.h Matrix.h NumberField.h AlgebraicNumber.h AlgebraicNumber_in_O_pO.h Polynomial.h QuadraticCharacters.h Quotient.h RootConfig.h SieveConfig.h VeryLong.h VeryLongModular.h SparseMatrix.h crt.h dickman.h blockLanczos.h lip.h lippar.h lll.h mod.h mt19937int.h pow.h pselect.h root.h timings.h Siever.h squfof.h PointerHashTable.h graph.h RelationManager.h discriminant.h gcd.h convert.h qs.h PolynomialOptimizer.h SieveUtils.h
ALL_CPPS = FactorBase.cpp Ideal.cpp LatticeSiever.cpp LongModular.cpp NumberField.cpp AlgebraicNumber.cpp AlgebraicNumber_in_O_pO.cpp Polynomial.cpp SparseMatrix.cpp QuadraticCharacters.cpp VeryLong.cpp VeryLongModular.cpp bl.cpp blockLanczos.cpp calcroot.cpp dickman.cpp discriminant.cpp lll.cpp lsieve.cpp pselect.cpp root.cpp skewed.cpp timings.cpp sieve.cpp Siever.cpp filter.cpp squfof.cpp RelationManager.cpp buildMatrix.cpp stest.cpp testhnf.cpp testan.cpp MatrixUnitTest.cpp convert.cpp convertmain.cpp testMatrix.cpp testfactor.cpp qs.cpp processDependencies.cpp AlgebraicNumberUnitTest.cpp BitOperationsUnitTest.cpp LongModularUnitTest.cpp VeryLongUnitTest.cpp VeryLongModularUnitTest.cpp CombinationsUnitTest.cpp HashTableUnitTest.cpp QuotientUnitTest.cpp LatticeSieverUnitTest.cpp PolynomialOptimizer.cpp SieveUtils.cpp ef.cpp MemoryMappedFileUnitTest.cpp SparseMatrixUnitTest.cpp

ALL_CS = lip.c mt19937int.c 
ALL_SRCS = $(ALL_CPPS) $(ALL_CS) $(ALL_HEADERS)
ALL_OBJS = $(ALL_CPPS:.cpp=.o) $(C_SRCS:.c=.o)

QS_CPP_SRCS = qs.cpp SparseMatrix.cpp timings.cpp VeryLongModular.cpp VeryLong.cpp squfof.cpp blockLanczos.cpp

C_SRCS = mt19937int.c lip.c

include $(QS_CPP_SRCS:.cpp=.d)
include $(C_SRCS:.c=.d)

QS_OBJS = $(QS_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

qs : $(OUTPUTDIR) $(OUTPUTDIR)/qs
	

$(OUTPUTDIR)/qs : $(QS_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/qs.exe $(QS_OBJS) $(LIBS)

cleanqs: 
	rm -rf $(OUTPUTDIR)/qs.exe $(QS_OBJS)

TESTHNF_CPP_SRCS = testhnf.cpp timings.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp squfof.cpp qs.cpp SparseMatrix.cpp

C_SRCS = mt19937int.c lip.c

include $(TESTHNF_CPP_SRCS:.cpp=.d)
include $(C_SRCS:.c=.d)

TESTHNF_OBJS = $(TESTHNF_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

testhnf : $(OUTPUTDIR) $(OUTPUTDIR)/testhnf
	

$(OUTPUTDIR)/testhnf : $(TESTHNF_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/testhnf.exe $(TESTHNF_OBJS) $(LIBS)

cleantesthnf: 
	rm -rf $(OUTPUTDIR)/testhnf.exe $(TESTHNF_OBJS)

TESTFACTOR_CPP_SRCS = testfactor.cpp VeryLong.cpp VeryLongModular.cpp squfof.cpp qs.cpp SparseMatrix.cpp timings.cpp blockLanczos.cpp

C_SRCS = mt19937int.c lip.c

include $(TESTFACTOR_CPP_SRCS:.cpp=.d)
include $(C_SRCS:.c=.d)

TESTFACTOR_OBJS = $(TESTFACTOR_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

testfactor : $(OUTPUTDIR) $(OUTPUTDIR)/testfactor
	

$(OUTPUTDIR)/testfactor : $(TESTFACTOR_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/testfactor.exe $(TESTFACTOR_OBJS) $(LIBS)

cleantestfactor: 
	rm -rf $(OUTPUTDIR)/testfactor.exe $(TESTFACTOR_OBJS)


POLYNOMIAL_UT_CPP_SRCS = PolynomialUnitTest.cpp Polynomial.cpp VeryLong.cpp VeryLongModular.cpp LongModular.cpp qs.cpp squfof.cpp timings.cpp SparseMatrix.cpp

include $(POLYNOMIAL_UT_CPP_SRCS:.cpp=.d)

POLYNOMIAL_UT_OBJS = $(POLYNOMIAL_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

PolynomialUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/PolynomialUnitTest
	

$(OUTPUTDIR)/PolynomialUnitTest : $(POLYNOMIAL_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/PolynomialUnitTest.exe $(POLYNOMIAL_UT_OBJS) $(TEST_LIBS)

cleanpolynomialut: 
	rm -rf $(OUTPUTDIR)/PolyinomialUnitTest.exe $(POLYNOMIAL_UT_OBJS)


POLYNOMIAL_OPTIMIZER_UT_CPP_SRCS = PolynomialOptimizerUnitTest.cpp PolynomialOptimizer.cpp pselect.cpp Polynomial.cpp VeryLong.cpp VeryLongModular.cpp LongModular.cpp qs.cpp squfof.cpp timings.cpp SparseMatrix.cpp discriminant.cpp dickman.cpp

include $(POLYNOMIAL_OPTIMIZER_UT_CPP_SRCS:.cpp=.d)

POLYNOMIAL_OPTIMIZER_UT_OBJS = $(POLYNOMIAL_OPTIMIZER_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

PolynomialOptimizerUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/PolynomialOptimizerUnitTest
	

$(OUTPUTDIR)/PolynomialOptimizerUnitTest : $(POLYNOMIAL_OPTIMIZER_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/PolynomialOptimizerUnitTest.exe $(POLYNOMIAL_OPTIMIZER_UT_OBJS) $(TEST_LIBS)

cleanpolynomialoptimizerut: 
	rm -rf $(OUTPUTDIR)/PolyinomialOptimizerUnitTest.exe $(POLYNOMIAL_OPTIMIZER_UT_OBJS)


LATTICE_SIEVER_UT_CPP_SRCS = LatticeSieverUnitTest.cpp LatticeSiever.cpp FactorBase.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp timings.cpp lll.cpp squfof.cpp qs.cpp SparseMatrix.cpp SieveUtils.cpp convert.cpp

include $(LATTICE_SIEVER_UT_CPP_SRCS:.cpp=.d)

LATTICE_SIEVER_UT_OBJS = $(LATTICE_SIEVER_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

LatticeSieverUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/LatticeSieverUnitTest
	

$(OUTPUTDIR)/LatticeSieverUnitTest : $(LATTICE_SIEVER_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/LatticeSieverUnitTest.exe $(LATTICE_SIEVER_UT_OBJS) $(TEST_LIBS)

cleanlatticesieverut: 
	rm -rf $(OUTPUTDIR)/LatticeSieverUnitTest.exe $(LATTICE_SIEVER_UT_OBJS)

MATRIX_UT_CPP_SRCS = MatrixUnitTest.cpp 

include $(MATRIX_UT_CPP_SRCS:.cpp=.d)

MATRIX_UT_OBJS = $(MATRIX_UT_CPP_SRCS:.cpp=.o) 

MatrixUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/MatrixUnitTest
	

$(OUTPUTDIR)/MatrixUnitTest : $(MATRIX_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/MatrixUnitTest.exe $(MATRIX_UT_OBJS) $(TEST_LIBS)

cleanmatrixut: 
	rm -rf $(OUTPUTDIR)/MatrixUnitTest.exe $(MATRIX_UT_OBJS)

SPARSE_MATRIX_UT_CPP_SRCS = SparseMatrixUnitTest.cpp SparseMatrix.cpp

include $(SPARSE_MATRIX_UT_CPP_SRCS:.cpp=.d)

SPARSE_MATRIX_UT_OBJS = $(SPARSE_MATRIX_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

SparseMatrixUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/SparseMatrixUnitTest
	

$(OUTPUTDIR)/SparseMatrixUnitTest : $(SPARSE_MATRIX_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/SparseMatrixUnitTest.exe $(SPARSE_MATRIX_UT_OBJS) $(TEST_LIBS)

cleansparsematrixut: 
	rm -rf $(OUTPUTDIR)/SparseMatrixUnitTest.exe $(SPARSE_MATRIX_UT_OBJS)

COMBINATIONS_UT_CPP_SRCS = CombinationsUnitTest.cpp 

include $(COMBINATIONS_UT_CPP_SRCS:.cpp=.d)

COMBINATIONS_UT_OBJS = $(COMBINATIONS_UT_CPP_SRCS:.cpp=.o) 

CombinationsUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/CombinationsUnitTest
	

$(OUTPUTDIR)/CombinationsUnitTest : $(COMBINATIONS_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/CombinationsUnitTest.exe $(COMBINATIONS_UT_OBJS) $(TEST_LIBS)

cleancombinationsut: 
	rm -rf $(OUTPUTDIR)/CombinationsUnitTest.exe $(COMBINATIONS_UT_OBJS)

MMF_UT_CPP_SRCS = MemoryMappedFileUnitTest.cpp 

include $(MMF_UT_CPP_SRCS:.cpp=.d)

MMF_UT_OBJS = $(MMF_UT_CPP_SRCS:.cpp=.o) 

MemoryMappedFileUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/MemoryMappedFileUnitTest
	

$(OUTPUTDIR)/MemoryMappedFileUnitTest : $(MMF_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/MemoryMappedFileUnitTest.exe $(MMF_UT_OBJS) $(TEST_LIBS)

cleanmmfut: 
	rm -rf $(OUTPUTDIR)/MemoryMappedFileUnitTest.exe $(MMF_UT_OBJS)

HASHTABLE_UT_CPP_SRCS = HashTableUnitTest.cpp 

include $(HASHTABLE_UT_CPP_SRCS:.cpp=.d)

HASHTABLE_UT_OBJS = $(HASHTABLE_UT_CPP_SRCS:.cpp=.o) 

HashTableUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/HashTableUnitTest
	

$(OUTPUTDIR)/HashTableUnitTest : $(HASHTABLE_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/HashTableUnitTest.exe $(HASHTABLE_UT_OBJS) $(TEST_LIBS)

cleanhashtableut: 
	rm -rf $(OUTPUTDIR)/HashTableUnitTest.exe $(HASHTABLE_UT_OBJS)

QUOTIENT_UT_CPP_SRCS = QuotientUnitTest.cpp 

include $(QUOTIENT_UT_CPP_SRCS:.cpp=.d)

QUOTIENT_UT_OBJS = $(QUOTIENT_UT_CPP_SRCS:.cpp=.o) 

QuotientUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/QuotientUnitTest
	

$(OUTPUTDIR)/QuotientUnitTest : $(QUOTIENT_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/QuotientUnitTest.exe $(QUOTIENT_UT_OBJS) $(TEST_LIBS)

cleanquotienteut: 
	rm -rf $(OUTPUTDIR)/QuotientUnitTest.exe $(QUOTIENT_UT_OBJS)

AN_UT_CPP_SRCS = AlgebraicNumberUnitTest.cpp NumberField.cpp Ideal.cpp AlgebraicNumber.cpp AlgebraicNumber_in_O_pO.cpp FactorBase.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp discriminant.cpp lll.cpp root.cpp timings.cpp squfof.cpp SparseMatrix.cpp qs.cpp

AN_UT_OBJS = $(AN_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(AN_UT_CPP_SRCS:.cpp=.d)

AlgebraicNumberUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/AlgebraicNumberUnitTest
	

$(OUTPUTDIR)/AlgebraicNumberUnitTest : $(AN_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/AlgebraicNumberUnitTest.exe $(AN_UT_OBJS) $(TEST_LIBS)

cleananut: 
	rm -rf $(OUTPUTDIR)/AlgebraicNumberUnitTest.exe $(AN_UT_OBJS)

BO_UT_CPP_SRCS = BitOperationsUnitTest.cpp 

BO_UT_OBJS = $(BO_UT_CPP_SRCS:.cpp=.o)

include $(BO_UT_CPP_SRCS:.cpp=.d)

BitOperationsUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/BitOperationsUnitTest
	

$(OUTPUTDIR)/BitOperationsUnitTest : $(BO_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/BitOperationsUnitTest.exe $(BO_UT_OBJS) $(TEST_LIBS) 

cleanbout: 
	rm -rf $(OUTPUTDIR)/BitOperationsUnitTest.exe $(BO_UT_OBJS)

LM_UT_CPP_SRCS = LongModularUnitTest.cpp LongModular.cpp

LM_UT_OBJS = $(LM_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(LM_UT_CPP_SRCS:.cpp=.d)

LongModularUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/LongModularUnitTest
	

$(OUTPUTDIR)/LongModularUnitTest : $(LM_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/LongModularUnitTest.exe $(LM_UT_OBJS) $(TEST_LIBS)

cleanlmut: 
	rm -rf $(OUTPUTDIR)/LongModularUnitTest.exe $(LM_UT_OBJS)

VL_UT_CPP_SRCS = VeryLongUnitTest.cpp VeryLong.cpp squfof.cpp VeryLongModular.cpp LongModular.cpp qs.cpp timings.cpp SparseMatrix.cpp

VL_UT_OBJS = $(VL_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(VL_UT_CPP_SRCS:.cpp=.d)

VeryLongUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/VeryLongUnitTest
	

$(OUTPUTDIR)/VeryLongUnitTest : $(VL_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/VeryLongUnitTest.exe $(VL_UT_OBJS) $(TEST_LIBS)

cleanvlut: 
	rm -rf $(OUTPUTDIR)/VeryLongUnitTest.exe $(VL_UT_OBJS)

VLM_UT_CPP_SRCS = VeryLongModularUnitTest.cpp VeryLong.cpp VeryLongModular.cpp LongModular.cpp squfof.cpp qs.cpp timings.cpp SparseMatrix.cpp

VLM_UT_OBJS = $(VLM_UT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(VLM_UT_CPP_SRCS:.cpp=.d)

VeryLongModularUnitTest : $(OUTPUTDIR) $(OUTPUTDIR)/VeryLongModularUnitTest
	

$(OUTPUTDIR)/VeryLongModularUnitTest : $(VLM_UT_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/VeryLongModularUnitTest.exe $(VLM_UT_OBJS) $(TEST_LIBS)

cleanvlmut: 
	rm -rf $(OUTPUTDIR)/VeryLongModularUnitTest.exe $(VLM_UT_OBJS)

SKEWED_CPP_SRCS = skewed.cpp timings.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp discriminant.cpp pselect.cpp PolynomialOptimizer.cpp dickman.cpp squfof.cpp qs.cpp SparseMatrix.cpp

C_SRCS = mt19937int.c lip.c

include $(SKEWED_CPP_SRCS:.cpp=.d)
include $(C_SRCS:.c=.d)

SKEWED_OBJS = $(SKEWED_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

skewed : $(OUTPUTDIR) $(OUTPUTDIR)/skewed
	

$(OUTPUTDIR)/skewed : $(SKEWED_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/skewed $(SKEWED_OBJS) $(LIBS)

cleanskewed: 
	rm -rf $(OUTPUTDIR)/skewed.exe $(SKEWED_OBJS)

EF_CPP_SRCS = ef.cpp timings.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp discriminant.cpp pselect.cpp PolynomialOptimizer.cpp dickman.cpp squfof.cpp qs.cpp SparseMatrix.cpp

C_SRCS = mt19937int.c lip.c

include $(EF_CPP_SRCS:.cpp=.d)
include $(C_SRCS:.c=.d)

EF_OBJS = $(EF_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

ef : $(OUTPUTDIR) $(OUTPUTDIR)/ef
	

$(OUTPUTDIR)/ef : $(EF_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/ef $(EF_OBJS) $(LIBS)

cleanef: 
	rm -rf $(OUTPUTDIR)/ef.exe $(EF_OBJS)

LSIEVE_CPP_SRCS = lsieve.cpp LatticeSiever.cpp FactorBase.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp timings.cpp lll.cpp squfof.cpp qs.cpp SparseMatrix.cpp

LSIEVE_OBJS = $(LSIEVE_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(LSIEVE_CPP_SRCS:.cpp=.d)

lsieve : $(OUTPUTDIR) $(OUTPUTDIR)/lsieve
	

$(OUTPUTDIR)/lsieve : $(LSIEVE_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/lsieve $(LSIEVE_OBJS) $(LIBS)

cleanlsieve: 
	rm -rf $(OUTPUTDIR)/lsieve.exe $(LSIEVE_OBJS)

BM_CPP_SRCS = buildMatrix.cpp SparseMatrix.cpp timings.cpp RelationManager.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp QuadraticCharacters.cpp LongModular.cpp squfof.cpp qs.cpp convert.cpp

BM_OBJS = $(BM_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(BM_CPP_SRCS:.cpp=.d)

buildMatrix : $(OUTPUTDIR) $(OUTPUTDIR)/buildMatrix
	

$(OUTPUTDIR)/buildMatrix : $(BM_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/buildMatrix $(BM_OBJS) $(LIBS)

cleanbm:
	rm -rf $(OUTPUTDIR)/buildMatrix.exe $(BM_OBJS)

TM_CPP_SRCS = testMatrix.cpp SparseMatrix.cpp timings.cpp 

TM_OBJS = $(TM_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(TM_CPP_SRCS:.cpp=.d)

testMatrix : $(OUTPUTDIR) $(OUTPUTDIR)/testMatrix
	
$(OUTPUTDIR)/testMatrix : $(TM_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/testMatrix.exe $(TM_OBJS) $(LIBS)

cleantm:
	rm -rf $(OUTPUTDIR)/testMatrix.exe $(TM_OBJS)

BL_CPP_SRCS = bl.cpp SparseMatrix.cpp blockLanczos.cpp timings.cpp

BL_OBJS = $(BL_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(BL_CPP_SRCS:.cpp=.d)

bl : $(OUTPUTDIR) $(OUTPUTDIR)/bl
	

$(OUTPUTDIR)/bl : $(BL_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/bl $(BL_OBJS) $(LIBS)

cleanbl:
	rm -rf $(OUTPUTDIR)/bl.exe $(BL_OBJS)

CALCROOT_CPP_SRCS = calcroot.cpp AlgebraicNumber.cpp AlgebraicNumber_in_O_pO.cpp NumberField.cpp Ideal.cpp FactorBase.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp discriminant.cpp lll.cpp root.cpp timings.cpp squfof.cpp qs.cpp SparseMatrix.cpp

CALCROOT_OBJS = $(CALCROOT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(CALCROOT_CPP_SRCS:.cpp=.d)

calcroot : $(OUTPUTDIR) $(OUTPUTDIR)/calcroot
	

$(OUTPUTDIR)/calcroot : $(CALCROOT_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/calcroot $(CALCROOT_OBJS) $(LIBS)

cleancalcroot:
	rm -rf $(OUTPUTDIR)/calcroot.exe $(CALCROOT_OBJS)

TESTAN_CPP_SRCS = testan.cpp NumberField.cpp Ideal.cpp AlgebraicNumber.cpp AlgebraicNumber_in_O_pO.cpp FactorBase.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp discriminant.cpp lll.cpp root.cpp timings.cpp squfof.cpp SparseMatrix.cpp qs.cpp

TESTAN_OBJS = $(TESTAN_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(TESTAN_CPP_SRCS:.cpp=.d)

testan : $(OUTPUTDIR) $(OUTPUTDIR)/testan
	

$(OUTPUTDIR)/testan : $(TESTAN_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/testan $(TESTAN_OBJS) $(LIBS)

cleantestan :
	rm -rf $(OUTPUTDIR)/testan.exe $(TESTAN_OBJS)

SIEVE_CPP_SRCS = sieve.cpp Siever.cpp FactorBase.cpp LongModular.cpp VeryLong.cpp VeryLongModular.cpp Polynomial.cpp QuadraticCharacters.cpp discriminant.cpp timings.cpp lll.cpp squfof.cpp convert.cpp qs.cpp SparseMatrix.cpp

SIEVE_OBJS = $(SIEVE_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(SIEVE_CPP_SRCS:.cpp=.d)

sieve : $(OUTPUTDIR) $(OUTPUTDIR)/sieve
	

$(OUTPUTDIR)/sieve : $(SIEVE_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/sieve $(SIEVE_OBJS) $(LIBS)

cleansieve: 
	rm -rf $(OUTPUTDIR)/sieve.exe $(SIEVE_OBJS)

FILTER_CPP_SRCS = filter.cpp RelationManager.cpp SparseMatrix.cpp timings.cpp convert.cpp

FILTER_OBJS = $(FILTER_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(FILTER_CPP_SRCS:.cpp=.d)

filter : $(OUTPUTDIR) $(OUTPUTDIR)/filter
	

$(OUTPUTDIR)/filter : $(FILTER_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/filter $(FILTER_OBJS) $(LIBS)

cleanfilter: 
	rm -rf $(OUTPUTDIR)/filter.exe $(FILTER_OBJS)

PD_CPP_SRCS = processDependencies.cpp SparseMatrix.cpp timings.cpp convert.cpp

PD_OBJS = $(PD_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(PD_CPP_SRCS:.cpp=.d)

processDependencies : $(OUTPUTDIR) $(OUTPUTDIR)/processDependencies
	

$(OUTPUTDIR)/processDependencies : $(PD_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/processDependencies.exe $(PD_OBJS) $(LIBS)

cleanpd:
	rm -rf $(OUTPUTDIR)/processDependencies.exe $(PD_OBJS)

CONVERT_CPP_SRCS = convertmain.cpp convert.cpp

CONVERT_OBJS = $(CONVERT_CPP_SRCS:.cpp=.o) 

include $(CONVERT_CPP_SRCS:.cpp=.d)

convert : $(OUTPUTDIR) $(OUTPUTDIR)/convert
	

$(OUTPUTDIR)/convert : $(CONVERT_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/convert $(CONVERT_OBJS) $(LIBSNOGMP)

cleanconvert: 
	rm -rf $(OUTPUTDIR)/convert.exe $(CONVERT_OBJS)

processDependenciespl : $(OUTPUTDIR) $(OUTPUTDIR)/processDependencies.pl
	

$(OUTPUTDIR)/processDependencies.pl : processDependencies.pl
	cp processDependencies.pl $(OUTPUTDIR)
	chmod +x $(OUTPUTDIR)/processDependencies.pl

cleanpdpl:
	rm -rf $(OUTPUTDIR)/processDepenencies.pl

postProcessing : $(OUTPUTDIR) $(OUTPUTDIR)/postProcessing.pl
	

$(OUTPUTDIR)/postProcessing.pl : postProcessing.pl
	sed 's/vcbin/gbin/' < postProcessing.pl > $(OUTPUTDIR)/postProcessing.pl
	chmod +x $(OUTPUTDIR)/postProcessing.pl

cleanpp:
	rm -rf $(OUTPUTDIR)/postProcessing.pl

DRT_CPP_SRCS = DenseRowTest.cpp SparseMatrix.cpp timings.cpp

DRT_OBJS = $(DRT_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(DRT_CPP_SRCS:.cpp=.d)

DenseRowTest : $(OUTPUTDIR) $(OUTPUTDIR)/DenseRowTest
	
$(OUTPUTDIR)/DenseRowTest : $(DRT_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/DenseRowTest $(DRT_OBJS) $(LIBS)

cleandrt:
	rm -rf $(OUTPUTDIR)/DenseRowTest.exe $(DRT_OBJS)

GEN_CPP_SRCS = gen.cpp timings.cpp

GEN_OBJS = $(GEN_CPP_SRCS:.cpp=.o)

include $(GEN_CPP_SRCS:.cpp=.d)

gen : $(OUTPUTDIR) $(OUTPUTDIR)/gen
	
$(OUTPUTDIR)/gen : $(GEN_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/gen.exe $(GEN_OBJS) $(LIBS)

cleangen:
	rm -rf $(OUTPUTDIR)/gen.exe $(GEN_OBJS)

INV_CPP_SRCS = inv.cpp timings.cpp

INV_OBJS = $(INV_CPP_SRCS:.cpp=.o)

include $(INV_CPP_SRCS:.cpp=.d)

inv : $(OUTPUTDIR) $(OUTPUTDIR)/inv
	
$(OUTPUTDIR)/inv : $(INV_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/inv.exe $(INV_OBJS) $(LIBS)

cleaninv:
	rm -rf $(OUTPUTDIR)/inv.exe $(INV_OBJS)

TESTMMF_CPP_SRCS = testmmf.cpp

TESTMMF_OBJS = $(TESTMMF_CPP_SRCS:.cpp=.o)

include $(TESTMMF_CPP_SRCS:.cpp=.d)

testmmf : $(OUTPUTDIR) $(OUTPUTDIR)/testmmf
	
$(OUTPUTDIR)/testmmf : $(TESTMMF_OBJS) 
	$(CC) $(CCFLAGS) -o $(OUTPUTDIR)/testmmf.exe $(TESTMMF_OBJS) $(LIBS)

cleantestmmf :
	rm -rf $(OUTPUTDIR)/testmmf.exe $(TESTMMF_OBJS)

STEST_CPP_SRCS = stest.cpp VeryLong.cpp VeryLongModular.cpp squfof.cpp qs.cpp SparseMatrix.cpp timings.cpp

STEST_OBJS = $(STEST_CPP_SRCS:.cpp=.o) $(C_SRCS:.c=.o)

include $(STEST_CPP_SRCS:.cpp=.d)

stest : $(OUTPUTDIR) $(OUTPUTDIR)/stest
	
$(OUTPUTDIR)/stest : $(STEST_OBJS) 
	g++ $(CCFLAGS) -o $(OUTPUTDIR)/stest $(STEST_OBJS) $(LIBS)

cleanstest: 
	rm -rf $(OUTPUTDIR)/stest.exe $(STEST_OBJS)

clean:
	rm -rf $(OUTPUTDIR)/*.exe *.o *.d

