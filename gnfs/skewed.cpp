#include "mt19937int.h"
#include "pselect.h"
#include <time.h>
#include <string.h>

int main(int argc, char** argv)
{

   sgenrand(time(0));
#if 0
   VeryLong stage10("10742788291266565907178411279942116612663921794753294588877817210355464150980121879033832926235281090750672083504941996433143425558334401855808989426892463");
   VeryLong RSA576("188198812920607963838697239461650439807163563379417382700763356422988859715234665485319060606504743045317388011303396716199692321205734031879550656996221305168759307650257059");
   VeryLong jim("62941029491309");
   VeryLong N = RSA576;
#endif
   skewed_polynomial_selection();

   return 0;
}
