#include <iostream>
#include "Polynomial.h"
#include "Quotient.h"
#include "VeryLong.h"
#include <vector>
#include "pow.h"

// Code taken from LiDIA
// The following two functions were contributed
//  by Roland Dreier (dreier@math.berkeley.edu)

// Calculate resultant of aa and bb via subresultant algorithm.
// Algorithm cribbed verbatim from Algorithm 3.3.7 of H. Cohen's "A
// course in computational algebraic number theory."  (Even the
// variables are named pretty much the same as in his book!)
//
// Author      : Roland Dreier (dreier@math.berkeley.edu)
//

VeryLong
resultant(const Polynomial <VeryLong> &aa,
          const Polynomial <VeryLong> &bb)
{

    const VeryLong zero(0L);
    const VeryLong one(1L);
    // Return zero if either polynomial is zero.
    if ((aa == Polynomial<VeryLong>(0L)) || (bb == Polynomial<VeryLong>(0L)))
    {
        return zero;
    }
    // otherwise...

    // Initialization steps:

    VeryLong acont, bcont;
    Polynomial <VeryLong> a, b;
    bool neg = false;

    // Maybe skip one reduction by making sure deg(a) >= deg(b).
    if (aa.deg() >= bb.deg())
    {
        a = aa;
        b = bb;
    }
    else
    {
        a = bb;
        b = aa;

        // Possibly change sign!!
        neg = ((a.deg() % 2) && (b.deg() % 2));
    }
    //cout << "resultant of " << a << " and " << b << endl;

    acont = a.content();
    bcont = b.content();

    a.make_primitive();
    b.make_primitive();

    VeryLong g = 1L;
    VeryLong h = 1L;

    VeryLong t;
    VeryLong pow_temp;

    //t = exp(VeryLong(acont), VeryLong((long int)b.deg()));
    t = pow<VeryLong, long int>(acont, (long int)b.deg());
    //pow_temp = exp(VeryLong(bcont), VeryLong((long int)a.deg()));
    pow_temp = pow<VeryLong, long int>(bcont, (long int)a.deg());
    t = t * pow_temp;

    // Loop over pseudo-division and reduction steps:

    long int delta;
    Polynomial <VeryLong> r;

    do
    {
        delta = a.deg() - b.deg();

        if ((a.deg() % 2) && (b.deg() % 2))
        {
            neg = !neg;
        }

        //cout << "about to calculate remainder(" << a << "," << b << ")" << endl;

        r = remainder(a, b);
        a = b;

        //pow_temp = exp(h, VeryLong(delta));
        pow_temp = pow<VeryLong, long int>(h, delta);
        //cout << "1. r = " << r << ", pow_temp = " << pow_temp << ", g = " << g << ", h = " << h << endl;
        b = r / pow_temp;
        b = b / g;

        //cout << "b = " << b << endl;

        //cout << "Top coeff of a = " << a.coefficient(a.deg()) << endl;
        g = a.coefficient(a.deg());
        //pow_temp = exp(g, VeryLong(delta--));
        pow_temp = pow<VeryLong, long int>(g, delta--);
        //cout << "2. r = " << r << ", pow_temp = " << pow_temp << ", g = " << g << ", h = " << h << endl;

        if (delta<=0)
        {
            //h = exp(h, VeryLong(-delta));
            h = pow<VeryLong, long int>(h, -delta);
            h = h * pow_temp;

        }
        else
        {
            //h = exp(h, VeryLong(delta));
            h = pow<VeryLong, long int>(h, delta);
            h = pow_temp / h;
        }
        //cout << "3. r = " << r << ", pow_temp = " << pow_temp << ", g = " << g << ", h = " << h << endl;
    }
    while (b.deg() > 0);

    // Finish up:

    VeryLong temp = b.coefficient(b.deg());
    //cout << "temp = " << temp << endl;
    //pow_temp = exp(temp, VeryLong((long int)a.deg()));
    pow_temp = pow<VeryLong, long int>(temp, (long int)a.deg());

    delta = a.deg()-1;
    if (delta <= 0)
    {
        //h = exp(h, VeryLong(-delta));
        h = pow<VeryLong, long int>(h, -delta);
        h = h * pow_temp;
    }
    else
    {
        //h = exp(h,VeryLong(delta));
        h = pow<VeryLong, long int>(h,delta);
        h = pow_temp / h;
    }
    if (neg)
        t = -one * t;
    VeryLong res = t*h;
    //cout << "Resultant = " << res << endl;
    return res;
}

//
// Calculate discriminant of a polynomial using the formula:
//   disc(a) = (-1)^(d*(d-1)/2) * resultant(a, a') / lead_coeff(a)
// where d = deg(f).
//
// Rather than raising -1 to a power, just look at d mod 4.
// Remember, d*(d-1)/2 is even iff d = 0 or 1 mod 4.
//
// Author      : Roland Dreier (dreier@math.berkeley.edu)
//
VeryLong
discriminant(const Polynomial <VeryLong> &a)
{
//	cout << "discriminant of " << a << endl;
    if (a.deg() <= 0)
    {
        return(VeryLong(0L));
    }

    //VeryLong multiplier = exp(a.coefficient(a.deg()), VeryLong(a.deg()*2L - 3L));
    VeryLong multiplier = pow<VeryLong, long int>(a.coefficient(a.deg()), a.deg()*2L - 3L);
    if ((int(a.deg()) % 4) <= 1)
    {
        return(resultant(a, a.derivative()) / a.coefficient(a.deg()));
    }
    else
    {
        return(-(resultant(a, a.derivative()) / a.coefficient(a.deg())));
    }
}
