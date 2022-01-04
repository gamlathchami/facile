// below functions are taken from jStat library(https://github.com/jstat/jstat)
// and are slightly modified to fit this project.

/* 
Copyright (c) 2013 jStat

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
*/

function betacf(x, a, b) {
    let fpmin = 1e-30;
    let m = 1;
    let qab = a + b;
    let qap = a + 1;
    let qam = a - 1;
    let c = 1;
    let d = 1 - qab * x / qap;
    let m2, aa, del, h;
    
    // These q's will be used in factors that occur in the coefficients
    if (Math.abs(d) < fpmin)
    d = fpmin;
    d = 1 / d;
    h = d;
    
    for (; m <= 100; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        // One step (the even one) of the recurrence
        d = 1 + aa * d;
        if (Math.abs(d) < fpmin)
        d = fpmin;
        c = 1 + aa / c;
        if (Math.abs(c) < fpmin)
        c = fpmin;
        d = 1 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        // Next step of the recurrence (the odd one)
        d = 1 + aa * d;
        if (Math.abs(d) < fpmin)
        d = fpmin;
        c = 1 + aa / c;
        if (Math.abs(c) < fpmin)
        c = fpmin;
        d = 1 / d;
        del = d * c;
        h *= del;
        if (Math.abs(del - 1.0) < 3e-7)
        break;
    }
    
    return h;
};


function gammaln(x) {
    let j = 0;
    let cof = [
        76.18009172947146, -86.50532032941677, 24.01409824083091,
        -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
    ];
    let ser = 1.000000000190015;
    let xx, y, tmp;
    tmp = (y = xx = x) + 5.5;
    tmp -= (xx + 0.5) * Math.log(tmp);
    for (; j < 6; j++)
    ser += cof[j] / ++y;
    return Math.log(2.5066282746310005 * ser / xx) - tmp;
};


function ibeta(x, a, b) {
    // Factors in front of the continued fraction.
    let bt = (x === 0 || x === 1) ?  0 :
    Math.exp(gammaln(a + b) - gammaln(a) -
    gammaln(b) + a * Math.log(x) + b *
    Math.log(1 - x));
    if (x < 0 || x > 1)
    return false;
    if (x < (a + 1) / (a + b + 2))
    // Use continued fraction directly.
    return bt * betacf(x, a, b) / a;
    // else use continued fraction after making the symmetry transformation.
    return 1 - bt * betacf(1 - x, b, a) / b;
};


function ibetainv(p, a, b) {
    let EPS = 1e-8;
    let a1 = a - 1;
    let b1 = b - 1;
    let j = 0;
    let lna, lnb, pp, t, u, err, x, al, h, w, afac;
    if (p <= 0)
    return 0;
    if (p >= 1)
    return 1;
    if (a >= 1 && b >= 1) {
        pp = (p < 0.5) ? p : 1 - p;
        t = Math.sqrt(-2 * Math.log(pp));
        x = (2.30753 + t * 0.27061) / (1 + t* (0.99229 + t * 0.04481)) - t;
        if (p < 0.5)
        x = -x;
        al = (x * x - 3) / 6;
        h = 2 / (1 / (2 * a - 1)  + 1 / (2 * b - 1));
        w = (x * Math.sqrt(al + h) / h) - (1 / (2 * b - 1) - 1 / (2 * a - 1)) *
        (al + 5 / 6 - 2 / (3 * h));
        x = a / (a + b * Math.exp(2 * w));
    } else {
        lna = Math.log(a / (a + b));
        lnb = Math.log(b / (a + b));
        t = Math.exp(a * lna) / a;
        u = Math.exp(b * lnb) / b;
        w = t + u;
        if (p < t / w)
        x = Math.pow(a * w * p, 1 / a);
        else
        x = 1 - Math.pow(b * w * (1 - p), 1 / b);
    }
    afac = -gammaln(a) - gammaln(b) + gammaln(a + b);
    for(; j < 10; j++) {
        if (x === 0 || x === 1)
        return x;
        err = ibeta(x, a, b) - p;
        t = Math.exp(a1 * Math.log(x) + b1 * Math.log(1 - x) + afac);
        u = err / t;
        x -= (t = u / (1 - 0.5 * Math.min(1, u * (a1 / x - b1 / (1 - x)))));
        if (x <= 0)
        x = 0.5 * (x + t);
        if (x >= 1)
        x = 0.5 * (x + t + 1);
        if (Math.abs(t) < EPS * x && j > 0)
        break;
    }
    return x;
};


function centralFInv(x, df1, df2) {
    return df2 / (df1 * (1 / ibetainv(x, df1 / 2, df2 / 2) - 1));
}


function tukeyQinv(p, c, v) {
    let p0 = 0.322232421088;
    let q0 = 0.993484626060e-01;
    let p1 = -1.0;
    let q1 = 0.588581570495;
    let p2 = -0.342242088547;
    let q2 = 0.531103462366;
    let p3 = -0.204231210125;
    let q3 = 0.103537752850;
    let p4 = -0.453642210148e-04;
    let q4 = 0.38560700634e-02;
    let c1 = 0.8832;
    let c2 = 0.2368;
    let c3 = 1.214;
    let c4 = 1.208;
    let c5 = 1.4142;
    let vmax = 120.0;
  
    let ps = 0.5 - 0.5 * p;
    let yi = Math.sqrt(Math.log(1.0 / (ps * ps)));
    let t = yi + (((( yi * p4 + p3) * yi + p2) * yi + p1) * yi + p0)
       / (((( yi * q4 + q3) * yi + q2) * yi + q1) * yi + q0);
    if (v < vmax) t += (t * t * t + t) / v / 4.0;
    let q = c1 - c2 * t;
    if (v < vmax) q += -c3 / v + c4 * t / v;
    return t * (q * Math.log(c - 1.0) + c5);
}


function erf(x) {
    let cof = [-1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2,
               -9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
               4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
               1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
               6.529054439e-9, 5.059343495e-9, -9.91364156e-10,
               -2.27365122e-10, 9.6467911e-11, 2.394038e-12,
               -6.886027e-12, 8.94487e-13, 3.13092e-13,
               -1.12708e-13, 3.81e-16, 7.106e-15,
               -1.523e-15, -9.4e-17, 1.21e-16,
               -2.8e-17];
    let j = cof.length - 1;
    let isneg = false;
    let d = 0;
    let dd = 0;
    let t, ty, tmp, res;
  
    if (x < 0) {
      x = -x;
      isneg = true;
    }
  
    t = 2 / (2 + x);
    ty = 4 * t - 2;
  
    for(; j > 0; j--) {
      tmp = d;
      d = ty * d - dd + cof[j];
      dd = tmp;
    }
  
    res = t * Math.exp(-x * x + 0.5 * (cof[0] + ty * d) - dd);
    return isneg ? res - 1 : 1 - res;
}

function normalCdf(x, mean, std) {
    return 0.5 * (1 + erf((x - mean) / Math.sqrt(2 * std * std)));
}

function tukeyWprob(w, rr, cc) {
    let nleg = 12;
    let ihalf = 6;
  
    let C1 = -30;
    let C2 = -50;
    let C3 = 60;
    let bb   = 8;
    let wlar = 3;
    let wincr1 = 2;
    let wincr2 = 3;
    let xleg = [
      0.981560634246719250690549090149,
      0.904117256370474856678465866119,
      0.769902674194304687036893833213,
      0.587317954286617447296702418941,
      0.367831498998180193752691536644,
      0.125233408511468915472441369464
    ];
    let aleg = [
      0.047175336386511827194615961485,
      0.106939325995318430960254718194,
      0.160078328543346226334652529543,
      0.203167426723065921749064455810,
      0.233492536538354808760849898925,
      0.249147045813402785000562436043
    ];
  
    let qsqz = w * 0.5;
  
    // if w >= 16 then the integral lower bound (occurs for c=20)
    // is 0.99999999999995 so return a value of 1.
  
    if (qsqz >= bb)
      return 1.0;
  
    // find (f(w/2) - 1) ^ cc
    // (first term in integral of hartley's form).
  
    let pr_w = 2 * normalCdf(qsqz, 0, 1, 1, 0) - 1; // erf(qsqz / M_SQRT2)
    // if pr_w ^ cc < 2e-22 then set pr_w = 0
    if (pr_w >= Math.exp(C2 / cc))
      pr_w = Math.pow(pr_w, cc);
    else
      pr_w = 0.0;
  
    // if w is large then the second component of the
    // integral is small, so fewer intervals are needed.
  
    let wincr;
    if (w > wlar)
      wincr = wincr1;
    else
      wincr = wincr2;
  
    // find the integral of second term of hartley's form
    // for the integral of the range for equal-length
    // intervals using legendre quadrature.  limits of
    // integration are from (w/2, 8).  two or three
    // equal-length intervals are used.
  
    // blb and bub are lower and upper limits of integration.
  
    let blb = qsqz;
    let binc = (bb - qsqz) / wincr;
    let bub = blb + binc;
    let einsum = 0.0;
  
    // integrate over each interval
  
    let cc1 = cc - 1.0;
    for (let wi = 1; wi <= wincr; wi++) {
      let elsum = 0.0;
      let a = 0.5 * (bub + blb);
  
      // legendre quadrature with order = nleg
  
      let b = 0.5 * (bub - blb);
  
      for (let jj = 1; jj <= nleg; jj++) {
        let j, xx;
        if (ihalf < jj) {
          j = (nleg - jj) + 1;
          xx = xleg[j-1];
        } else {
          j = jj;
          xx = -xleg[j-1];
        }
        let c = b * xx;
        let ac = a + c;
  
        // if exp(-qexpo/2) < 9e-14,
        // then doesn't contribute to integral
  
        let qexpo = ac * ac;
        if (qexpo > C3)
          break;
  
        let pplus = 2 * normalCdf(ac, 0, 1, 1, 0);
        let pminus= 2 * normalCdf(ac, w, 1, 1, 0);
  
        // if rinsum ^ (cc-1) < 9e-14,
        // then doesn't contribute to integral
  
        let rinsum = (pplus * 0.5) - (pminus * 0.5);
        if (rinsum >= Math.exp(C1 / cc1)) {
          rinsum = (aleg[j-1] * Math.exp(-(0.5 * qexpo))) * Math.pow(rinsum, cc1);
          elsum += rinsum;
        }
      }
      elsum *= (((2.0 * b) * cc) / Math.sqrt(2 * Math.PI));
      einsum += elsum;
      blb = bub;
      bub += binc;
    }
  
    // if pr_w ^ rr < 9e-14, then return 0
    pr_w += einsum;
    if (pr_w <= Math.exp(C1 / rr))
      return 0;
  
    pr_w = Math.pow(pr_w, rr);
    if (pr_w >= 1) // 1 was iMax was eps
      return 1;
    return pr_w;
}
  

function tukeyCdf(q, nmeans, df) {
    // Identical implementation as the R ptukey() function as of commit 68947
    let rr = 1;
    let cc = nmeans;

    let nlegq = 16;
    let ihalfq = 8;

    let eps1 = -30.0;
    let eps2 = 1.0e-14;
    let dhaf  = 100.0;
    let dquar = 800.0;
    let deigh = 5000.0;
    let dlarg = 25000.0;
    let ulen1 = 1.0;
    let ulen2 = 0.5;
    let ulen3 = 0.25;
    let ulen4 = 0.125;
    let xlegq = [
      0.989400934991649932596154173450,
      0.944575023073232576077988415535,
      0.865631202387831743880467897712,
      0.755404408355003033895101194847,
      0.617876244402643748446671764049,
      0.458016777657227386342419442984,
      0.281603550779258913230460501460,
      0.950125098376374401853193354250e-1
    ];
    let alegq = [
      0.271524594117540948517805724560e-1,
      0.622535239386478928628438369944e-1,
      0.951585116824927848099251076022e-1,
      0.124628971255533872052476282192,
      0.149595988816576732081501730547,
      0.169156519395002538189312079030,
      0.182603415044923588866763667969,
      0.189450610455068496285396723208
    ];

    if (q <= 0)
      return 0;

    // df must be > 1
    // there must be at least two values

    if (df < 2 || rr < 1 || cc < 2) return NaN;

    if (!Number.isFinite(q))
      return 1;

    if (df > dlarg)
      return tukeyWprob(q, rr, cc);

    // calculate leading constant

    let f2 = df * 0.5;
    let f2lf = ((f2 * Math.log(df)) - (df * Math.log(2))) - gammaln(f2);
    let f21 = f2 - 1.0;

    // integral is divided into unit, half-unit, quarter-unit, or
    // eighth-unit length intervals depending on the value of the
    // degrees of freedom.

    let ff4 = df * 0.25;
    let ulen;
    if      (df <= dhaf)  ulen = ulen1;
    else if (df <= dquar) ulen = ulen2;
    else if (df <= deigh) ulen = ulen3;
    else                  ulen = ulen4;

    f2lf += Math.log(ulen);

    // integrate over each subinterval

    let ans = 0.0;

    for (let i = 1; i <= 50; i++) {
      var otsum = 0.0;

      // legendre quadrature with order = nlegq
      // nodes (stored in xlegq) are symmetric around zero.

      let twa1 = (2 * i - 1) * ulen;

      for (let jj = 1; jj <= nlegq; jj++) {
        let j, t1;
        if (ihalfq < jj) {
          j = jj - ihalfq - 1;
          t1 = (f2lf + (f21 * Math.log(twa1 + (xlegq[j] * ulen))))
              - (((xlegq[j] * ulen) + twa1) * ff4);
        } else {
          j = jj - 1;
          t1 = (f2lf + (f21 * Math.log(twa1 - (xlegq[j] * ulen))))
              + (((xlegq[j] * ulen) - twa1) * ff4);
        }

        // if exp(t1) < 9e-14, then doesn't contribute to integral
        let qsqz;
        if (t1 >= eps1) {
          if (ihalfq < jj) {
            qsqz = q * Math.sqrt(((xlegq[j] * ulen) + twa1) * 0.5);
          } else {
            qsqz = q * Math.sqrt(((-(xlegq[j] * ulen)) + twa1) * 0.5);
          }

          // call wprob to find integral of range portion

          let wprb = tukeyWprob(qsqz, rr, cc);
          let rotsum = (wprb * alegq[j]) * Math.exp(t1);
          otsum += rotsum;
        }
        // end legendre integral for interval i
        // L200:
      }

      // if integral for interval i < 1e-14, then stop.
      // However, in order to avoid small area under left tail,
      // at least  1 / ulen  intervals are calculated.
      if (i * ulen >= 1.0 && otsum <= eps2)
        break;

      // end of interval i
      // L330:

      ans += otsum;
    }

    if (otsum > eps2) { // not converged
      throw new Error('tukey.cdf failed to converge');
    }
    if (ans > 1)
      ans = 1;
    return ans;
}


function tukeyInv(p, nmeans, df) {
    // Identical implementation as the R qtukey() function as of commit 68947
    let rr = 1;
    let cc = nmeans;

    let eps = 0.0001;
    let maxiter = 50;

    // df must be > 1 ; there must be at least two values
    if (df < 2 || rr < 1 || cc < 2) return NaN;

    if (p < 0 || p > 1) return NaN;
    if (p === 0) return 0;
    if (p === 1) return Infinity;

    // Initial value

    let x0 = tukeyQinv(p, cc, df);

    // Find prob(value < x0)

    let valx0 = tukeyCdf(x0, nmeans, df) - p;

    // Find the second iterate and prob(value < x1).
    // If the first iterate has probability value
    // exceeding p then second iterate is 1 less than
    // first iterate; otherwise it is 1 greater.

    let x1;
    if (valx0 > 0.0)
      x1 = Math.max(0.0, x0 - 1.0);
    else
      x1 = x0 + 1.0;
    let valx1 = tukeyCdf(x1, nmeans, df) - p;

    // Find new iterate

    let ans;
    for(let iter = 1; iter < maxiter; iter++) {
      ans = x1 - ((valx1 * (x1 - x0)) / (valx1 - valx0));
      valx0 = valx1;

      // New iterate must be >= 0

      x0 = x1;
      if (ans < 0.0) {
        ans = 0.0;
        valx1 = -p;
      }
      // Find prob(value < new iterate)

      valx1 = tukeyCdf(ans, nmeans, df) - p;
      x1 = ans;

      // If the difference between two successive
      // iterates is less than eps, stop

      let xabs = Math.abs(x1 - x0);
      if (xabs < eps)
        return ans;
    }

    throw new Error('tukeyInv failed to converge');
}


export { centralFInv, tukeyInv };
