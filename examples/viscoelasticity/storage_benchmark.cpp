// Storage and runtime benchmark for backward access to Maxwell history variables.
//
// Np independent material points carry the SA1 Maxwell update
//   eps_v <- A*eps_v + (1-A)*eps(t_{n+1})
// under prescribed strain histories. Three strategies give the backward sweep
// access to every forward state eps_v(t_n), n = N,...,0:
//   store_all   : floating-point forward, every state stored
//   checkpoint  : floating-point forward, uniform sqrt(N) checkpoints, segments recomputed
//   reversible  : fixed-point forward (Algorithm 3), states rebuilt by the inverse update,
//                 ancillary contents flushed every L steps
// A fourth run, explicit_inverse, applies the floating-point inverse without storage.
//
// The backward sweep "consumes" each state by accumulating a checksum, so every
// strategy performs the same work per recovered state.
//
// Build: g++ -O3 -march=native -ffp-contract=off -std=c++17 -I../../src storage_benchmark.cpp -o storage_benchmark

#include "types.h"
#include "Fixed.h"
#include "Rational.h"
#include "Dual.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

typedef Dual<Fixed_E> DualE;

struct Problem {
  int    Np;        // number of scalar history variables (material points x components)
  int    N;         // number of time steps
  double dt;        // time step
  double tau;       // relaxation time
  int    L;         // ancillary flush interval (steps)
  std::vector<int>    kind;   // loading history per point
  std::vector<double> scale;  // amplitude scaling per point
};

// Prescribed strain histories of Section 7, scaled per point
inline double strain(const Problem& P, int i, int n) {
  const double t = n*P.dt;
  const double s = P.scale[i];
  switch (P.kind[i]) {
    case 0:  return (n == 0) ? 0.0 : 0.1*s;                                  // step
    case 1:  return 0.1*s*t;                                                 // constant rate
    case 2:  return 0.1*s*std::sin(0.2*M_PI*t);                              // sinusoidal
    default: return s*std::max(-0.1,std::min(0.1,0.2*std::sin(0.2*M_PI*t))); // clipped sinusoid
  }
}

typedef std::chrono::steady_clock Clock;
inline double seconds_since(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

struct Result {
  std::string name;
  double t_forward = 0.0, t_backward = 0.0;
  long   bytes = 0;            // stored for backward access (beyond the working state)
  long   extra_forward_steps = 0;
  double checksum = 0.0;       // sum of consumed states
  double max_error = 0.0;      // max |recovered - forward| (explicit inverse only)
  long   flushes = 0;
  long   max_ancillary = 0;    // max |ancillary mantissa| reached
};

// ---------------------------------------------------------------- store_all
Result run_store_all(const Problem& P) {
  Result R; R.name = "store_all";
  const double A = std::exp(-P.dt/P.tau);
  std::vector<double> x(P.Np,0.0);
  std::vector<double> hist(size_t(P.N+1)*P.Np);
  auto t0 = Clock::now();
  std::memcpy(&hist[0],x.data(),sizeof(double)*P.Np);
  for (int n=0; n<P.N; n++) {
    double* row = &hist[size_t(n+1)*P.Np];
    for (int i=0; i<P.Np; i++) {
      x[i] = A*x[i] + (1.0-A)*strain(P,i,n+1);
      row[i] = x[i];
    }
  }
  R.t_forward = seconds_since(t0);
  t0 = Clock::now();
  double cs = 0.0;
  for (int n=P.N; n>=0; n--) {
    const double* row = &hist[size_t(n)*P.Np];
    for (int i=0; i<P.Np; i++) cs += row[i];
  }
  R.t_backward = seconds_since(t0);
  R.checksum = cs;
  R.bytes = long(hist.size()*sizeof(double));
  return R;
}

// --------------------------------------------------------------- checkpoint
Result run_checkpoint(const Problem& P) {
  Result R; R.name = "checkpoint";
  const double A = std::exp(-P.dt/P.tau);
  const int k = int(std::ceil(std::sqrt(double(P.N))));
  const int Nseg = (P.N + k - 1)/k;
  std::vector<double> x(P.Np,0.0);
  std::vector<double> ckpt(size_t(Nseg)*P.Np);   // states at n = 0, k, 2k, ...
  std::vector<double> buf(size_t(k)*P.Np);       // recomputed segment
  auto t0 = Clock::now();
  for (int n=0; n<P.N; n++) {
    if (n % k == 0) std::memcpy(&ckpt[size_t(n/k)*P.Np],x.data(),sizeof(double)*P.Np);
    for (int i=0; i<P.Np; i++) x[i] = A*x[i] + (1.0-A)*strain(P,i,n+1);
  }
  R.t_forward = seconds_since(t0);
  t0 = Clock::now();
  double cs = 0.0;
  for (int i=0; i<P.Np; i++) cs += x[i];              // state N is the working state
  for (int s=Nseg-1; s>=0; s--) {
    const int n0 = s*k, n1 = std::min(n0+k,P.N);      // segment covers states n0..n1-1
    std::vector<double> y(&ckpt[size_t(s)*P.Np],&ckpt[size_t(s)*P.Np]+P.Np);
    std::memcpy(&buf[0],y.data(),sizeof(double)*P.Np);
    for (int n=n0; n<n1-1; n++) {
      double* row = &buf[size_t(n-n0+1)*P.Np];
      for (int i=0; i<P.Np; i++) { y[i] = A*y[i] + (1.0-A)*strain(P,i,n+1); row[i] = y[i]; }
      R.extra_forward_steps++;
    }
    for (int n=n1-1; n>=n0; n--) {
      const double* row = &buf[size_t(n-n0)*P.Np];
      for (int i=0; i<P.Np; i++) cs += row[i];
    }
  }
  R.t_backward = seconds_since(t0);
  R.checksum = cs;
  R.bytes = long((ckpt.size() + buf.size())*sizeof(double));
  return R;
}

// --------------------------------------------------------- explicit_inverse
Result run_explicit_inverse(const Problem& P, const std::vector<double>& ref_final_row0) {
  Result R; R.name = "explicit_inverse";
  const double A = std::exp(-P.dt/P.tau);
  std::vector<double> x(P.Np,0.0);
  std::vector<double> fwd0(P.Np,0.0);   // forward state at n = 0 (all zeros)
  auto t0 = Clock::now();
  for (int n=0; n<P.N; n++)
    for (int i=0; i<P.Np; i++) x[i] = A*x[i] + (1.0-A)*strain(P,i,n+1);
  R.t_forward = seconds_since(t0);
  t0 = Clock::now();
  double cs = 0.0;
  for (int i=0; i<P.Np; i++) cs += x[i];
  for (int n=P.N; n>=1; n--) {
    for (int i=0; i<P.Np; i++) { x[i] = (x[i] - (1.0-A)*strain(P,i,n))/A; cs += x[i]; }
  }
  R.t_backward = seconds_since(t0);
  for (int i=0; i<P.Np; i++) R.max_error = std::max(R.max_error,std::fabs(x[i]-fwd0[i]));
  R.checksum = cs;
  (void)ref_final_row0;
  return R;
}

// --------------------------------------------------------------- reversible
// Forward (Algorithm 3): squeeze by p/q, then add the truncated inhomogeneous term.
// Backward: subtract the same term, then apply the inverse squeeze q/p.
// The ancillary contents are flushed every L steps and restored in reverse order.
Result run_reversible(const Problem& P, std::vector<Integer>* fwd_mantissa) {
  Result R; R.name = "reversible";
  const double A = std::exp(-P.dt/P.tau);
  const Rational A_rat(A);
  const double one_minus_A = 1.0 - Real(A_rat);
  std::vector<DualE> x(P.Np,DualE(Fixed_E(0.0),Fixed_E(0.0)));
  std::vector<std::vector<Integer> > records;
  long max_anc = 0;
  int counter = 0;
  auto t0 = Clock::now();
  for (int n=0; n<P.N; n++) {
    for (int i=0; i<P.Np; i++) {
      DualE y = x[i]*A_rat;
      y = y + DualE(Fixed_E(strain(P,i,n+1)*one_minus_A),Fixed_E(0.0));
      x[i] = y;
      max_anc = std::max(max_anc,long(std::labs(long(y.second.mantissa))));
    }
    if (fwd_mantissa) for (int i=0; i<P.Np; i++) fwd_mantissa->push_back(x[i].first.mantissa);
    if (++counter == P.L) {
      counter = 0;
      std::vector<Integer> rec(P.Np);
      for (int i=0; i<P.Np; i++) { rec[i] = x[i].second.mantissa; x[i].second = Fixed_E(0.0); }
      records.push_back(std::move(rec));
    }
  }
  R.t_forward = seconds_since(t0);
  R.flushes = long(records.size());
  R.max_ancillary = max_anc;
  for (auto& rec : records) R.bytes += long(rec.size()*sizeof(Integer));
  R.bytes += long(P.Np*sizeof(Integer));   // persistent ancillary mantissa
  t0 = Clock::now();
  double cs = 0.0;
  long mismatches = 0;
  for (int i=0; i<P.Np; i++) cs += Real(x[i].first);
  for (int n=P.N; n>=1; n--) {
    if (counter == 0 && !records.empty()) {
      counter = P.L;
      std::vector<Integer>& rec = records.back();
      for (int i=0; i<P.Np; i++) x[i].second.mantissa = rec[i];
      records.pop_back();
    }
    counter--;
    for (int i=0; i<P.Np; i++) {
      DualE y = x[i] - DualE(Fixed_E(strain(P,i,n)*one_minus_A),Fixed_E(0.0));
      x[i] = y/A_rat;
      cs += Real(x[i].first);
    }
    if (fwd_mantissa && n-1 >= 1) {
      const Integer* ref = &(*fwd_mantissa)[size_t(n-2)*P.Np];
      for (int i=0; i<P.Np; i++) mismatches += (x[i].first.mantissa != ref[i]);
    }
  }
  R.t_backward = seconds_since(t0);
  for (int i=0; i<P.Np; i++) mismatches += (x[i].first.mantissa != 0) + (x[i].second.mantissa != 0);
  R.checksum = cs;
  R.max_error = double(mismatches);   // number of mismatched mantissas (verification run)
  return R;
}

int main(int argc, char** argv) {
  Problem P;
  P.Np  = (argc > 1) ? std::atoi(argv[1]) : 10368;
  P.N   = (argc > 2) ? std::atoi(argv[2]) : 15000;
  P.L   = (argc > 3) ? std::atoi(argv[3]) : 1000;
  const int repeats = (argc > 4) ? std::atoi(argv[4]) : 5;
  P.dt  = 1.0e-3;
  P.tau = 0.3;
  P.kind.resize(P.Np);
  P.scale.resize(P.Np);
  for (int i=0; i<P.Np; i++) {
    P.kind[i]  = i % 4;
    P.scale[i] = 0.5 + double((i*7919) % P.Np)/double(P.Np);   // in [0.5, 1.5)
  }

  // Verification run of the reversible update: every mantissa compared at every step
  std::vector<Integer> fwd;
  fwd.reserve(size_t(P.N)*P.Np);
  Result V = run_reversible(P,&fwd);
  std::printf("verify reversible: mismatched_mantissas=%.0f flushes=%ld max_ancillary=%ld int32_max=%d\n",
              V.max_error,V.flushes,V.max_ancillary,std::numeric_limits<Integer>::max());
  fwd.clear(); fwd.shrink_to_fit();

  std::vector<double> dummy;
  std::printf("method,repeat,Np,N,L,t_forward,t_backward,bytes,extra_forward_steps,checksum,max_error,flushes,max_ancillary\n");
  for (int r=0; r<repeats; r++) {
    Result rs[4] = { run_store_all(P), run_checkpoint(P), run_reversible(P,nullptr), run_explicit_inverse(P,dummy) };
    for (auto& R : rs) {
      std::printf("%s,%d,%d,%d,%d,%.6f,%.6f,%ld,%ld,%.17g,%.6e,%ld,%ld\n",
                  R.name.c_str(),r,P.Np,P.N,P.L,R.t_forward,R.t_backward,R.bytes,
                  R.extra_forward_steps,R.checksum,R.max_error,R.flushes,R.max_ancillary);
    }
  }
  return 0;
}
