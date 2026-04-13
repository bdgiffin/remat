#ifndef PASS_PHASE_H
#define PASS_PHASE_H

enum class PassPhase {
  Forward = 0,
  Backward = 1,
  BackwardAdjoint = 2
};

inline const char* pass_phase_name(PassPhase phase) {
  switch (phase) {
    case PassPhase::Forward:  return "Forward";
    case PassPhase::Backward:  return "Backward";
    case PassPhase::BackwardAdjoint: return "BackwardAdjoint";
  }
  return "Unknown";
}

inline bool is_reverse_phase(PassPhase phase) {
  return (phase == PassPhase::Backward) || (phase == PassPhase::BackwardAdjoint);
}

#endif // PASS_PHASE_H
