#include "random_backend.h"

namespace {
struct RandomBackendState {
  bool is_stub = false;
  const char* name = "external_random";
};

RandomBackendState& random_backend_state()
{
  static RandomBackendState state;
  return state;
}
} // namespace

bool gulls_random_is_stub()
{
  return random_backend_state().is_stub;
}

const char* gulls_random_backend_name()
{
  return random_backend_state().name;
}

void gulls_register_random_stub_backend(const char* backend_name)
{
  RandomBackendState& state = random_backend_state();
  state.is_stub = true;
  state.name = backend_name ? backend_name : "random_stub";
}
