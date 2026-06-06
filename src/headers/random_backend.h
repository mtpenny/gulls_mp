#ifndef RANDOM_BACKEND_H
#define RANDOM_BACKEND_H

bool gulls_random_is_stub();
const char* gulls_random_backend_name();
void gulls_register_random_stub_backend(const char* backend_name);

#endif
