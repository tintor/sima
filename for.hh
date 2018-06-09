#ifndef __FOR_H__
#define __FOR_H__

#define FOR(I, N) for (decltype(N) I = 0; I < N; I++)
#define FOR_EACH(A, B) for (auto& A : B)
#define FOR_EACH_EDGE(A, B, V) for (auto *B = V.begin(), *A = V.begin() + 2; B < V.begin() + 3; A = B++)

#endif
