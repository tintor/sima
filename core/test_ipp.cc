#include <catch.hpp>

#include <stdio.h>
#include <ipp.h>
#include <arrayfire.h>

#include <core/format.h>
#include <valarray>

#define PRINT_INFO(feature, text) printf("  %-30s= ", #feature); \
    printf("%c\t%c\t", (cpuFeatures & feature) ? 'Y' : 'N', (enabledFeatures & feature) ? 'Y' : 'N'); \
    printf( #text "\n")

TEST_CASE("basic_ipp", "[ipp]") {
    const       IppLibraryVersion *libVersion;
    IppStatus   status;
    Ipp64u      cpuFeatures, enabledFeatures;

    ippInit();                      /* Initialize Intel(R) IPP library */
    libVersion = ippGetLibVersion();/* Get Intel(R) IPP library version info */
    printf("%s %s\n", libVersion->Name, libVersion->Version);

    status = ippGetCpuFeatures(&cpuFeatures, 0);/* Get CPU features and features enabled with selected library level */
    if (ippStsNoErr != status) return;
    enabledFeatures = ippGetEnabledCpuFeatures();
    printf("Features supported: by CPU\tby Intel(R) IPP\n");
    printf("------------------------------------------------\n");
    PRINT_INFO(ippCPUID_MMX,        Intel(R) Architecture MMX technology supported);
    PRINT_INFO(ippCPUID_SSE,        Intel(R) Streaming SIMD Extensions);
    PRINT_INFO(ippCPUID_SSE2,       Intel(R) Streaming SIMD Extensions 2);
    PRINT_INFO(ippCPUID_SSE3,       Intel(R) Streaming SIMD Extensions 3);
    PRINT_INFO(ippCPUID_SSSE3,      Supplemental Streaming SIMD Extensions 3);
    PRINT_INFO(ippCPUID_MOVBE,      Intel(R) MOVBE instruction);
    PRINT_INFO(ippCPUID_SSE41,      Intel(R) Streaming SIMD Extensions 4.1);
    PRINT_INFO(ippCPUID_SSE42,      Intel(R) Streaming SIMD Extensions 4.2);
    PRINT_INFO(ippCPUID_AVX,        Intel(R) Advanced Vector Extensions instruction set);
    PRINT_INFO(ippAVX_ENABLEDBYOS,  Intel(R) Advanced Vector Extensions instruction set is supported by OS);
    PRINT_INFO(ippCPUID_AES,        Intel(R) AES New Instructions);
    PRINT_INFO(ippCPUID_CLMUL,      Intel(R) CLMUL instruction);
    PRINT_INFO(ippCPUID_RDRAND,     Intel(R) RDRAND instruction);
    PRINT_INFO(ippCPUID_F16C,       Intel(R) F16C new instructions);
    PRINT_INFO(ippCPUID_AVX2,       Intel(R) Advanced Vector Extensions 2 instruction set);
    PRINT_INFO(ippCPUID_ADCOX,      Intel(R) ADOX/ADCX new instructions);
    PRINT_INFO(ippCPUID_RDSEED,     Intel(R) RDSEED instruction);
    PRINT_INFO(ippCPUID_PREFETCHW,  Intel(R) PREFETCHW instruction);
    PRINT_INFO(ippCPUID_SHA,        Intel(R) SHA new instructions);
    PRINT_INFO(ippCPUID_AVX512F,    Intel(R) Advanced Vector Extensions 512 Foundation instruction set);
    PRINT_INFO(ippCPUID_AVX512CD,   Intel(R) Advanced Vector Extensions 512 CD instruction set);
    PRINT_INFO(ippCPUID_AVX512ER,   Intel(R) Advanced Vector Extensions 512 ER instruction set);
    PRINT_INFO(ippCPUID_AVX512PF,   Intel(R) Advanced Vector Extensions 512 PF instruction set);
    PRINT_INFO(ippCPUID_AVX512BW,   Intel(R) Advanced Vector Extensions 512 BW instruction set);
    PRINT_INFO(ippCPUID_AVX512VL,   Intel(R) Advanced Vector Extensions 512 VL instruction set);
    PRINT_INFO(ippCPUID_AVX512VBMI, Intel(R) Advanced Vector Extensions 512 Bit Manipulation instructions);
    PRINT_INFO(ippCPUID_MPX,        Intel(R) Memory Protection Extensions);
    PRINT_INFO(ippCPUID_AVX512_4FMADDPS,    Intel(R) Advanced Vector Extensions 512 DL floating-point single precision);
    PRINT_INFO(ippCPUID_AVX512_4VNNIW,      Intel(R) Advanced Vector Extensions 512 DL enhanced word variable precision);
    PRINT_INFO(ippCPUID_KNC,        Intel(R) Xeon Phi(TM) Coprocessor);
    PRINT_INFO(ippCPUID_AVX512IFMA, Intel(R) Advanced Vector Extensions 512 IFMA (PMADD52) instruction set);
    PRINT_INFO(ippAVX512_ENABLEDBYOS,       Intel(R) Advanced Vector Extensions 512 is supported by OS);
}

float add(const float* a, const float* b, int n) {
    float s = 0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

float add(std::valarray<float>& a, std::valarray<float> b) {
    return (a * b).sum();
}

float add(af::array a, const af::array b) {
    return af::sum(a * b).scalar<float>();
}

float add_ipp(const float* a, const float* b, int n) {
    float m;
    ippsDotProd_32f(b, a, n, &m);
    return m;
}

template<typename Fn>
void Bench(string_view name, Fn fn) {
    print("%s", name);
    float result;
    for (int i = 0; i < 5; i++) {
        print(" %h", Duration([&]() { result = fn(); }));
    }
    println(" [%s]", result);
}

TEST_CASE("vector benchmark", "[vb][.]") {
    for (int n = 1; n <= 1000000000lu; n *= 10) {
        for (int i = 0; i < 1; i++) {
            println();
            println("n %s", n);

            std::vector<float> a(n, 2), b(n, 3);
            Bench("base", [&]() { return add(a.data(), b.data(), n); });

            std::valarray<float> va(2, n), vb(3, n);
            Bench("valarray", [&]() { return add(va, vb); });

            if (n < 1000000000 && (af::getAvailableBackends() & 1) == 1) {
                af::setBackend(AF_BACKEND_CPU);
                af::array fa = af::constant(2.0f, n), fb = af::constant(3.0f, n);
                fa.eval();
                fb.eval();
                Bench("af::array cpu", [&]() { return add(fa, fb); });
            }

            if (n < 1000000000 && (af::getAvailableBackends() & 2) == 2) {
                af::setBackend(AF_BACKEND_CUDA);
                af::array fa = af::constant(2.0f, n), fb = af::constant(3.0f, n);
                fa.eval();
                fb.eval();
                Bench("af::array cuda", [&]() { return add(fa, fb); });
            }

            if (n < 1000000000 && (af::getAvailableBackends() & 4) == 4) {
                af::setBackend(AF_BACKEND_OPENCL);
                af::array fa = af::constant(2.0f, n), fb = af::constant(3.0f, n);
                fa.eval();
                fb.eval();
                Bench("af::array opencl", [&]() { return add(fa, fb); });
            }

            std::vector<float> ia(n, 2), ib(n, 3);
            Bench("ipp", [&]() { return add_ipp(ia.data(), ib.data(), n); });
        }
    }
}
