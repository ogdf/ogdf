FROM docker.io/ogdf/clang:18
ARG gitref=master
RUN apt-get update && apt-get install -y libunwind-dev ninja-build
RUN git clone https://github.com/ogdf/ogdf.git \
 && cd ogdf \
 && git checkout $gitref \
 && mkdir build \
 && cd build \
 && cmake .. \
        -G "Ninja Multi-Config" \
        -DCMAKE_CONFIGURATION_TYPES="Debug;Release" \
        -DBUILD_SHARED_LIBS=ON \
        -DOGDF_MEMORY_MANAGER=POOL_NTS \
        -DOGDF_INCLUDE_CGAL=ON \
        -DOGDF_USE_ASSERT_EXCEPTIONS=ON \
        -DOGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE=ON_LIBUNWIND \
        -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE \
        -DCMAKE_POLICY_DEFAULT_CMP0069=NEW \
        -DOGDF_ARCH="x86-64"  \
        -DCMAKE_C_COMPILER_LAUNCHER=ccache \
        -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
 && cmake --build . --config Debug --target install -j6 \
 && cmake --build . --config Release --target install -j6
RUN pip install --break-system-packages ogdf-python
