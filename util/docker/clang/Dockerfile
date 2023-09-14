ARG distro=debian
ARG release=latest
FROM $distro:$release
ARG llvmver=15
RUN apt-get update \
 && apt-get install -y \
      software-properties-common \
      gnupg \
      apt-transport-https \
      ca-certificates \
      wget \
 && wget https://apt.llvm.org/llvm.sh \
 && chmod +x llvm.sh \
 && ./llvm.sh $llvmver all \
 && ln -sf /usr/bin/clang++-$llvmver /usr/bin/clang++ \
 && ln -sf /usr/bin/clang-cpp-$llvmver /usr/bin/clang-cpp \
 && ln -sf /usr/bin/clang-format-$llvmver /usr/bin/clang-format \
 && ln -sf /usr/bin/clang-tidy-$llvmver /usr/bin/clang-tidy
