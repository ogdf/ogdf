ARG distro=debian
ARG release=latest
FROM $distro:$release
ARG llvmver=15
# only use llvm.sh if the right clang is not already in the distro repos
# sometimes we also need to run llvm.sh twice to make it work
RUN apt-get update \
 && apt-get install -y \
      software-properties-common \
      gnupg \
      apt-transport-https \
      ca-certificates \
      wget \
      curl \
      python3-pip \
      git \
      gcovr \
 && { apt-cache search clang-$llvmver | grep clang-$llvmver; } || { \
    wget https://apt.llvm.org/llvm.sh \
 && chmod +x llvm.sh \
 && ./llvm.sh $llvmver || ./llvm.sh $llvmver; } \
 && apt-get install -y clang-${llvmver} lldb-${llvmver} lld-${llvmver} clangd-${llvmver} clang-tidy-${llvmver} \
    clang-format-${llvmver} clang-tools-${llvmver} llvm-${llvmver}-dev llvm-${llvmver}-tools libomp-${llvmver}-dev \
    libc++-${llvmver}-dev libc++abi-${llvmver}-dev libclang-common-${llvmver}-dev libclang-${llvmver}-dev \
    libclang-cpp${llvmver}-dev libunwind-${llvmver}-dev \
 && ln -sf /usr/bin/clang++-$llvmver /usr/bin/clang++ \
 && ln -sf /usr/bin/clang-cpp-$llvmver /usr/bin/clang-cpp \
 && ln -sf /usr/bin/clang-format-$llvmver /usr/bin/clang-format \
 && ln -sf /usr/bin/clang-tidy-$llvmver /usr/bin/clang-tidy \
 && ln -sf /usr/bin/llvm-profdata-${llvmver} /usr/bin/llvm-profdata \
 && ln -sf /usr/bin/llvm-cov-$llvmver /usr/bin/llvm-cov

# Bypass pip externally managed warning to install sonar output converter,
# see https://stackoverflow.com/questions/75608323/
RUN mkdir -p ~/.config/pip/ \
 && printf "[global]\nbreak-system-packages = true" > ~/.config/pip/pip.conf
RUN pip install git+https://github.com/N-Coder/clang-tidy-converter.git

# Install clang-tidy-cache
RUN wget https://raw.githubusercontent.com/matus-chochlik/ctcache/main/clang-tidy-cache \
 && mv clang-tidy-cache /usr/bin/clang-tidy-cache \
 && chmod +x /usr/bin/clang-tidy-cache
