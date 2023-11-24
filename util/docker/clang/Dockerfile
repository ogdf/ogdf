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
      curl \
      python3-pip \
      git \
      gcovr \
 && wget https://apt.llvm.org/llvm.sh \
 && chmod +x llvm.sh \
 && ./llvm.sh $llvmver all \
 && ln -sf /usr/bin/clang++-$llvmver /usr/bin/clang++ \
 && ln -sf /usr/bin/clang-cpp-$llvmver /usr/bin/clang-cpp \
 && ln -sf /usr/bin/clang-format-$llvmver /usr/bin/clang-format \
 && ln -sf /usr/bin/clang-tidy-$llvmver /usr/bin/clang-tidy

# Bypass pip externally managed warning to install sonar output converter,
# see https://stackoverflow.com/questions/75608323/
RUN mkdir -p ~/.config/pip/ \
 && printf "[global]\nbreak-system-packages = true" > ~/.config/pip/pip.conf
RUN pip install git+https://github.com/N-Coder/clang-tidy-converter.git

# Enable clang-tidy-cache and gcov
RUN wget https://raw.githubusercontent.com/matus-chochlik/ctcache/main/clang-tidy-cache \
 && echo -e "#!""/bin/bash\n$(realpath clang-tidy-cache) $(which clang-tidy) \"\$@\"" >> clang-tidy \
 && echo -e "#!""/bin/bash\nllvm-cov-15 gcov \"\$@\"" >> gcov \
 && chmod +x clang-tidy-cache clang-tidy gcov
