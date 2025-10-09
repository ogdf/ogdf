ARG compiler=gcc
ARG version=13
FROM $compiler:$version

ARG ld_lib_path=
ENV LD_LIBRARY_PATH=${ld_lib_path}

RUN apt-get -y update \
 && apt-get -y install graphviz unzip wget git time ccache cmake doxygen python3-pip python3-venv \
 && apt-get -y upgrade ca-certificates \
 && mkdir -p ~/.config/pip/ \
 && printf "[global]\nbreak-system-packages = true" >> ~/.config/pip/pip.conf \
 && pip3 install pipx \
 && pipx ensurepath

# CGAL
ARG CGAL_INSTALL=false
RUN [ "$CGAL_INSTALL" = "false" ] || { \
  apt-get -y install libgmp-dev libmpfr-dev libboost-thread-dev && \
  wget -q -O CGAL.zip https://github.com/CGAL/cgal/releases/download/v$CGAL_INSTALL/CGAL-$CGAL_INSTALL-library.zip && \
  unzip -q CGAL.zip && \
  rm CGAL.zip && \
  mkdir CGAL-$CGAL_INSTALL/build && \
  cd CGAL-$CGAL_INSTALL/build && \
  cmake -DCMAKE_BUILD_TYPE=Release .. && \
  make --silent install ;}
