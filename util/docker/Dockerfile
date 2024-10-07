ARG compiler=gcc
ARG version=13
FROM $compiler:$version

RUN apt -y update \
 && apt -y install graphviz unzip wget git time ccache cmake doxygen python3-pip python3-venv \
 && apt -y upgrade ca-certificates \
 && pip3 install pipx --break-system-packages \
 && pipx ensurepath

# CGAL
ARG CGAL_INSTALL=false
RUN [ "$CGAL_INSTALL" = "false" ] || { \
  apt -y install libgmp-dev libmpfr-dev libboost-thread-dev && \
  wget -q -O CGAL.zip https://github.com/CGAL/cgal/releases/download/v5.5.2/CGAL-5.5.2-library.zip && \
  unzip -q CGAL.zip && \
  rm CGAL.zip && \
  mkdir CGAL-5.5.2/build && \
  cd CGAL-5.5.2/build && \
  cmake -DCMAKE_BUILD_TYPE=Release .. && \
  make --silent install ;}
