FROM nektarpp/nektar-dev:v5.5.0

ARG INSTALL_PREFIX=/usr/local

USER root
COPY .  /root/tmp
COPY example $INSTALL_PREFIX/share/doc/nektar++/driftwave-solver

RUN cd /root/tmp && mkdir build && cd build && \
    cmake -DNektar++_DIR=$INSTALL_PREFIX/lib64/nektar++/cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX/bin .. && \
    make install && \
    cd /root && chmod -R u+w /root/tmp/.git/objects && rm -Rf /root/tmp

USER nektar
