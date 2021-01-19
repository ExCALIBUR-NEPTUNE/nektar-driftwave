FROM nektarpp/nektar-dev:2510247d

USER root
COPY .  /root/tmp

RUN cd /root/tmp && mkdir build && \
    cmake -DNektar++_DIR=/usr/local/lib64/nektar++/cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make install && cd /root && rm -Rf /root/tmp

USER nektar
