ARG PIXI_IMAGE=ghcr.io/prefix-dev/pixi:0.71.1@sha256:cc1401e1553164ad032908e1f6a1dbbaa8c8e067ce4347ae831c43a066bbc66c
ARG RUNTIME_IMAGE=docker.io/rockylinux/rockylinux:8.10@sha256:f5529992e67440c1a4ae7788244d4381c6909159a88eacd95b7523ae47ced82e

FROM ${PIXI_IMAGE} AS pixi-bin

FROM ${RUNTIME_IMAGE} AS stack-builder

ARG TARGETARCH

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN test -z "${TARGETARCH}" || test "${TARGETARCH}" = "amd64"

RUN dnf install -y \
        bash \
        binutils \
        bzip2 \
        ca-certificates \
        curl \
        diffutils \
        file \
        findutils \
        gzip \
        make \
        patch \
        perl \
        tar \
        xz \
    && dnf clean all \
    && rm -rf /var/cache/dnf

COPY --from=pixi-bin /usr/local/bin/pixi /usr/local/bin/pixi

WORKDIR /opt/cppstack

COPY pixi.toml pixi.lock /opt/cppstack/

RUN pixi install --locked

ENV CPP_STACK_PREFIX=/opt/cppstack/.pixi/envs/default
ENV PATH=/opt/cppstack/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LD_LIBRARY_PATH=/opt/cppstack/.pixi/envs/default/lib

RUN test "$(x86_64-conda-linux-gnu-g++ -dumpfullversion)" = "12.3.0" \
    && pkg-config --exact-version=4.6.4 cfitsio \
    && pkg-config --exact-version=3.3.11 fftw3

FROM stack-builder AS pmi2-builder

ARG BUILD_JOBS=4
ARG SLURM_VERSION=25.11.2
ARG SLURM_URL=https://download.schedmd.com/slurm/slurm-25.11.2.tar.bz2
ARG SLURM_SHA256=9bfd844f746c1268d6e0cb14bda586b548b1c4b12f8c72181397d0cec9ce5d39

ENV PMI2_PREFIX=/opt/pmi2

RUN test "${SLURM_VERSION}" = "25.11.2" \
    && mkdir -p /build/slurm \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${SLURM_URL}" --output /build/slurm/slurm.tar.bz2 \
    && printf '%s  %s\n' \
        "${SLURM_SHA256}" /build/slurm/slurm.tar.bz2 \
        | sha256sum --check -

RUN tar -xjf /build/slurm/slurm.tar.bz2 -C /build/slurm \
    && mkdir -p /build/slurm/build \
    && cd /build/slurm/build \
    && CC=x86_64-conda-linux-gnu-gcc \
       CXX=x86_64-conda-linux-gnu-g++ \
       AR=x86_64-conda-linux-gnu-ar \
       NM=x86_64-conda-linux-gnu-nm \
       OBJDUMP=x86_64-conda-linux-gnu-objdump \
       RANLIB=x86_64-conda-linux-gnu-ranlib \
       STRIP=x86_64-conda-linux-gnu-strip \
       "/build/slurm/slurm-${SLURM_VERSION}/configure" \
        --prefix="${PMI2_PREFIX}" \
        --without-munge \
        --without-shared-libslurm \
        --disable-pam \
        --with-hdf5=no \
    && make -C contribs/pmi2 -j"${BUILD_JOBS}" \
    && make -C contribs/pmi2 install \
    && test -r "${PMI2_PREFIX}/include/slurm/pmi2.h" \
    && test -r "${PMI2_PREFIX}/lib/libpmi2.so"

FROM pmi2-builder AS openmpi-builder

ARG BUILD_JOBS=4
ARG OPENMPI_VERSION=4.1.8
ARG OPENMPI_URL=https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.8.tar.bz2
ARG OPENMPI_SHA256=466f68e3132a1dc02710cc2011fafced8336d98359fa2dae4dddcfd5719f12a9

ENV OPENMPI_PREFIX=/opt/openmpi-4.1.8

RUN test "${OPENMPI_VERSION}" = "4.1.8" \
    && mkdir -p /build/openmpi \
    && curl --fail --location --retry 5 --connect-timeout 30 \
        "${OPENMPI_URL}" --output /build/openmpi/openmpi.tar.bz2 \
    && printf '%s  %s\n' \
        "${OPENMPI_SHA256}" /build/openmpi/openmpi.tar.bz2 \
        | sha256sum --check -

RUN tar -xjf /build/openmpi/openmpi.tar.bz2 -C /build/openmpi \
    && mkdir -p /build/openmpi/build \
    && cd /build/openmpi/build \
    && CC=x86_64-conda-linux-gnu-gcc \
       CXX=x86_64-conda-linux-gnu-g++ \
       AR=x86_64-conda-linux-gnu-ar \
       NM=x86_64-conda-linux-gnu-nm \
       OBJDUMP=x86_64-conda-linux-gnu-objdump \
       RANLIB=x86_64-conda-linux-gnu-ranlib \
       STRIP=x86_64-conda-linux-gnu-strip \
       CFLAGS="-O2 -fPIC -mtune=generic" \
       CXXFLAGS="-O2 -fPIC -mtune=generic" \
       LDFLAGS="-Wl,-rpath,${PMI2_PREFIX}/lib" \
       "/build/openmpi/openmpi-${OPENMPI_VERSION}/configure" \
        --prefix="${OPENMPI_PREFIX}" \
        --with-pmi="${PMI2_PREFIX}" \
        --with-slurm \
        --disable-mpi-fortran \
        --disable-oshmem \
        --disable-dlopen \
        --without-cuda \
        --without-libfabric \
        --without-psm \
        --without-psm2 \
        --without-ucx \
        --without-verbs \
        --with-wrapper-ldflags="-Wl,-rpath,${OPENMPI_PREFIX}/lib -Wl,-rpath,${PMI2_PREFIX}/lib" \
    && make -j"${BUILD_JOBS}" \
    && make install \
    && "${OPENMPI_PREFIX}/bin/mpicxx" --showme:version \
        | grep -F "Open MPI 4.1.8" \
    && "${OPENMPI_PREFIX}/bin/ompi_info" --all \
        | grep -F "MCA ess: pmi" \
    && "${OPENMPI_PREFIX}/bin/ompi_info" --all \
        | grep -F "MCA schizo: slurm"

FROM ${RUNTIME_IMAGE} AS runtime

ARG CPP_UID=1000
ARG CPP_GID=1000

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN dnf install -y \
        bash \
        binutils \
        ca-certificates \
        file \
        findutils \
        hostname \
        procps-ng \
        shadow-utils \
    && dnf clean all \
    && rm -rf /var/cache/dnf \
    && if getent group "${CPP_GID}" >/dev/null; then \
           existing_group="$(getent group "${CPP_GID}" | cut -d: -f1)"; \
           groupmod --new-name cpp "${existing_group}"; \
       else \
           groupadd --gid "${CPP_GID}" cpp; \
       fi \
    && if getent passwd "${CPP_UID}" >/dev/null; then \
           existing_user="$(getent passwd "${CPP_UID}" | cut -d: -f1)"; \
           usermod --login cpp --home /home/cpp --move-home \
               --gid "${CPP_GID}" --shell /bin/bash "${existing_user}"; \
       else \
           useradd --uid "${CPP_UID}" --gid "${CPP_GID}" \
               --create-home --shell /bin/bash cpp; \
       fi \
    && mkdir -p \
        /workspace/src_pipe \
        /data/archive/science \
        /data/archive/dqmask \
        /data/catalogs/AstroDir \
        /data/catalogs/ExtSrcDir \
        /data/catalogs/ExtCatInput \
        /data/calib/FlatDir \
        /data/DataProcess \
        /data/output/rearr \
        /data/output/expolist \
        /data/output/fdout \
    && chown -R cpp:cpp /workspace/src_pipe /data

COPY --from=stack-builder \
    /opt/cppstack/.pixi/envs/default \
    /opt/cppstack/.pixi/envs/default
COPY --from=pmi2-builder /opt/pmi2 /opt/pmi2
COPY --from=openmpi-builder /opt/openmpi-4.1.8 /opt/openmpi-4.1.8

ENV CPP_IMAGE_ID=gxx12.3-openmpi4.1.8-pmi2
ENV CPP_STACK_PREFIX=/opt/cppstack/.pixi/envs/default
ENV PMI2_PREFIX=/opt/pmi2
ENV OPENMPI_PREFIX=/opt/openmpi-4.1.8
ENV PATH=/opt/openmpi-4.1.8/bin:/opt/cppstack/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV CPATH=/opt/pmi2/include:/opt/cppstack/.pixi/envs/default/include:/opt/cppstack/.pixi/envs/default/include/eigen3
ENV LIBRARY_PATH=/opt/openmpi-4.1.8/lib:/opt/pmi2/lib:/opt/cppstack/.pixi/envs/default/lib
ENV LD_LIBRARY_PATH=/opt/openmpi-4.1.8/lib:/opt/pmi2/lib:/opt/cppstack/.pixi/envs/default/lib
ENV PKG_CONFIG_PATH=/opt/openmpi-4.1.8/lib/pkgconfig:/opt/cppstack/.pixi/envs/default/lib/pkgconfig:/opt/cppstack/.pixi/envs/default/share/pkgconfig

RUN ln -s "${CPP_STACK_PREFIX}/bin/x86_64-conda-linux-gnu-g++" /usr/local/bin/g++ \
    && ln -s "${CPP_STACK_PREFIX}/bin/x86_64-conda-linux-gnu-gcc" /usr/local/bin/gcc \
    && test "$(g++ -dumpfullversion)" = "12.3.0" \
    && mpicxx --showme:version | grep -F "Open MPI 4.1.8" \
    && ompi_info --all | grep -F "MCA ess: pmi" \
    && ompi_info --all | grep -F "MCA schizo: slurm" \
    && pkg-config --exact-version=4.6.4 cfitsio \
    && pkg-config --exact-version=3.3.11 fftw3 \
    && test -r "${PMI2_PREFIX}/include/slurm/pmi2.h" \
    && test -r "${PMI2_PREFIX}/lib/libpmi2.so"

LABEL org.opencontainers.image.title="Fourier_Quad cpp_Standard portable HPC runtime"
LABEL org.opencontainers.image.description="G++ 12.3.0, OpenMPI 4.1.8, Slurm PMI2 client, and the cpp_Standard scientific stack"
LABEL cpppipeline.image.id="gxx12.3-openmpi4.1.8-pmi2"
LABEL cpppipeline.gxx.version="12.3.0"
LABEL cpppipeline.openmpi.version="4.1.8"
LABEL cpppipeline.pmi2.client="slurm-25.11.2"
LABEL cpppipeline.slurm.direct-launch="true"

USER cpp
WORKDIR /workspace/src_pipe

CMD ["/bin/bash"]
