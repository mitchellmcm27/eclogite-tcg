# Copyright (C) 2018 ENKI Project
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of Thermocodegen
#
# Thermocodegen is free software: you can redistribute it and/or modify
# it under the terms of the MIT License
#
# Thermocodegen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# License for more details.
#
# You should have received a copy of the MIT license
# along with Thermocodegen
#
# Builds a Docker image with a full installation of the TCG SLB Database

ARG BASEIMAGE=registry.gitlab.com/cianwilson/thermocodegen:jammy-merge-upgrades

FROM ${BASEIMAGE} as db-onbuild
MAINTAINER cwilson <cwilson@carnegiescience.edu>

ARG FLAGS=""
ARG DIR="tcg_slb_database/"

USER root
ONBUILD WORKDIR /src

# install database
ONBUILD ENV DATABASE_HOME /usr/local/tcg_slb_database
# This leaves the sources inside the container. This is a limitation of
# Docker. There is some trickery in the intermediate and db
# containers that can be used to remove this source if needed, see
# below.
# Note that we copy the source from outside the file because this is a private
# repo.
ONBUILD ADD $DIR /src/tcg_slb_database
ONBUILD RUN \
    cd tcg_slb_database/database && \
    for f in `find reactions -type f -name '*.rxml'`; do FILENAMES="${FILENAMES} $(realpath $f)"; done && \
    for f in $FILENAMES; do tcg_buildrx $f -i $DATABASE_HOME $FLAGS; done

ONBUILD WORKDIR /root

########################################

FROM db-onbuild as intermediate

########################################

FROM ${BASEIMAGE} as db

# This layer manually copies the build artifacts from intermediate into
# dev-env to make the final image. This is a workaround for a well known
# limitation of Docker that you cannot cleanup after an ADD operation.
# This reduces the container size as the /src folder no
# longer exists in the final image.
COPY --from=intermediate /usr/local /usr/local

# Real by default.
# Note that because we inherit from $BASEIMAGE we do not inherit these ENV from
# db-onbuild so this must be repeated here.
ENV DATABASE_HOME /usr/local/tcg_slb_database

