FROM registry.gitlab.com/enki-portal/thermocodegen:tf-focal

# Install Julia
RUN pip install jill -U
ENV PATH="/usr/local/bin:${PATH}"
RUN jill install 1.9 --upstream Official --confirm
RUN julia -e 'using Pkg; Pkg.add(["StatGeochem","Plots","JSON"]);'

RUN pip install --upgrade matplotlib julia

# Install Perple_X
# The Perple_X Julia interface requires all files to be in a specific directory
RUN mkdir ~/resources
RUN git clone -n https://github.com/jadconnolly/Perple_X.git ~/resources/perplex-stable \
    && cd ~/resources/perplex-stable \
    && git checkout 1aeec2f4f5d31762ecc8a5abbcb6338046406306 \
    && cd src \
    && make -j${nproc}

RUN cd ~/resources/perplex-stable\
    && cp src/* .\
    && cp -r datafiles/* .\
    && cp optionfiles/* .

# R Stuff
#RUN apt-get install -y r-base
#RUN R -e "install.packages(c('tidyverse', 'ggthemes', 'colorspace', 'latex2exp', 'latex2exp'))"

# Build eclogite reactions
RUN cd tcg_slb\
    && rm --rf database/reactions\
    && scripts/generate_reactions_eclogite -v 21\
    && scripts/build_reactions