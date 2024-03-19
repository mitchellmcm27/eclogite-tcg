FROM registry.gitlab.com/enki-portal/thermocodegen:focal

# Install Julia and the JSON package
RUN pip install jill -U
ENV PATH="/usr/local/bin:${PATH}"
RUN jill install 1.9 --upstream Official --confirm
RUN julia -e 'using Pkg; Pkg.add(["JSON"]);'
RUN pip install --upgrade matplotlib julia

# The Perple_X Julia interface requires all files to be in a specific directory
RUN mkdir ~/resources

# Install Perple_x v7.0.10
RUN git clone -n https://github.com/jadconnolly/Perple_X.git ~/resources/perplex-stable \
    && cd ~/resources/perplex-stable \
    && git checkout 1aeec2f4f5d31762ecc8a5abbcb6338046406306 \
    && cd src \
    && make -j${nproc}

# Copy files
RUN cd ~/resources/perplex-stable\
    && cp src/* .\
    && cp -r datafiles/* .\
    && cp optionfiles/* .

# R Stuff
#RUN apt-get install -y r-base
#RUN R -e "install.packages(c('tidyverse', 'ggthemes', 'colorspace', 'latex2exp', 'latex2exp'))"

# Clone eclogitization repo
RUN cd shared && git clone https://gitlab.com/mitchellmcm27/eclogite-tcg.git
RUN chmod +x shared/eclogite-tcg/tcg_slb/scripts/generate_reactions_eclogite
RUN chmod +x shared/eclogite-tcg/tcg_slb/scripts/build_reactions

# Build eclogite reactions
RUN cd shared/eclogite-tcg/tcg_slb\
    && rm -rf database/reactions/*.rxml\
    && scripts/generate_reactions_eclogite -v 21\
    && scripts/build_reactions