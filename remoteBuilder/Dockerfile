FROM tihonav/dampesw:crmc---geant4_10_5_p01---DmpSoftware-6-0-10

#### Add custom scripts
# Set DAMPE ENVS
ADD setDAMPE.sh /usr/sbin
ADD buildProject.sh /usr/sbin

# Add info to the bashrc
RUN echo ". /usr/sbin/setDAMPE.sh" >> /root/.bashrc

#### Remove old git and install the latest version
RUN yum remove -y git \
    && mkdir /tmp \
    && cd /tmp \
    && wget https://mirrors.edge.kernel.org/pub/software/scm/git/git-2.24.0.tar.gz \
    && tar xf git-2.24.0.tar.gz \
    && cd git-2.24.0 \
    && make configure \
    && ./configure --prefix=/usr/local \
    && make all \
    && make install
#Cleaning git installation
RUN rm -rf /tmp/*

RUN export PATH=/usr/local:$PATH

CMD ["/usr/sbin/buildProject.sh"]
