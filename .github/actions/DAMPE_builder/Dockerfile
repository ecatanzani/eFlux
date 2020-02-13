FROM tihonav/dampesw:crmc---geant4_10_5_p01---DmpSoftware-6-0-10

#### Add custom scripts
# Set DAMPE ENVS
ADD setDAMPE.sh /usr/sbin

# Add info to the bashrc
RUN echo ". /usr/sbin/setDAMPE.sh" >> /root/.bashrc

COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]