
FROM	ubuntu:12.04

RUN		apt-get update
RUN		apt-get install -y wget libxml-simple-perl curl

RUN		cd /opt && wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.5a/genetorrent-download_3.8.5a-ubuntu2.94-12.04_amd64.deb
RUN		cd /opt && wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.5a/genetorrent-common_3.8.5a-ubuntu2.94-12.04_amd64.deb
RUN		cd /opt && wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.5a/genetorrent-upload_3.8.5a-ubuntu2.94-12.04_amd64.deb
RUN		apt-get install -y libcurl3 libxqilla6 python
RUN		apt-get install -y libboost-program-options1.48.0 libboost-system1.48.0  libboost-filesystem1.48.0 libboost-regex1.48.0

RUN		cd /opt && dpkg --install genetorrent-download_3.8.5a-ubuntu2.94-12.04_amd64.deb genetorrent-common_3.8.5a-ubuntu2.94-12.04_amd64.deb genetorrent-upload_3.8.5a-ubuntu2.94-12.04_amd64.deb
