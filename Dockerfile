FROM ubuntu


RUN     apt-get update
RUN     apt-get install -y curl wget g++ make python libboost-dev libboost-thread-dev libboost-system-dev zlib1g-dev ncurses-dev
WORKDIR /opt

RUN     curl -O http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.12.tar.gz
RUN     tar xvzf tophat-2.0.12.tar.gz
RUN     wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
RUN     tar xvjf samtools-0.1.19.tar.bz2
RUN     cd /opt/samtools-0.1.19 && make
RUN     cp /opt/samtools-0.1.19/*.a /usr/local/lib/
RUN     mkdir /usr/local/include/bam
RUN     cp /opt/samtools-0.1.19/*.h /usr/local/include/bam/
RUN     cd /opt/tophat-2.0.12 && ./configure --prefix=/opt --with-boost-libdir=/usr/lib/x86_64-linux-gnu/ && make && make install

RUN     wget https://github.com/alexdobin/STAR/archive/STAR_2.4.0f1.tar.gz
RUN     tar xvzf STAR_2.4.0f1.tar.gz
RUN     cd /opt/STAR-STAR_2.4.0f1 && make
RUN     cp bin/Linux_x86_64/STAR /usr/local/bin/
