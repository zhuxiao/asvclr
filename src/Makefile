ABPOA_PREFIX=./abPOA
CXXFLAGS = -O0 -g -Wall -fmessage-length=0 -std=c++11 -I$(ABPOA_PREFIX)/include

OBJS = asvclr_main.o Paras.o Chrome.o Genome.o Block.o Base.o \
       events.o Region.o misAlnReg.o LocalAssembly.o \
       alnDataLoader.o RefSeqLoader.o FastaSeqLoader.o \
       clipAlnDataLoader.o varCand.o covLoader.o clipReg.o \
       blatAlnTra.o Thread.o util.o meminfo.o sv_sort.o \
       genotyping.o

LIBS +=-L$(ABPOA_PREFIX)/lib -lhts -lpthread -labpoa -lz

TARGET = asvclr

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

clean:
	rm -f $(OBJS)
	
clean-all: clean
	rm -f $(TARGET)
