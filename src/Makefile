CXXFLAGS = -O0 -g -Wall -fmessage-length=0 -std=c++11

OBJS = cindel_lr.o Paras.o Chrome.o Genome.o Block.o Base.o \
       events.o Region.o misAlnReg.o LocalAssembly.o \
       alnDataLoader.o RefSeqLoader.o FastaSeqLoader.o \
       clipAlnDataLoader.o varCand.o covLoader.o clipReg.o \
       util.o

LIBS += -lhts -lpthread

TARGET = cindel_lr

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

clean:
	rm -f $(OBJS)
	
clean-all: clean
	rm -f $(TARGET)
