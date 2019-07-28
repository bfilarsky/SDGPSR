SRCDIR = ../src
OBJDIR = obj
src = $(wildcard $(SRCDIR)/*.cpp) \
	  $(wildcard *.cpp)
obj = $(addprefix $(OBJDIR)/,$(notdir $(src:.cpp=.o))) $(OBJDIR)/exampleMain.o

CXXFLAGS = -O2 -I ..
LDFLAGS = -lfftw3 -lpthread

SDGPSR : $(obj) | sampleData/iqData.dat
	$(CXX) -o $@ $^ $(LDFLAGS)
	
sampleData/iqData.dat :
	cat sampleData/iqData0.dat sampleData/iqData1.dat sampleData/iqData2.dat sampleData/iqData3.dat sampleData/iqData4.dat sampleData/iqData5.dat sampleData/iqData6.dat sampleData/iqData7.dat > sampleData/iqData.dat

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS)
	
$(OBJDIR)/%.o : %.cpp | $(OBJDIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJDIR) :
	mkdir -p $(OBJDIR)
	
.PHONY: clean
clean:
	rm -fr $(OBJDIR) SDGPSR *.bin
