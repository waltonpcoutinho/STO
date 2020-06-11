#makefile

SYSTEM = x86-64_linux
LIBFORMAT = static_pic 

####diretorios com as libs do cplex
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio127/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio127/concert
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

####diretorios com as libs do ADOL-C
ADOLLIBDIR    = /home/walton/adolc_base/lib64/

####diretorios com as libs do AMPL
AMPLLIBDIR    = /opt/AMPL/amplapi/lib/

####diretorio com as libs do CMA-ES
#CMAESLIBDIR   = /usr/local/lib/

#### define o compilador
CPPC = g++
############################

#### opcoes de compilacao e includes
CCOPT = $(BITS_OPTION) -O0 -g -fPIC -fexceptions -DNDEBUG -DIL_STD -std=gnu++11 -pthread#-std=c++0x #-03 for optimization
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
ADOLINC       = /home/walton/adolc_base/include/
AMPLINC       = /opt/AMPL/amplapi/include/
CMAESINC      = /usr/local/include/libcmaes/
EIGENINC      = /usr/include/eigen3/
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(ADOLINC) -I$(AMPLINC) -I$(EIGENINC) -I$(CMAESINC)
#############################

#### flags do linker
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -L$(ADOLLIBDIR) -ladolc -L$(AMPLLIBDIR) -lampl -fopenmp -lm -lpthread -pthread #-L$(CMAESLIBDIR) -lcmaes -lm -lpthread
############################

#### diretorios com os source files e com os objs files
SRCDIR = src
OBJDIR = obj
SRCPOPDIR = src/ILS
OBJPOPDIR = obj/ILS
#############################

#### lista de todos os srcs e todos os objs
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
SRCSPOP = $(wildcard $(SRCPOPDIR)/*.cpp)
OBJSPOP = $(patsubst $(SRCPOPDIR)/%.cpp, $(OBJPOPDIR)/%.o, $(SRCSPOP))
#############################

#### regra principal, gera o executavel
exeGTO: $(OBJS) $(OBJSPOP)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(OBJSPOP) -o $@ $(CCLNFLAGS)
############################

#inclui os arquivos de dependencias
-include $(OBJS:.o=.d)
-include $(OBJSPOP:.o=.d)

#regra para cada arquivo objeto: compila e gera o arquivo de dependencias do arquivo objeto
#cada arquivo objeto depende do .c e dos headers (informacao dos header esta no arquivo de dependencias gerado pelo compiler)
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC) -std=c++0x  -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

$(OBJPOPDIR)/%.o: $(SRCPOPDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC) -std=c++0x  -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

#delete objetos e arquivos de dependencia
clean:
	@echo "\033[33mcleaning obj directory and bin \033[0m"
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d
	@rm -f $(OBJPOPDIR)/*.o $(OBJPOPDIR)/*.d
	@rm -f exeGTO


rebuild: clean exeGTO

