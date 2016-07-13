include ../ConfigInfo.txt
include ../MakeConfig

RECO := backprojection

INCLUDE_FLAGS := $(IMOD_INCLUDE_FLAGS) $(TXBR_INCLUDE_FLAGS)

all: $(RECO)

$(RECO): run_backprojection.c backprojectionGPUsupport.o backprojection.2.o calculateBox.o backprojectionGPUsupport.o util/txbrutil.o util/txbrtext.o 
	$(CC) $(INCLUDE_FLAGS) run_backprojection.c backprojectionGPUsupport.o backprojection.2.o calculateBox.o txbrutil.o txbrtext.o -o run_backprojection $(IMOD_LINK_FLAGS) $(TIFF_LINK_FLAGS) $(CUDA_LIB_FLAGS)

backprojection.2.o: backprojection.2.c
	$(CC) $(INCLUDE_FLAGS) -c backprojection.2.c

backprojectionGPUsupport.o: backprojectionGPUsupport.cu
	$(NVCC) -shared --compiler-options '-fPIC' -c backprojectionGPUsupport.cu $(CUDA_INCLUDE_FLAGS) $(CUDA_LIB_FLAGS) $(INCLUDE_FLAGS)
	
calculateBox.o: calculateBox.c 
	$(CC) $(INCLUDE_FLAGS) -c calculateBox.c

util/txbrutil.o: $(TXBR)/txbr/util/util_C/txbrutil.c
	$(CC) $(INCLUDE_FLAGS) -c $(TXBR)/txbr/util/util_C/txbrutil.c
	
util/txbrtext.o: $(TXBR)/txbr/txbr/utilities/util_C/txbrtext.c
	$(CC) $(INCLUDE_FLAGS) -c $(TXBR)/txbr/txbr/utilities/util_C/txbrtext.c

clean:
	find . -type f -name "*.o"  -exec rm "{}" \;
	find . -type f -name "*.so" -exec rm "{}" \;
	find . -type f -name "*.linkinfo" -exec rm "{}" \;
	find . -type f -name "*~"  -exec rm "{}" \;
	find . -type f -name "core*"  -exec rm "{}" \;
	find . -type f -name "run_$(RECO)" -exec rm "{}" \;
	find . -type f -name "deviceQuery" -exec rm "{}" \;



